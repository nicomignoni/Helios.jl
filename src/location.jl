using Pkg.Artifacts, Serialization

const ALTITUDE_STEP = 28.0
const ALTITUDE_MIN = 450.0
const ALTITUDE_DATA = open(artifact"data/altitude.jld", "r") do io;
    deserialize(io)
end

"""
    ALBEDO

Ground surface albedo, also known as the reflection coefficient, represents the degree
of reflectivity of the surrounding enviroment [payne1972albedo](@cite).
"""
const ALBEDO = (
    sea = 0.06,
    urban = 0.18,
    asphalt = 0.12,
    concrete = 0.30,
    soil = 0.17,
    grass = 0.20,
    sand = 0.40,
    snow = 0.65,
    aluminum = 0.85,
    copper = 0.74,
    fresh_grass = 0.26,
    fresh_snow = 0.75,
    fresh_steel = 0.35,
    dirty_steel = 0.08,
)

"""
    Location(latitude, longitude, altitude, temperature, pressure)

Represents a physical geographic location together with ambient conditions.
"""
struct Location{T<:Real}
    latitude::T
    longitude::T
    altitude::T
    temperature::T
    pressure::T
end

"""
    Location(latitude, longitude, altitude)

Create a `Location` from geographic coordinates and ambient conditions.
"""
function Location(
    latitude,
    longitude;
    altitude = altitude(latitude, longitude),
    temperature = 25.0, 
    pressure = altitude2pressure(altitude)
)
    return Location(promote(latitude, longitude, altitude, temperature, pressure)...)
end

# TODO: embed the 0.0 for no-data directly in ALTITUDE_DATA 
"""
        altitude(latitude, longitude)

    Look up location altitude from low-resolution altitude map obtained by aggredating 
    mutliple open data sources by [Mapzen](https://www.mapzen.com). See 
    https://github.com/pvlib/pvlib-python/blob/main/pvlib/location.py
"""
function altitude(latitude, longitude)
    alt = interpolated_value(ALTITUDE_DATA, latitude, longitude)
    return alt ≈ 255.0 ? 0.0 : ALTITUDE_STEP*alt - ALTITUDE_MIN
end
