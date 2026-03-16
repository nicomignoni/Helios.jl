using Pkg.Artifacts, Serialization

const ALTITUDE_STEP = 28.0
const ALTITUDE_MIN = 450.0
const ALTITUDE_DATA = open(artifact"data/altitude.jld", "r") do io;
    deserialize(io)
end

"""
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
