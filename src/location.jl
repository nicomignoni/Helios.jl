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

# Arguments
- `latitude`: latitude (degrees)
- `longitude`: longitude (degrees)
- `altitude`: altitude above sea level (meters)

Ambient pressure is calculated using [`altitude2pressure`](@ref). Default temperature is
25° (celsius).
"""
function Location(
    latitude,
    longitude,
    altitude;
    temperature = 25.0, 
    pressure = altitude2pressure(altitude)
)
    latitude, longitude, altitude, temperature, pressure =
        promote(latitude, longitude, altitude, temperature, pressure)

    return Location(latitude, longitude, altitude, temperature, pressure)
end
