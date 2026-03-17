using Dates

const SOLAR_CONSTANT=1366.1 # W/m^2

"""
    Irradiance(dni, dhi, ghi)

Container for solar irradiance components.

# Fields
- `dni`: Direct Normal Irradiance (W/m²).
- `dhi`: Diffuse Horizontal Irradiance (W/m²).
- `ghi`: Global Horizontal Irradiance (W/m²).
"""
struct Irradiance{T<:Real}
    dni::T
    dhi::T
    ghi::T
end

"""
    day_angle(datetime::DateTime, offset::Int)

Calculates the day angle [degrees] for the Earth's orbit around the Sun. 

For the Spencer method, `offset=1`; for the ASCE method, `offset=0`.
"""
function day_angle(datetime::DateTime, offset)
    return 360.0(dayofyear(datetime) - offset) / 365.0
end

"""
    extraterrestrial_irradiance_spencer1971(datetime::DateTime)

Calculates the extraterrestrial irradiance [W/m^2], as proposed in 
[spencer1971fourier](@cite).
"""
function extraterrestrial_irradiance_spencer1971(datetime::DateTime)
    da = day_angle(datetime, 1.0)
    R = 1.00011 + 0.034221cosd(da) + 0.00128sind(da) + 0.000719cosd(2da) + 7.7e-05sind(2da)
    return R * SOLAR_CONSTANT
end

"""
    extraterrestrial_irradiance_asce(datetime::DateTime)

Calculates the extraterrestrial irradiance [W/m^2] using the ASCE evapotranspiration 
equation [walter2000asce](@cite).
"""
function extraterrestrial_irradiance_asce(datetime::DateTime)
    da = day_angle(datetime, 0.0)
    R = 1 + 0.033cosd(da) 
    return R * SOLAR_CONSTANT
end

"""
    extraterrestrial_irradiance_nrel(heliocentric_radius)

Calculates the extraterrestrial irradiance [W/m^2] as implemented by the NREL SPA 
[reda2004solar](@cite). 
"""
function extraterrestrial_irradiance_nrel(heliocentric_radius)
    return SOLAR_CONSTANT / heliocentric_radius^2
end
