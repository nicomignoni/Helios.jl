using Dates

const SOLAR_CONSTANT=1366.1 # W/m^2

"""
    day_angle(datetime::DateTime, offset::Int)

Calculates the day angle [deg] for the Earth's orbit around the Sun. 

For the Spencer method, `offset=1`; for the ASCE method, `offset=0`.
"""
function day_angle(datetime::DateTime, offset::Int)
    doy = dayofyear(datetime)
    return 360.0(doy - offset) / 365.0
end

"""
    extraterrestrial_irradiance_spencer1971(datetime::DateTime)

Calculates the extraterrestrial irradiance [W/m^2], as proposed in 
[spencer1971fourier](@cite).

```math
    \\begin{align*}
        $VN_EXTRATERRESTRIAL_IRRADIANCE &= $VN_SOLAR_CONSTANT(1.00011 + 0.034221\\cos 
        $VN_DAY_ANGLE + \\\\
        & \\quad 0.00128\\sin $VN_DAY_ANGLE + 0.000719\\cos 2$VN_DAY_ANGLE + 
        7.7e-05\\sin 2$VN_DAY_ANGLE)
    \\end{align*}
```
"""
function extraterrestrial_irradiance_spencer1971(datetime::DateTime)
    da = day_angle(datetime, 1.0)
    R = 1.00011 + 0.034221cosd(da) + 0.00128sind(da) + 0.000719cosd(2da) + 7.7e-05sind(2da)
    return R * SOLAR_CONSTANT
end

"""
    extraterrestrial_irradiance_asce(day_angle, solar_constant)

Calculates the extraterrestrial irradiance using the ASCE evapotranspiration equation 
[walter2000asce](@cite).

```math
    $VN_EXTRATERRESTRIAL_IRRADIANCE = $VN_SOLAR_CONSTANT(1 + 0.033\\cos($VN_DAY_ANGLE))
```
- ``$VN_DAY_ANGLE`` corresponds to `day_angle` [deg]
- ``$VN_SOLAR_CONSTANT`` corresponds to `solar_constant` [W/m^2]
"""
function extraterrestrial_irradiance_asce(day_angle, solar_constant)
    R = 1 + 0.033cosd(day_angle) 
    return R * solar_constant
end

"""
    extraterrestrial_irradiance_nrel(heliocentric_radius, solar_constant)

Calculates the extraterrestrial irradiance as implemented by the NREL SPA 
[reda2004solar](@cite). 

```math
    $VN_EXTRATERRESTRIAL_IRRADIANCE = \\frac{$VN_SOLAR_CONSTANT}{$VN_HELIOCENTRIC_RADIUS^2}
```
- ``$VN_HELIOCENTRIC_RADIUS`` corresponds to `heliocentric_radius` [AU]
- ``$VN_SOLAR_CONSTANT`` corresponds to `solar_constant` [W/m^2]
"""
function extraterrestrial_irradiance_nrel(heliocentric_radius, solar_constant)
    return solar_constant / heliocentric_radius^2
end
