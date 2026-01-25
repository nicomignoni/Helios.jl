"""
    clearsky_ineichen(
        sun_apparent_elevation,
        observer_altitude,
        absolute_airmass,
        linke_turbidity,
        extraterrestial_radiation,
        perez_enhancement::Bool
    )

Returns the global horizontal irradiance (GHI), direct normal irradiance (DNI), and diffuse 
horizontal (DHI), all in [W/m^2], following the Ineichen/Perez clear sky model 
[ineichen2002new, perez2002new](@cite).

A report on clear sky models found the Ineichen/Perez model to have excellent performance 
with a minimal input data set [stein2012global](@cite). 

# Arguments
- `sun_apparent_elevation` - The apparent (refraction corrected) Sun's 
    elevation angle. See [`topocentric_apparent_elevation`](@ref).
- `observer_altitude` - [m] Altitude above sea level. 
- `absolute_airmass` - Pressure corrected airmass. See [`absolute_airmass`](@ref).
- `linke_turbidity` - Linke Turbidity. See 
- `extraterrestial_radiation` - [W/m^2] Extraterrestrial irradiance.
- `perez_enhancement::Bool` - Controls if the Perez enhancement factor should be applied. 
    Setting to `true` 
    [may produce spurious results](https://github.com/pvlib/pvlib-python/issues/435) 
    for times when the Sun is near the horizon and the airmass is high.
"""
function clearsky_ineichen(
    sun_apparent_elevation,
    observer_altitude,
    absolute_airmass,
    linke_turbidity,
    extraterrestial_radiation,
    perez_enhancement::Bool
)
    sin_elev = max(sind(sun_apparent_elevation), 0.0)

    tl = linke_turbidity

    fh1 = exp(-observer_altitude / 8000.0)
    fh2 = exp(-observer_altitude / 1250.0)
    cg1 = 5.09e-05observer_altitude + 0.868
    cg2 = 3.92e-05observer_altitude + 0.0387

    ghi = exp(-cg2 * absolute_airmass * (fh1 + fh2 * (tl - 1)))

    # https://github.com/pvlib/pvlib-python/issues/435
    if perez_enhancement
        ghi *= exp(0.01absolute_airmass^1.8)
    end

    # use maximum to map airmass nans to 0s. multiply and divide by tl to
    # reinsert tl nans
    ghi = cg1 * extraterrestial_radiation * sin_elev * tl / tl * max(ghi, 0.0)

    # From [1] (Following [2] leads to 0.664 + 0.16268 / fh1)
    # See https://github.com/pvlib/pvlib-python/pull/808
    b = 0.664 + 0.163 / fh1
    # BncI = "normal beam clear sky radiation"
    bnci = b * exp(-0.09absolute_airmass * (tl - 1))
    bnci = extraterrestial_radiation * max(bnci, 0.0)

    # "empirical correction" SE 73, 157 & SE 73, 312.
    bnci_2 = (1 - (0.1 - 0.2exp(-tl)) / (0.1 + 0.882 / fh1)) / sin_elev
    bnci_2 = ghi * min(max(bnci_2, 0.0), 1e20)

    dni = min(bnci, bnci_2)
    dhi = ghi - dni * sin_elev

    return dni, dhi, ghi
end

"""
    clearsky_haurwitz(sun_apparent_elevation)

Implements the Haurwitz clear sky model for global horizontal irradiance (GHI) as presented 
in [haurwitz1945insolation, haurwitz1946insolation](@cite). 

A report on clear sky models found the Haurwitz model to have the best performance in terms 
of average monthly error among models which require only the Sun's elevation 
[stein2012global](@cite).

# Arguments
- `sun_apparent_elevation` - The apparent (refraction corrected) Sun's 
    elevation angle. See [`topocentric_apparent_elevation`](@ref).
"""
function clearsky_haurwitz(sun_apparent_elevation)
    sin_elev = sind(sun_apparent_elevation)
    ghi = sin_elev < 0.0 ? 0.0 : 1098.0sin_elev * exp(-0.059/sin_elev)
    return ghi
end

"""
    clearsky_simplified_solis(
        sun_apparent_elevation,
        aod700,
        precipitable_water,
        pressure,
        extraterrestial_radiation
    )

Calculate the clear sky direct normal irradiance (DNI), and diffuse horizontal irradiance 
(DHI) according to the simplified Solis model. 

Reference [ineichen2008broadband](@cite) describes the accuracy of the model as being 15, 
20, and 18 W/m^2 for the beam, global, and diffuse components, respectively. 
Reference [ineichen2016validation](@cite) provides comparisons with other clear sky models.

# Arguments
- `sun_apparent_elevation` - The apparent (refraction corrected) Sun's 
    elevation angle. See [`topocentric_apparent_elevation`](@ref).
- `aod700` - The aerosol optical depth at 700 nm. Algorithm derived for values 
    between 0 and 0.45.
- `precipitable_water` - [cm] The precipitable water of the atmosphere. Algorithm 
    derived for values between 0.2 and 10 cm. Values less than 0.2 will be assumed to be 
    equal to 0.2.
- `pressure` - [Pa] The atmospheric pressure. Algorithm derived for altitudes between 
    sea level and 7000 m, or 101325 and 41000 Pascals.
- `extraterrestial_radiation` - [W/m^2] Extraterrestrial irradiance.
"""
function clearsky_simplified_solis(
    sun_apparent_elevation,
    aod700,
    precipitable_water,
    pressure,
    extraterrestial_radiation
)
    # algorithm fails for pw < 0.2
    w = maximum(precipitable_water, 0.2)

    log_w = log(w)
    log_scaled_p = log(pressure / ATMOSPHERIC_PRESSURE)

    I₀₀ = 1.08w^0.0051
    I₀₁ = 0.97w^0.032
    I₀₂ = 0.12w^0.56
    I₀ = extraterrestial_radiation * 
         (I₀₀*aod700^2 + I₀₁*aod700 + I₀₀ + 0.071log_scaled_p)

    tb₁ = 1.82 + 0.056log_w + 0.0071log_w^2
    tb₀ = 0.33 + 0.045log_w + 0.0096log_w^2
    tbₚ = 0.0089w + 0.13
    τb = tb₁*aod700 + tb₀ + tbₚ*log_scaled_p

    b₁ = 0.00925aod700^2 + 0.0148aod700 - 0.0172
    b₀ = -0.7565aod700^2 + 0.5057aod700 + 0.4557
    b = b₁ * log_w + b₀

    tg₁ = 1.24 + 0.047log_w + 0.0061log_w^2
    tg₀ = 0.27 + 0.043log_w + 0.0090log_w^2
    tgₚ = 0.0079w + 0.1
    τg = tg₁*aod700 + tg₀ + tgₚ*log_scaled_p

    g = -0.0147log_w - 0.3079aod700^2 + 0.2846aod700 + 0.3798

    if aod700 < 0.05
        td₄ = 86.0w - 13800.0
        td₃ = -3.11w + 79.4
        td₂ = -0.23w + 74.8
        td₁ = 0.092w - 8.86
        td₀ = 0.0042w + 3.12
        tdₚ = -0.83(1.0 + aod700)^(-17.2)
    else
        td₄ = -0.21w + 11.6
        td₃ = 0.27w - 20.7
        td₂ = -0.134w + 15.5
        td₁ = 0.0554w - 5.71
        td₀ = 0.0057w + 2.94
        tdₚ = -0.71(1.0 + aod700)^(-15.0)
    end

    τd = td₄*aod700^4 + td₃*aod700^3 + td₂*aod700^2 +
         td₁*aod700 + td₀ + tdₚ*log_scaled_p
    
    dₚ = 1 / (18.0 + 152.0aod700)
    d = -0.337aod700^2 + 0.63aod700 + 0.116 + dₚ*log_scaled_p
    
    # this prevents the creation of nans at night instead of 0s
    sin_elev = maximum(1.0e-30, sind(sun_apparent_elevation))
    
    dni = I₀ * exp(-τb / sin_elev^b)
    ghi = I₀ * exp(-τg / sin_elev^g) * sin_elev
    dhi = I₀ * exp(-τd / sin_elev^d)
    
    return dhi, ghi, dhi
end

