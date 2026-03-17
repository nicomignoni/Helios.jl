using Pkg.Artifacts, Serialization

# The linke-turbidity.jld contains a single 2160 x 4320 x 12 array of type uint8. 
# To determine the Linke urbidity for a position on the Earth's surface 
# for a given month do the following: LT = LinkeTurbidity(LatitudeIndex, LongitudeIndex, 
# month).

# The nodes of the grid are 5' (1/12=0.0833[arcdeg]) apart.
# From Section 8 of Aerosol optical depth and Linke turbidity climatology
# http://www.meteonorm.com/images/uploads/downloads/ieashc36_report_TL_AOD_climatologies.pdf
# 1st row: 89.9583 S, 2nd row: 89.875 S
# 1st column: 179.9583 m^2, 2nd column: 179.875 m^2
const LINKE_TURBIDITY_METEOTEST = open(artifact"data/linke-turbidity.jld", "r") do io;
    deserialize(io)
end

const ATMOSPHERIC_PRESSURE = 1013.25 # mbar
const ABSOLUTE_ZERO = 273.15 # kelvin

"""
    aod_bb_hulstrom1980(aod380, aod500)

Computes the approximate broadband aerosol optical depth (AOD), as proposed in 
[bird1980direct, bird1981review](@cite).

It correlates broadband AOD using two wavelengths, 380 nm and 500 nm, denoted with 
`aod380` and `aod500`, respectively.
"""
function aod_bb_hulstrom1980(aod380, aod500)
    return 0.27583aod380 + 0.35aod500
end

"""
    linke_turbidity_meteotest(location::Location, datetime::DateTime)

Computes the Linke turbidity based on SoDa [sodapro, remund2003worldwide](@cite).

It performs a [bilinear interpolation](https://en.wikipedia.org/wiki/Bilinear_interpolation) 
of the Linke turbidity array, for the corrsponding latitude, longitude, and month values. 
In such an array, rows represent global latitudes, columns represent global 
longitudes, and depth (third dimension) represents months of the year.
"""
function linke_turbidity_meteotest(location::Location, datetime::DateTime)
    lt = interpolated_value(
            LINKE_TURBIDITY_METEOTEST, 
            location.latitude, 
            location.longitude,
            month(datetime)
        )

    # The values within the Linke turbidity array are scaled by 20
    return lt / 20
end

"""
    linke_turbidity_kasten1996(absolute_airmass, precipitable_water, aod_bb)

Computes the Linke turbidity using Kasten pyrheliometric formula.

Note that broadband aerosol optical depth (AOD) can be approximated by AOD measured at 700 
nm according to Molineaux [molineaux1998equivalence](@cite). Bird and Hulstrom offer an 
alternate approximation using broadband aerosol optical depth (`aod_bb`),
measured at 380 nm and 500 nm. Based on original implementation by Armel Oumbe. Note that 
these calculations are only valid for `absolute_airmass` less than 5 and 
`precipitable_water` less than 5 cm.
"""
function linke_turbidity_kasten1996(absolute_airmass, precipitable_water, aod_bb)
    δ_cda = -0.101 + 0.235absolute_airmass^(-0.16)
    δ_w = 0.112absolute_airmass^(-0.55) * precipitable_water^0.34
    return -(9.4 + 0.9absolute_airmass) * 
            log(exp(-absolute_airmass * (δ_cda + δ_w + aod_bb))) / absolute_airmass
end

"""
    pressure2altitude(pressure)

Computes the local altitude from site pressure. 

The following assumptions [psas2004altitudepressure](@cite) hold:
- Base pressure: 101325 Pa
- Temperature at zero altitude: 288.15 K
- Gravitational acceleration: 9.80665 m/s^2
- Lapse rate: -6.5E-3 K/m
- Gas constant for air: 287.053 J/(kg K)
- Relative Humidity: 0%
"""
function pressure2altitude(pressure)
    return 44331.5 - 4946.62(0.01pressure)^0.190263
end

"""
    altitude2pressure(altitude)

Computes the site pressure from local altitude.

It follows the same assumption of [`pressure2altitude`](@ref).
"""
function altitude2pressure(altitude)
    return ((44331.514 - altitude) / 11880.516)^(1 / 0.1902632)
end

"""
    absolute_airmass(relative_airmass, pressure)

Computes the absolute airmass [gueymard1993critical](@cite).
"""
function absolute_airmass(relative_airmass, pressure)
    return relative_airmass * pressure / ATMOSPHERIC_PRESSURE
end

# Relative airmass ------------------------------------------------------------------------
"""
    relative_airmass_simple(solpos::SolarPosition)

Computes the relative airmass as secant of the Sun's apparent zenith.
"""
function relative_airmass_simple(solpos::SolarPosition)
    return secd(apparent_zenith(solpos))
end

"""
    relative_airmass_kasten1966(solpos::SolarPosition)

Computes the relative airmass as proposed in [kasten1965new](@cite).
"""
function relative_airmass_kasten1966(solpos::SolarPosition)
    app_zenith = apparent_zenith(solpos)
    return 1.0 / (cosd(app_zenith) + 
           0.15((93.885 - app_zenith)^(- 1.253)))
end

"""
    relative_airmass_youngirvine1967(solpos::SolarPosition)

Computes the relative airmass as proposed in [young1967multicolor](@cite).
"""
function relative_airmass_youngirvine1967(solpos::SolarPosition)
    sec_zenith = 1.0 / cosd(zenith(solpos))
    return sec_zenith * (1 - 0.0012(sec_zenith^2 - 1))
end

"""
    relative_airmass_kastenyoung1989(solpos::SolarPosition)

Computes the relative airmass as proposed in [kasten1989revised](@cite). 
"""
function relative_airmass_kastenyoung1989(solpos::SolarPosition)
    app_zenith = min(96.07995, apparent_zenith(solpos))
    return 1.0 / (cosd(app_zenith) + 0.50572((96.07995 - app_zenith)^(- 1.6364)))
end

"""
    relative_airmass_gueymard1993(solpos::SolarPosition)

Computes the relative airmass as proposed in 
[gueymard1993critical, gueymard1993development](@cite).
"""
function relative_airmass_gueymard1993(solpos::SolarPosition)
    app_zenith = apparent_zenith(solpos)
    return 1.0 / (cosd(app_zenith) + 0.00176759app_zenith*
           ((94.37515 - app_zenith)^(-1.21563)))
end

"""
    relative_airmass_young1994(solpos::SolarPosition)

Computes the relative airmass as proposed in [young1994air](@cite).
"""
function relative_airmass_young1994(solpos::SolarPosition)
    _zenith = zenith(solpos)
    return (1.002432cosd(_zenith)^2 + 0.148386cosd(_zenith) + 0.0096467) /
           (cosd(_zenith)^3 + 0.149864cosd(_zenith)^2 + 
            0.0102963cosd(_zenith) + 0.000303978)
end

"""
    relative_airmass_pickering2002(solpos::SolarPosition)

Computes the relative airmass as proposed in [pickering2002southern](@cite).
"""
function relative_airmass_pickering2002(solpos::SolarPosition)
    app_zenith = apparent_zenith(solpos)
    return 1.0 / sin(90.0 - app_zenith + 244.0 / (165.0 + 47.0(90.0 - app_zenith)^1.1))
end

"""
    relative_airmass_gueymard2003(solpos::SolarPosition)

Computes the relative airmass as proposed in [gueymard2003direct, gueymard2019clear](@cite). 
"""
function relative_airmass_gueymard2003(solpos::SolarPosition)
    app_zenith = apparent_zenith(solpos)
    return 1.0 / (cosd(app_zenith) + 0.48353app_zenith^0.095846 / (96.741 - app_zenith)^1.754)
end

# Vapor -----------------------------------------------------------------------------------
"""
    apparent_vapor_scale_height_gueymard94(air_temperature)

Computes the apparent water vapor scale height [km] as proposed in 
[gueymard1994analysis; Eq. 4](@cite).
"""
function apparent_vapor_scale_height_gueymard94(air_temperature)
    T_kelvin = air_temperature + ABSOLUTE_ZERO
    T_normlized = T_kelvin / ABSOLUTE_ZERO
    return 0.4976 + 1.5265T_normalized + exp(13.6897T_normalized - 14.9188T_normalized^3)
end

"""
    surface_vapor_density_gueymard94(
        air_temperature,
        relative_humidity,
        apparent_vapor_scale_height
    )

Computes the surface water vapor density [g/m^3] as proposed in 
[gueymard1994analysis; Eq. 2](@cite).
"""
function surface_vapor_density_gueymard94(
    air_temperature, 
    relative_humidity,
    apparent_vapor_scale_height
)
    return 216.7relative_humidity * apparent_vapor_scale_height / air_temperature 
end

"""
    saturation_vapor_pressure_gueymard93(air_temperature)

Computes the saturation water vapor pressure [mbar] as proposed in 
[gueymard1993critical; Eq. 1](@cite).
"""
function saturation_vapor_pressure_gueymard93(air_temperature)
    T_kelvin = air_temperature + ABSOLUTE_ZERO
    return exp(22.330 - 49.140(100/T_kelvin) - 10.922(100/T_kelvin)^2 - 
               0.39015T_kelvin/100)
end

"""
    saturation_vapor_pressure_oke2018(dew_temperature)

Computes the saturation water vapor pressure [mbar] as proposed in [oke2018guide](@cite).
"""
function saturation_vapor_pressure_oke2018(dew_temperature)
    return 6.112exp(17.62dew_temperature / (243.12 + dew_temperature)) 
end

"""
    vapor_pressure_oke2018(air_temperature)

Computes the water vapor pressure [mbar] as proposed in [oke2018guide](@cite). 
"""
function vapor_pressure_oke2018(air_temperature)
    return 6.112exp(17.62air_temperature / (243.12 + air_temperature))
end

# Precipitable water and humidity ---------------------------------------------------------
"""
    precipitable_water_gueymard94(saturation_vapor_pressure, surface_vapor_density)

Computes precipitable water [cm] as proposed in [gueymard1994analysis; Eq. 1](@cite).

The accuracy of this method is approximately 20% for moderate precipitable water (1-3 cm) 
and less accurate otherwise. 
"""
function precipitable_water_gueymard94(
    saturation_vapor_pressure,
    surface_vapor_density
)
    return max(0.1saturation_vapor_pressure * surface_vapor_density, 0.1)
end

"""
    relative_humidity_oke2018(vapor_pressure, saturation_vapor_pressure)

Computes the relative humidity [%] as proposed in [oke2018guide](@cite). 
"""
function relative_humidity_oke2018(vapor_pressure, saturation_vapor_pressure)
    return 100saturation_vapor_pressure / vapor_pressure
end
