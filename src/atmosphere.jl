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

# Atmospheric pressure at sea level
const ATMOSPHERIC_PRESSURE = 1013.25 # mbar

const ABSOLUTE_ZERO = 273.15 # kelvin

# Bounds for latitude and longitude
const LATITUDE_MAX = 90.0
const LATITUDE_MIN = -90.0
const LONGITUDE_MAX = 180.0
const LONGITUDE_MIN = -180.0

const LATITUDE_RANGE = LATITUDE_MAX - LATITUDE_MIN
const LONGITUDE_RANGE = LONGITUDE_MAX - LONGITUDE_MIN

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
    month_index = month(datetime)
    i = size(LINKE_TURBIDITY_METEOTEST, 1) * 
        (location.latitude - LATITUDE_MIN) / LATITUDE_RANGE
    j = size(LINKE_TURBIDITY_METEOTEST, 2) * 
        (location.longitude - LONGITUDE_MIN) / LONGITUDE_RANGE

    i⁻, i⁺ = floor(i) |> Int, ceil(i) |> Int 
    j⁻, j⁺ = floor(j) |> Int, ceil(j) |> Int 

    # Index position between rounded row and column as [0, 1] fraction
    ν, τ = i - i⁻, j - j⁻

    # Bilinear interpolation of row-column data
    lt = ν * τ * LINKE_TURBIDITY_METEOTEST[i⁻, j⁻, month_index] + 
         ν * (1 - τ) * LINKE_TURBIDITY_METEOTEST[i⁻, j⁺, month_index] + 
         (1 - ν) * τ * LINKE_TURBIDITY_METEOTEST[i⁺, j⁻, month_index] + 
         (1 - ν) * (1 - τ) * LINKE_TURBIDITY_METEOTEST[i⁺, j⁺, month_index] 

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
```math
    $VN_ALTITUDE = 44331.5 - 4946.62\\left(\\frac{$VN_PRESSURE}{100}\\right)^{0.190263}
```
where ``$VN_ALTITUDE`` is the altitude and ``$VN_PRESSURE`` is the `pressure`.
"""
function pressure2altitude(pressure)
    return 44331.5 - 4946.62(0.01pressure)^0.190263
end

"""
    altitude2pressure(altitude)

Computes the site pressure from local altitude.

It follows the same assumption of [`pressure2altitude`](@ref) an it is calculated as 
```math
    $VN_PRESSURE = \\left(\\frac{44331.514 - $VN_ALTITUDE}{11880.516}\\right)^{
    \\frac{1}{0.1902632}}
```
where ``$VN_PRESSURE`` is the pressure and ``$VN_ALTITUDE`` is the `altitude`.
"""
function altitude2pressure(altitude)
    return ((44331.514 - altitude) / 11880.516)^(1 / 0.1902632)
end

"""
    absolute_airmass(relative_airmass, pressure)

Computes the absolute airmass [gueymard1993critical](@cite).

```math 
    $VN_ABSOLUTE_AIRMASS = $VN_RELATIVE_AIRMASS \\frac{$VN_PRESSURE}
    {$VN_ATMOSPHERIC_PRESSURE} 
```
where  ``$VN_ABSOLUTE_AIRMASS`` is the absolute airmass, ``$VN_RELATIVE_AIRMASS`` is the 
`relative_airmass`, ``$VN_PRESSURE`` is the `pressure`, and ``$VN_ATMOSPHERIC_PRESSURE`` 
is the atmospheric pressure.
"""
function absolute_airmass(relative_airmass, pressure)
    return relative_airmass * pressure / ATMOSPHERIC_PRESSURE
end

# Relative airmass ------------------------------------------------------------------------
"""
    relative_airmass_simple(apparent_zenith)
    relative_airmass_simple(s::SolarPosition)

Computes the relative airmass as secant of the Sun's apparent zenith.

```math
$VN_RELATIVE_AIRMASS = \\sec $VN_TOPOCENTRIC_APPARENT_ZENITH
```
where ``$VN_RELATIVE_AIRMASS`` is the relative airmass and 
``$VN_TOPOCENTRIC_APPARENT_ZENITH`` is the Sun's apparent zenith, i.e., the complementary 
to 90° of the [`topocentric_apparent_elevation`](@ref).
"""
function relative_airmass_simple(apparent_zenith)
    return secd(apparent_zenith)
end

relative_airmass_simple(s::SolarPosition) = relative_airmass_simple(s.apparent_zenith)

"""
    relative_airmass_kasten1966(apparent_zenith)
    relative_airmass_kasten1966(s::SolarPosition)

Computes the relative airmass as proposed in [kasten1965new](@cite).
 
```math
    $VN_RELATIVE_AIRMASS = \\frac{1}{\\cos $VN_TOPOCENTRIC_APPARENT_ZENITH + 0.15(93.885 - 
    $VN_TOPOCENTRIC_APPARENT_ZENITH)^{-1.253})}
```
where ``$VN_TOPOCENTRIC_APPARENT_ZENITH`` is the Sun's apparent_zenith zenith, i.e., 
the complementary to 90° of the [`topocentric_apparent_elevation`](@ref).
"""
function relative_airmass_kasten1966(apparent_zenith)
    return 1.0 / (cosd(apparent_zenith) + 0.15((93.885 - apparent_zenith)^(- 1.253)))
end

relative_airmass_kasten1966(s::SolarPosition) =
    relative_airmass_kasten1966(s.apparent_zenith)

"""
    relative_airmass_youngirvine1967(zenith)
    relative_airmass_youngirvine1967(s::SolarPosition)

Computes the relative airmass as proposed in [young1967multicolor](@cite).

```math
    $VN_RELATIVE_AIRMASS = \\sec $VN_TOPOCENTRIC_ZENITH (1 - 0.0012(
    \\sec $VN_TOPOCENTRIC_ZENITH^2 - 1))
```
where ``$VN_TOPOCENTRIC_ZENITH`` is the Sun's zenith, i.e., the complementary to 90° of the
[`topocentric_elevation`](@ref).
"""
function relative_airmass_youngirvine1967(zenith)
    sec_zenith = 1.0 / cosd(zenith)
    return sec_zenith * (1 - 0.0012(sec_zenith^2 - 1))
end

relative_airmass_youngirvine1967(s::SolarPosition) = 
    relative_airmass_youngirvine1967(s.zenith)

"""
    relative_arimass_kastenyoung1989(apparent_zenith)
    relative_airmass_kastenyoung1989(s::SolarPosition)

Computes the relative airmass as proposed in [kasten1989revised](@cite). 

```math
    $VN_RELATIVE_AIRMASS = \\frac{1}{\\cosd $VN_TOPOCENTRIC_APPARENT_ZENITH + 
    0.50572((96.07995 - $VN_TOPOCENTRIC_APPARENT_ZENITH)^{- 1.6364})}
```
where ``$VN_TOPOCENTRIC_APPARENT_ZENITH`` is the Sun's apparent zenith, i.e., the 
complementary to 90° of the [`topocentric_apparent_elevation`](@ref).
"""
function relative_airmass_kastenyoung1989(apparent_zenith)
    apparent_zenith = min(96.07995, apparent_zenith)
    return 1.0 / (cosd(apparent_zenith) + 0.50572((96.07995 - apparent_zenith)^(- 1.6364)))
end

relative_airmass_kastenyoung1989(s::SolarPosition) = 
    relative_airmass_kastenyoung1989(s.apparent_zenith)

"""
    relative_airmass_gueymard1993(apparent_zenith)
    relative_airmass_gueymard1993(s::SolarPosition)

Computes the relative airmass as proposed in 
[gueymard1993critical, gueymard1993development](@cite).

```math
    $VN_RELATIVE_AIRMASS = \\frac{1}{\\cos $VN_TOPOCENTRIC_APPARENT_ZENITH 
    + 0.00176759 $VN_TOPOCENTRIC_APPARENT_ZENITH((94.37515 - 
    $VN_TOPOCENTRIC_APPARENT_ZENITH)^{-1.21563})}
```
where ``$VN_TOPOCENTRIC_APPARENT_ZENITH`` is the Sun's apparent zenith, i.e., the 
complementary to 90° of the [`topocentric_apparent_elevation`](@ref).
"""
function relative_airmass_gueymard1993(apparent_zenith)
    return 1.0 / (cosd(apparent_zenith) + 0.00176759apparent_zenith*
           ((94.37515 - apparent_zenith)^(-1.21563)))
end

relative_airmass_gueymard1993(s::SolarPosition) =
    relative_airmass_gueymard1993(s.apparent_zenith)

"""
    relative_airmass_young1994(zenith)
    relative_airmass_young1994(s::SolarPosition)

Computes the relative airmass as proposed in [young1994air](@cite).

```math
    $VN_RELATIVE_AIRMASS = \\frac{1.002432\\cos^2 $VN_TOPOCENTRIC_ZENITH + 
    0.148386\\cos $VN_TOPOCENTRIC_ZENITH + 0.0096467}{\\cos^3 $VN_TOPOCENTRIC_ZENITH + 
    0.149864\\cos^2 $VN_TOPOCENTRIC_ZENITH + 0.0102963\\cos $VN_TOPOCENTRIC_ZENITH + 
    0.000303978}
```
where ``$VN_TOPOCENTRIC_ZENITH`` is the Sun's zenith, i.e., the 
complementary to 90° of the [`topocentric_elevation`](@ref).
"""
function relative_airmass_young1994(zenith)
    return (1.002432cosd(zenith)^2 + 0.148386cosd(zenith) + 0.0096467) /
           (cosd(zenith_rad)^3 + 0.149864cosd(zenith)^2 + 0.0102963cosd(zenith) + 0.000303978)
end

relative_airmass_young1994(s::SolarPosition) = 
    relative_airmass_young1994(s.zenith)

"""
    relative_airmass_pickering2002(apparent_zenith)
    relative_airmass_pickering2002(s::SolarPosition)

Computes the relative airmass as proposed in [pickering2002southern](@cite).

```math
    $VN_RELATIVE_AIRMASS = \\frac{1}{\\sin(90 - $VN_TOPOCENTRIC_APPARENT_ZENITH + 
    \\frac{244}{165 + 47(90 - $VN_TOPOCENTRIC_APPARENT_ZENITH)^{1.1}})
```
where ``$VN_TOPOCENTRIC_APPARENT_ZENITH`` is the Sun's apparent zenith, i.e., the 
complementary to 90° of the [`topocentric_apparent_elevation`](@ref).
"""
function relative_airmass_pickering2002(apparent_zenith)
    return 1.0 / sin(90 - apparent_zenith + 244.0 / (165 + 47.0(90 - apparent_zenith)^1.1))
end

relative_airmass_pickering2002(s::SolarPosition) = 
    relative_airmass_pickering2002(s.apparent_zenith)

"""
    relative_airmass_gueymard2003(apparent_zenith)
    relative_airmass_gueymard2003(s::SolarPosition)

Computes the relative airmass as proposed in [gueymard2003direct, gueymard2019clear](@cite). 

```math
    $VN_RELATIVE_AIRMASS = \\frac{1}{\\cos $VN_TOPOCENTRIC_APPARENT_ZENITH + 
    0.48353 $VN_TOPOCENTRIC_APPARENT_ZENITH^{0.095846} / (96.741 - 
    $VN_TOPOCENTRIC_APPARENT_ZENITH)^{1.754}}
```
where ``$VN_TOPOCENTRIC_APPARENT_ZENITH`` is the Sun's apparent zenith, i.e., the 
complementary to 90° of the [`topocentric_apparent_elevation`](@ref).
"""
function relative_airmass_gueymard2003(sun_apparent_elevation)
    return 1.0 / (cosd(apparent_zenith) + 0.48353apparent_zenith^0.095846 / 
           (96.741 - apparent_zenith)^1.754)
end

relative_airmass_gueymard2003(s::SolarPosition) =
    relative_airmass_gueymard2003(s.apparent_zenith)

# Vapor -----------------------------------------------------------------------------------
"""
    apparent_vapor_scale_height_gueymard94(air_temperature)

Computes the apparent water vapor scale height [km] as proposed in 
[gueymard1994analysis; Eq. 4](@cite).

```math
    $VN_APPARENT_VAPOR_SCALE_HEIGHT = 
    0.4976 + 1.5265\\frac{$VN_AIR_TEMPERATURE + $VN_ABSOLUTE_ZERO}{$VN_ABSOLUTE_ZERO} + 
    \\exp(13.6897\\frac{$VN_AIR_TEMPERATURE + $VN_ABSOLUTE_ZERO}{$VN_ABSOLUTE_ZERO} - 
    14.9188\\frac{$VN_AIR_TEMPERATURE + $VN_ABSOLUTE_ZERO}{$VN_ABSOLUTE_ZERO}^3)
```
- ``$VN_AIR_TEMPERATURE`` is the `air_temperature` [°C] 
- ``$VN_ABSOLUTE_ZERO`` is the absolute zero [°C]
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

```math
    $VN_SURFACE_VAPOR_DENSITY = \\frac{216.7 $VN_RELATIVE_HUMIDITY 
    $VN_APPARENT_VAPOR_SCALE_HEIGHT}{$VN_AIR_TEMPERATURE}
```
- ``$VN_RELATIVE_HUMIDITY`` is the `relative_humidity` [%]
- ``$VN_APPARENT_VAPOR_SCALE_HEIGHT`` is the `apparent_vapor_scale` [km], see 
    [`apparent_vapor_scale_height_gueymard94`](@ref)
- ``$VN_AIR_TEMPERATURE`` is the `air_temperature` [°C]
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

```math
    $VN_SATURATION_VAPOR_PRESSURE = 
    \\exp\\left(22.330 - 49.140 \\frac{100}{$VN_AIR_TEMPERATURE + $VN_ABSOLUTE_ZERO} 
    - 10.922 \\left(\\frac{100}{$VN_AIR_TEMPERATURE + $VN_ABSOLUTE_ZERO}\\right)^2 - 
    0.39015 \\frac{$VN_AIR_TEMPERATURE + $VN_ABSOLUTE_ZERO}{100} \\right)
```
- ``$VN_AIR_TEMPERATURE`` is the `air_temperature` [°C].
- ``$VN_ABSOLUTE_ZERO`` is the absolute zero [°C].
"""
function saturation_vapor_pressure_gueymard93(air_temperature)
    T_kelvin = air_temperature + ABSOLUTE_ZERO
    return exp(22.330 - 49.140(100/T_kelvin) - 10.922(100/T_kelvin)^2 - 
               0.39015T_kelvin/100)
end

"""
    saturation_vapor_pressure_oke2018(dew_temperature)

Computes the saturation water vapor pressure [mbar] as proposed in [oke2018guide](@cite).

```math
    $VN_SATURATION_VAPOR_PRESSURE = 
    6.112\\exp(17.62$VN_DEW_TEMPERATURE / (243.12 + $VN_DEW_TEMPERATURE))
```
where ``$VN_DEW_TEMPERATURE`` is the `dew_temperature` [°C].
"""
function saturation_vapor_pressure_oke2018(dew_temperature)
    return 6.112exp(17.62dew_temperature / (243.12 + dew_temperature)) 
end

"""
    vapor_pressure_oke2018(air_temperature)

Computes the water vapor pressure [mbar] as proposed in [oke2018guide](@cite). 

```math
    $VN_VAPOR_PRESSURE = 6.112\\exp\\left(
    \\frac{17.62$VN_AIR_TEMPERATURE}{243.12 + $VN_AIR_TEMPERATURE})\\right)
```
where ``$VN_AIR_TEMPERATURE`` is the `air_temperature` [°C].
"""
function vapor_pressure_oke2018(air_temperature)
    return 6.112exp(17.62air_temperature / (243.12 + air_temperature))
end

# Precipitable water and humidity ---------------------------------------------------------
"""
    precipitable_water_gueymard94(
        saturation_vapor_pressure,
        surface_vapor_density
    )

Computes precipitable water [cm] as proposed in [gueymard1994analysis; Eq. 1](@cite).

The accuracy of this method is approximately 20% for moderate precipitable water (1-3 cm) 
and less accurate otherwise. 
```math
        $VN_PRECIPITABLE_WATER = 0.1 $VN_SURFACE_VAPOR_DENSITY 
        $VN_SATURATION_VAPOR_PRESSURE
```
- ``$VN_SURFACE_VAPOR_DENSITY`` is the `surface_vapor_density` [km], 
    see [`surface_vapor_density_gueymard94`](@ref)
- ``$VN_SATURATION_VAPOR_PRESSURE`` is the saturation water vapor pressure [mbar], 
    see [`saturation_vapor_pressure_oke2018`](@ref), 
    [`saturation_vapor_pressure_gueymard93`](@ref)
"""
function precipitable_water_gueymard94(
    saturation_vapor_pressure,
    surface_vapor_density
)
    return max(0.1saturation_vapor_pressure * surface_vapor_density, 0.1)
end

"""
    relative_humidity_oke2018(
        vapor_pressure,
        saturation_vapor_pressure
    )

Computes the relative humidity [%] as proposed in [oke2018guide](@cite). 

```math
    $VN_RELATIVE_HUMIDITY = 100 \\frac{$VN_SATURATION_VAPOR_PRESSURE}{$VN_VAPOR_PRESSURE}
```
- ``$VN_SATURATION_VAPOR_PRESSURE`` is the saturation water vapor pressure [mbar], 
    see [`saturation_vapor_pressure_oke2018`](@ref), 
    [`saturation_vapor_pressure_gueymard93`](@ref) 
- ``$VN_VAPOR_PRESSURE`` is the water vapor pressure [mbar], see 
    [`vapor_pressure_oke2018`](@ref)
"""
function relative_humidity_oke2018(vapor_pressure, saturation_vapor_pressure)
    return 100saturation_vapor_pressure / vapor_pressure
end
