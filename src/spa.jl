using Pkg.Artifacts, Serialization, Dates

const EARTH_PERIODIC_TERMS = open(artifact"data/earth-periodic-terms.jld", "r") do io
    deserialize(io)
end

# Limits an angle (in degrees) between 0° and 360°
mod360(angle) = mod(angle, 360)

"""
    delta_T(datetime::DateTime)

Computes the difference between Terrestrial Dynamical Time (TD) and Universal Time (UT).

It is evaluated as a piecewise polynomial, as per 
[NASA](https://eclipse.gsfc.nasa.gov/SEcat5/deltatpoly.html), which in turn refer to 
[morrison2004historical, stephenson1986atlas](@cite). Note that ``\\Delta T`` is unknown 
for years before `-1999` and after `3000`. Values could be calculated for such intervals, 
although they are not intended to be used for these years.
"""
function delta_T(datetime::DateTime)
    year, month = yearmonth(datetime) 
    y = year + (month - 0.5) / 12

    if year < -500
        return -20 + 32((y - 1820) / 100)^2
    elseif -500 <= year < 500
        t = y / 100
        return 10583.6 - 1014.41t + 33.78311t^2 - 5.952053t^3 -
               0.1798452t^4 + 0.022174192t^5 + 0.0090316521t^6
    elseif 500 <= year < 1600
        t = (y - 1000) / 100
        return 1574.2 - 556.01t + 71.23472t^2 + 0.319781t^3 -
               0.8503463t^4 - 0.005050998t^5 + 0.0083572073t^6
    elseif 1600 <= year < 1700
        t = y - 1600
        return 120 - 0.9808t - 0.01532t^2 + t^3 / 7129
    elseif 1700 <= year < 1800
        t = y - 1700
        return 8.83 + 0.1603t - 0.0059285t^2 + 0.00013336t^3
               t^4 / 1174000
    elseif 1800 <= year < 1860
        t = y - 1800
        return 13.72 - 0.332447t + 0.0068612t^2 + 0.0041116t^3 -
               0.00037436t^4 + 0.0000121272t^5 - 0.0000001699t^6 + 
               0.000000000875t^7
    elseif 1860 <= year < 1900
        t = y - 1860
        return 7.62 + 0.5737t - 0.251754t^2 + 0.01680668t^3 - 
               0.0004473624t^4 + t^5 / 233174
    elseif 1900 <= year < 1920
        t = y - 1900
        return -2.79 + 1.494119t - 0.0598939*t^2 + 0.0061966*t^3 - 
                0.000197*t^4
    elseif 1920 <= year < 1941
        t = y - 1920
        return 21.20 + 0.84493t - 0.076100t^2 + 0.0020936t^3
    elseif 1941 <= year < 1961
        t = y - 1950
        return 29.07 + 0.407t - t^2 / 233 + t^3 / 2547
    elseif 1961 <= year < 1986
        t = y - 1975
        return 45.45 + 1.067t - t^2 / 260 - t^3 / 718
    elseif 1986 <= year < 2005
        t = y - 2000
        return 63.86 + 0.3345t - 0.060374t^2 + 0.0017275t^3 + 
               0.000651814t^4 + 0.00002373599t^5
    elseif 2005 <= year < 2050
        t = y - 2000
        return 62.92 + 0.32217t + 0.005589t^2
    elseif 2050 <= year < 2150
        return -20 + 32((y - 1820)/100)^2 - 0.5628(2150 - y)
    elseif year >= 2150
        return -20 + 32((y-1820)/100)^2
    end
end

"""
    julian_ephemeris_day(julian_day, ΔT)

A continuous count of days measured in uniform Ephemeris Time (or its successors).

```math
$VN_JULIAN_EPHEMERIS_DAY = $VN_JULIAN_DAY + \\frac{\\Delta T}{86400}
```
where:
- ``$VN_JULIAN_DAY`` is the [Julian day](https://en.wikipedia.org/wiki/Julian_day), which
    can be computed using 
    [`Dates.datetime2julian`](https://docs.julialang.org/en/v1/stdlib/Dates/#Dates.datetime2julian)
- ``$VN_DELTA_T`` is the difference between Terrestrial Dynamical Time and Universal Time,
    which can be computed using [`delta_T`](@ref).  
""" 
function julian_ephemeris_day(julian_day, ΔT)
    return julian_day + ΔT / 86400.0
end

"""
    julian_century(julian_day)

A time interval of exactly 36,525 days (365.25 days × 100) used as a standard unit in 
astronomy.
 
```math
$VN_JULIAN_CENTURY = \\frac{$VN_JULIAN_DAY - 2451545}{36525}
```
where ``$VN_JULIAN_DAY`` is the [Julian day](https://en.wikipedia.org/wiki/Julian_day), 
which can be computed using 
[`Dates.datetime2julian`](https://docs.julialang.org/en/v1/stdlib/Dates/#Dates.datetime2julian).
"""
function julian_century(julian_day)
    return (julian_day - 2451545.0) / 36525.0
end

"""
    julian_ephemeris_century(julian_ephemeris_day)

A 36,525-day interval measured specifically in Ephemeris Time (or its modern dynamical time
scales) for high-precision astronomical modeling. 
 
```math
$VN_JULIAN_EPHEMERIS_CENTURY = \\frac{$VN_JULIAN_EPHEMERIS_DAY - 2451545}{36525}
``` 
where ``$VN_JULIAN_EPHEMERIS_DAY`` is the [`julian_ephemeris_day`](@ref).
"""
function julian_ephemeris_century(julian_ephemeris_day)
    return (julian_ephemeris_day − 2451545.0) / 36525.0
end

"""
    julian_ephemeris_millenium(julian_ephemeris_century)

A 1,000-year interval equal to 365,250 ephemeris days, defined within Ephemeris Time for 
long-term astronomical calculations. 

```math
$VN_JULIAN_EPHEMERIS_MILLENIUM = \\frac{$VN_JULIAN_EPHEMERIS_CENTURY}{10}
``` 
where ``$VN_JULIAN_EPHEMERIS_CENTURY`` is the [`julian_ephemeris_century`](@ref). 
"""
function julian_ephemeris_millenium(julian_ephemeris_century)
    return julian_ephemeris_century / 10.0
end

"""
    heliocentric_polynomial(julian_ephemeris_millenium, coefficients::AbstractVector)

Polynomial approximation for the heliocentric longitude, latitude, radius.
"""
function heliocentric_polynomial(julian_ephemeris_millenium, coefficients::AbstractVector)
    polynomial = 0
    jem = julian_ephemeris_millenium
    for (i, C) in enumerate(coefficients)
        polynomial += sum(C[:, 1] .* cos.(C[:, 2] .+ C[:, 3] .* jem)) * jem^(i - 1)
    end
    return polynomial / 1e8
end

"""
    heliocentric_longitude(julian_ephemeris_millenium)

A celestial object's angular distance north or south of the ecliptic plane as measured from 
the Sun. 

It is calculated as a polynomial expression, where the coefficients are collected in 
``\\mathcal{M} = \\{\\mathbf{M}_i\\}_{i=1}^6``, with 
``\\mathbf{M}_i \\in \\mathbb{R}^{n_i \\times 3}`` being the coeffient matrix for the 
``i``-th polynomal coefficient ``c_{\\mu, i}``. 

```math
    \\begin{align*}
        c_{\\mu, i} &= \\sum_{j = 1}^{n_i} M^{(j,1)}_i \\cos(M^{(j,2)}_i + M^{(j,3)}_i 
        $VN_JULIAN_EPHEMERIS_MILLENIUM) \\\\
        $VN_HELIOCENTRIC_LONGITUDE &= \\frac{180}{\\pi} \\frac{1}{1e8} \\sum_{i = 1}^6 
        c_{\\ell, i} $VN_JULIAN_EPHEMERIS_MILLENIUM^{i-1} \\mod 360
    \\end{align*}
```
where ``$VN_JULIAN_EPHEMERIS_MILLENIUM`` is the [`julian_ephemeris_millenium`](@ref). 
The ``180/\\pi`` term coverts ``$VN_HELIOCENTRIC_LONGITUDE`` from radians to degrees. 

The ``\\mod(\\cdot)`` operator ensures that the results remains in the `[0, 360]`
range. The coefficient matrices in ``\\mathcal{M}`` are collected in 
`SolarFunctions.EARTH_PERIODIC_TERMS.M`. 
"""
function heliocentric_longitude(julian_ephemeris_millenium)
    μ = heliocentric_polynomial(
            julian_ephemeris_millenium,
            EARTH_PERIODIC_TERMS.M
        )
    return μ |> rad2deg |> mod360
end

"""
    heliocentric_latitude(julian_ephemeris_millenium)

A celestial object's angular position around the Sun measured in the ecliptic plane from a 
defined reference direction.

It is calculated as a polynomial expression, where the coefficients are collected in 
``\\mathcal{L} = \\{\\mathbf{L}_i\\}_{i=1}^2``, with ``\\mathbf{L}_i \\in \\mathbb{R}^{n_i 
\\times 3}`` being the coeffient matrix for the ``i``-th polynomal coefficient 
``c_{\\lambda, i}``. 

```math
    \\begin{align*}
        c_{\\lambda, i} &= \\sum_{j = 1}^{n_i} L^{(j,1)}_i \\cos(L^{(j,2)}_i + L^{(j,3)}_i 
        $VN_JULIAN_EPHEMERIS_MILLENIUM) \\\\
        $VN_HELIOCENTRIC_LATITUDE &= \\frac{180}{\\pi} \\frac{1}{1e8} \\sum_{i = 1}^2 
        c_{b, i} $VN_JULIAN_EPHEMERIS_MILLENIUM^{i-1}
    \\end{align*}
```
where ``$VN_JULIAN_EPHEMERIS_MILLENIUM`` is the [`julian_ephemeris_millenium`](@ref). The 
``180/\\pi`` term coverts ``$VN_HELIOCENTRIC_LATITUDE`` from radians to degrees. 

The coefficient matrices in ``\\mathcal{L}`` are collected in 
`SolarFunctions.EARTH_PERIODIC_TERMS.L`. 
"""
function heliocentric_latitude(julian_ephemeris_millenium)
    λ = heliocentric_polynomial(
            julian_ephemeris_millenium,
            EARTH_PERIODIC_TERMS.L
        )
    return λ |> rad2deg
end

"""
    heliocentric_radius(julian_ephemeris_millenium)

The distance from the Sun to the celestial object in a chosen heliocentric coordinate 
system. In this case, it corresponds to the Earth-Sun distance.

It is calculated as a polynomial expression, where the coefficients are collected in 
``\\mathcal{R} = \\{\\mathbf{R}_i\\}_{i=1}^5``, with ``\\mathbf{R}_i \\in 
\\mathbb{R}^{n_i \\times 3}`` being the coeffient matrix for the ``i``-th polynomal 
coefficient ``c_{R, i}``. 

```math
    \\begin{align*}
        c_{R, i} &= \\sum_{j = 1}^{n_i} R^{(j,1)}_i \\cos(R^{(j,2)}_i + R^{(j,3)}_i 
        $VN_JULIAN_EPHEMERIS_MILLENIUM) \\\\
        $VN_HELIOCENTRIC_RADIUS &= \\frac{180}{\\pi} \\frac{1}{1e8} \\sum_{i = 1}^5 
        c_{r, i} $VN_JULIAN_EPHEMERIS_MILLENIUM^{i-1}
    \\end{align*}
```
where ``$VN_JULIAN_EPHEMERIS_MILLENIUM`` is the [`julian_ephemeris_millenium`](@ref). 

The coefficient matrices in ``\\mathcal{R}`` are collected in 
`SolarFunctions.EARTH_PERIODIC_TERMS.R`. 
"""
function heliocentric_radius(julian_ephemeris_millenium)
    return heliocentric_polynomial(
               julian_ephemeris_millenium, 
               EARTH_PERIODIC_TERMS.R
           )
end

"""
    geocentric_longitude(heliocentric_longitude)

A celestial object's ecliptic longitude as measured from Earth’s center, referenced to the 
ecliptic plane and the vernal equinox. 
 
```math
$VN_GEOCENTRIC_LONGITUDE = $VN_HELIOCENTRIC_LONGITUDE + 180 \\mod 360
```
where ``$VN_HELIOCENTRIC_LONGITUDE`` is the [`heliocentric_longitude`](@ref), and the 
``\\mod(\\cdot)`` operator ensures that the results remains in the `[0, 360]` range.
"""
function geocentric_longitude(heliocentric_longitude)
    return mod360(heliocentric_longitude + 180)
end

"""
    geocentric_latitude(heliocentric_latitude)

A celestial object's angular distance north or south of the ecliptic plane as measured from 
Earth’s center. 

```math
    $VN_GEOCENTRIC_LATITUDE = -$VN_HELIOCENTRIC_LATITUDE
```
where ``$VN_HELIOCENTRIC_LATITUDE`` is the [`heliocentric_latitude`](@ref). 
"""
function geocentric_latitude(heliocentric_latitude)
    return -heliocentric_latitude
end

"""
    geocentric_sun_ascension(
        apparent_sun_longitude,
        elliptic_obliquity,
        geocentric_latitude
    )

The angle measured eastward from the geocentric meridian, that locates a direction or point 
relative to Earth’s ellipsoid rather than its rotational axis. 

```math
    \\sigma_{\\text{geo}} = \\text{atan2}(\\sin$VN_APPARENT_SUN_LONGITUDE 
    \\cos$VN_ELLIPTIC_OBLIQUITY - \\tan$VN_GEOCENTRIC_LATITUDE\\sin$VN_ELLIPTIC_OBLIQUITY, 
    \\cos$VN_APPARENT_SUN_LONGITUDE)
```
where
- ``$VN_APPARENT_SUN_LONGITUDE`` corresponds to [`apparent_sun_longitude`](@ref) 
- ``$VN_ELLIPTIC_OBLIQUITY`` corresponds to [`elliptic_obliquity`](@ref)
- ``$VN_GEOCENTRIC_LATITUDE`` corresponds to [`geocentric_latitude`](@ref)
"""
function geocentric_sun_ascension(
    apparent_sun_longitude,
    elliptic_obliquity,
    geocentric_latitude
)
    s = sind(apparent_sun_longitude) * cosd(elliptic_obliquity) -
        tand(geocentric_latitude) * sind(elliptic_obliquity)
    c = cosd(apparent_sun_longitude)

    return atand(s, c) |> mod360
end

"""
    geocentric_sun_declination(
        apparent_sun_longitude,
        elliptic_obliquity,
        geocentric_latitude
    )

The Sun’s angular distance, north or south of Earth’s equatorial plane as measured from 
Earth’s center. 

```math
    $VN_GEOCENTRIC_SUN_DECLINATION = \\text{asin}(\\sin $VN_GEOCENTRIC_LATITUDE 
    \\cos $VN_ELLIPTIC_OBLIQUITY + \\cos $VN_GEOCENTRIC_LATITUDE \\sin $VN_ELLIPTIC_OBLIQUITY 
    \\sin $VN_APPARENT_SUN_LONGITUDE) \\mod 360
``` 
where:
- ``$VN_APPARENT_SUN_LONGITUDE`` corresponds to [`apparent_sun_longitude`](@ref)
- ``$VN_ELLIPTIC_OBLIQUITY`` corresponds to [`elliptic_obliquity`](@ref)
- ``$VN_GEOCENTRIC_LATITUDE`` corresponds to [`geocentric_latitude`](@ref)

The ``\\mod(\\cdot)`` operator ensures that the results remains in the `[0, 360]` range.
"""
function geocentric_sun_declination(
    apparent_sun_longitude,
    elliptic_obliquity,
    geocentric_latitude,
)
    return asind(
        sind(geocentric_latitude) * cosd(elliptic_obliquity) +
        cosd(geocentric_latitude) * sind(elliptic_obliquity) * sind(apparent_sun_longitude)
    )
end

"""
    mean_moon_elongation_from_sun(julian_ephemeris_century)

The average angular separation between the Moon and the Sun as measured along the ecliptic.

```math
    $VN_MEAN_MOON_ELONGATION_FROM_SUN = 297.85036 + 445267.111480$VN_JULIAN_EPHEMERIS_CENTURY 
    − 0.0019142$VN_JULIAN_EPHEMERIS_CENTURY^2 + \\frac{$VN_JULIAN_EPHEMERIS_CENTURY^3}{189474}
```
where ``$VN_JULIAN_EPHEMERIS_CENTURY`` is the [`julian_ephemeris_century`](@ref). 
"""
function mean_moon_elongation_from_sun(julian_ephemeris_century)
    jec = julian_ephemeris_century
    return 297.85036 + 445267.111480jec − 0.0019142jec ^ 2 + jec^3 / 189474.0
end

"""
    mean_sun_anomaly(julian_ephemeris_century)

The angular position of the Sun in its elliptical orbit, measured from perihelion and 
increasing uniformly in time. 

```math
    $VN_MEAN_SUN_ANOMALY = 357.52772 + 35999.050340$VN_JULIAN_EPHEMERIS_CENTURY 
    − 0.0001603$VN_JULIAN_EPHEMERIS_CENTURY^2 − \\frac{$VN_JULIAN_EPHEMERIS_CENTURY^3}
    {300000}
```
where ``$VN_JULIAN_EPHEMERIS_CENTURY`` is the [`julian_ephemeris_century`](@ref).
"""
function mean_sun_anomaly(julian_ephemeris_century)
    jec = julian_ephemeris_century
    return 357.52772 + 35999.050340jec − 0.0001603jec^2 − jec^3 / 300000.0
end

"""
    mean_moon_anomaly(julian_ephemeris_century)

The uniformly increasing angular position of the Moon in its elliptical orbit measured from 
perigee.

```math
    $VN_MEAN_MOON_ANOMALY = 134.96298 + 477198.867398$VN_JULIAN_EPHEMERIS_CENTURY +
    0.0086972$VN_JULIAN_EPHEMERIS_CENTURY^2 + \\frac{$VN_JULIAN_EPHEMERIS_CENTURY^3}{56250}
```
where ``$VN_JULIAN_EPHEMERIS_CENTURY`` is the [`julian_ephemeris_century`](@ref).
"""
function mean_moon_anomaly(julian_ephemeris_century)
    jec = julian_ephemeris_century
    return 134.96298 + 477198.867398jec + 0.0086972jec^2 + jec^3 / 56250.0
end

"""
    moon_latitude_argument(julian_ephemeris_century)

The angle from the Moon’s ascending node to its position measured along its orbit, using 
mean (unperturbed) orbital elements. 

```math
    $VN_MOON_LATITUDE_ARGUMENT = 93.27191 + 483202.017538$VN_JULIAN_EPHEMERIS_CENTURY - 
    0.0036825$VN_JULIAN_EPHEMERIS_CENTURY^2 + \\frac{$VN_JULIAN_EPHEMERIS_CENTURY^3}{327270}
```
where ``$VN_JULIAN_EPHEMERIS_CENTURY`` is the [`julian_ephemeris_century`](@ref).
"""
function moon_latitude_argument(julian_ephemeris_century)
    jec = julian_ephemeris_century
    return 93.27191 + 483202.017538jec - 0.0036825jec^2 + jec^3 / 327270.0
end

"""
    ascending_moon_longitude(julian_ephemeris_century)

The ecliptic longitude referenced to the mean equinox of the date, of the point where the 
Moon’s mean orbit crosses northward through the ecliptic. 

```math
    $VN_ASCENDING_MOON_LONGITUDE = 125.04452 - 1934.136261$VN_JULIAN_EPHEMERIS_CENTURY + 
    0.0020708$VN_JULIAN_EPHEMERIS_CENTURY^2 + \\frac{$VN_JULIAN_EPHEMERIS_CENTURY^3}{450000}
```
where ``$VN_JULIAN_EPHEMERIS_CENTURY`` is the [`julian_ephemeris_century`](@ref).
"""
function ascending_moon_longitude(julian_ephemeris_century)
    jec = julian_ephemeris_century
    return 125.04452 - 1934.136261jec + 0.0020708jec^2 + jec^3 / 450000.0
end

"""
    nutation_coefficients(julian_ephemeris_century)

Constructs the vector of weights for evaluating the [`nutation_longitude`](@ref) and 
[`nutation_obliquity`](@ref). 

```math
    \\begin{align*}
        w_i &= \\sum_{j=1}^5 X_j Y_{ij} \\\\
        \\mathbf{X} &= [$VN_MEAN_MOON_ELONGATION_FROM_SUN, $VN_MEAN_SUN_ANOMALY, 
        $VN_MEAN_MOON_ANOMALY, $VN_MOON_LATITUDE_ARGUMENT, $VN_ASCENDING_MOON_LONGITUDE]^\\top 
    \\end{align*}
```
The coefficients in matrix ``\\mathbf{Y} \\in \\mathbb{R}^{m \\times 5}`` are collected in 
`SolarFunctions.EARTH_PERIODIC_TERMS.Y`. 
"""
function nutation_coefficients(julian_ephemeris_century)
    X = [
        mean_moon_elongation_from_sun(julian_ephemeris_century),
        mean_sun_anomaly(julian_ephemeris_century),
        mean_moon_anomaly(julian_ephemeris_century),
        moon_latitude_argument(julian_ephemeris_century),
        ascending_moon_longitude(julian_ephemeris_century)
    ]
    return EARTH_PERIODIC_TERMS.Y * X
end

"""
    nutation_longitude(julian_ephemeris_century, nutation_coefficients::AbstractVector)

The small periodic variation in Earth’s ecliptic longitude caused by gravitational torques 
from the Moon and Sun. 

```math
    $VN_NUTATION_LONGITUDE = \\frac{1}{3.6e7} \\sum_{i=1}^n (\\Psi_{i,1} + \\Psi_{i,2} 
    $VN_JULIAN_EPHEMERIS_CENTURY) \\sin w_i
```
where:
- ``\\boldsymbol{\\Psi} \\in \\mathbb{R}^{m \\times 2}`` is the coefficient matrix, whose 
    terms are collected in `SolarFunctions.EARTH_PERIODIC_TERMS.ψ`
- ``$VN_JULIAN_EPHEMERIS_CENTURY`` is the [`julian_ephemeris_century`](@ref).
- ``w_i`` is the ``i``-th term of the [`nutation_coefficients`](@ref) vector
"""
function nutation_longitude(julian_ephemeris_century, nutation_coefficients::AbstractVector)
    jec = julian_ephemeris_century
    Ψ = EARTH_PERIODIC_TERMS.Ψ
    W = nutation_coefficients
    return sum((Ψ[:, 1] .+ Ψ[:, 2] * jec) .* sind.(W)) / 3.6e7
end

"""
    nutation_obliquity(julian_ephemeris_century, nutation_coefficients::AbstractVector)

The small periodic variation in Earth’s axial tilt (obliquity), ``$VN_NUTATION_OBLIQUITY``, 
resulting from gravitational perturbations. 

```math
    $VN_NUTATION_OBLIQUITY = \\frac{1}{3.6e7} \\sum_{i=1}^n (E_{i,1} + E_{i,2} 
    $VN_JULIAN_EPHEMERIS_CENTURY) \\cos w_i
```
where: 
- ``\\boldsymbol{E} \\in \\mathbb{R}^{m \\times 2}`` is the coefficient matrix, whose terms 
    are collected in `SolarFunctions.EARTH_PERIODIC_TERMS.ϵ`
- ``$VN_JULIAN_EPHEMERIS_CENTURY`` is the [`julian_ephemeris_century`](@ref).
- ``w_i`` is the ``i``-th term of the [`nutation_coefficients`](@ref) vector
"""
function nutation_obliquity(julian_ephemeris_century, nutation_coefficients::AbstractVector)
    jec = julian_ephemeris_century
    E = EARTH_PERIODIC_TERMS.E
    W = nutation_coefficients
    return sum((E[:, 1] .+ E[:, 2] * jec) .* cosd.(W)) / 3.6e7
end

"""
    mean_elliptic_obliquity(julian_ephemeris_millenium)

The angle between Earth’s equator and the mean (long-term averaged) ecliptic defined by an 
unperturbed elliptical Earth orbit. 

```math
    \\begin{align*}
        U &= \\frac{$VN_JULIAN_EPHEMERIS_MILLENIUM}{10} \\\\
        $VN_MEAN_ELLIPTIC_OBLIQUITY &= 84381.448 − 4680.93U − 155.0U^2 + 1999.25U^3 − 
        51.38U^4 − \\\\
        & \\qquad 249.67U^5 − 39.05U^6 + 7.12U^7 + 27.87U^8 + 5.79U^9 + 2.45U^{10}
    \\end{align*}
```
where ``$VN_JULIAN_EPHEMERIS_MILLENIUM`` is the [`julian_ephemeris_millenium`](@ref). 
"""
function mean_elliptic_obliquity(julian_ephemeris_millenium)
    U = julian_ephemeris_millenium / 10.0
    return 84381.448 − 4680.93U − 155.0U^2 + 1999.25U^3 − 51.38U^4 − 
           249.67U^5 − 39.05U^6 + 7.12U^7 + 27.87U^8 + 5.79U^9 + 2.45U^10
end

"""
    elliptic_obliquity(mean_elliptic_obliquity, nutation_obliquity)

The angle between Earth’s equator and the mean ecliptic defined as an ideal, unperturbed 
ellipse. 

```math
$VN_ELLIPTIC_OBLIQUITY = $VN_MEAN_ELLIPTIC_OBLIQUITY / 3600 + $VN_NUTATION_OBLIQUITY
```
where:
- ``$VN_MEAN_ELLIPTIC_OBLIQUITY`` corresponds to [`mean_elliptic_obliquity`](@ref)
- ``$VN_NUTATION_OBLIQUITY`` corresponds to [`nutation_obliquity`](@ref)
"""
function elliptic_obliquity(mean_elliptic_obliquity, nutation_obliquity)
    return mean_elliptic_obliquity / 3600.0 + nutation_obliquity
end

"""
    aberration_correction(heliocentric_radius)

An adjustment applied to celestial coordinates to account for the apparent displacement 
caused by Earth’s motion through space. 

```math
$VN_ABERRATION_CORRECTION = -\\frac{20.4898}{3600 $VN_HELIOCENTRIC_RADIUS}
```
where ``$VN_HELIOCENTRIC_RADIUS`` corresponds to the [`heliocentric_radius`](@ref) [AU].
"""
function aberration_correction(heliocentric_radius)
    return -20.4898 / (3600.0heliocentric_radius)
end

"""
    apparent_sun_longitude(geocentric_longitude, nutation_longitude, aberration_correction)

The Sun’s ecliptic longitude after applying corrections for nutation and aberration.

```math
$VN_APPARENT_SUN_LONGITUDE = $VN_GEOCENTRIC_LONGITUDE + $VN_NUTATION_LONGITUDE + 
$VN_ABERRATION_CORRECTION
```
where:
- ``$VN_GEOCENTRIC_LONGITUDE`` corresponds to [`geocentric_longitude`](@ref)
- ``$VN_NUTATION_LONGITUDE`` corresponds to [`nutation_longitude`](@ref)
- ``$VN_ABERRATION_CORRECTION`` corresponds to [`aberration_correction`](@ref)
"""
function apparent_sun_longitude(
    geocentric_longitude,
    nutation_longitude,
    aberration_correction
)
    return geocentric_longitude + nutation_longitude + aberration_correction
end

"""
    mean_sun_longitude(julian_ephemeris_millenium)

The Sun’s ecliptic longitude calculated from a uniformly moving fictitious Sun on the 
ecliptic. 

```math
    $VN_MEAN_SUN_LONGITUDE = 280.4664567 + 360007.6982779$VN_JULIAN_EPHEMERIS_MILLENIUM + 
    0.03032028$VN_JULIAN_EPHEMERIS_MILLENIUM^2 + \\frac{$VN_JULIAN_EPHEMERIS_MILLENIUM^3}{49931} 
    - \\frac{$VN_JULIAN_EPHEMERIS_MILLENIUM^4}{15300} 
    - \\frac{$VN_JULIAN_EPHEMERIS_MILLENIUM^5}{2000000}
```
where ``$VN_JULIAN_EPHEMERIS_MILLENIUM`` corresponds to the 
[`julian_ephemeris_millenium`](@ref).
"""
function mean_sun_longitude(julian_ephemeris_millenium)
    jem = julian_ephemeris_millenium
    return 280.4664567 + 360007.6982779jem + 
           0.03032028jem^2 + jem^3 / 49931 -
           jem^4 / 15300 - jem^5 / 2000000
end

"""
    mean_sidereal_greenwich_time(julian_day, julian_century)

The hour angle of the mean vernal equinox at the Greenwich meridian reflecting Earth’s 
rotation relative to the fixed stars without nutation effects. 

```math
    $VN_MEAN_SIDEREAL_GREENWICH_TIME = 280.46061837 + 360.98564736629($VN_JULIAN_DAY 
    − 2451545.0) + 0.000387933$VN_JULIAN_CENTURY^2 - \\frac{$VN_JULIAN_CENTURY^3}{38710000} 
    \\mod 360
```
where ``$VN_JULIAN_DAY`` is the [`Julian Day`](https://en.wikipedia.org/wiki/Julian_day), 
which can be computed using 
[`Dates.datetime2julian`](https://docs.julialang.org/en/v1/stdlib/Dates/#Dates.datetime2julian). 
The ``\\mod(\\cdot)`` operator ensures that the results remains in the `[0, 360]` 
range.
"""
function mean_sidereal_greenwich_time(julian_day, julian_century) 
    return mod360(
        280.46061837 + 360.98564736629 * (julian_day − 2451545.0) + 
        0.000387933julian_century^2 - julian_century^3 / 38710000.0
    )
end

"""
    apparent_sidereal_greenwich_time(
        mean_sidereal_greenwich_time,
        nutation_longitude, 
        elliptic_obliquity
    )

Greenwich sidereal time corrected for the effects of nutation giving the true rotation 
angle relative to the apparent equinox. 

```math
$VN_APPARENT_SIDEREAL_GREENWICH_TIME = $VN_MEAN_SIDEREAL_GREENWICH_TIME + 
$VN_NUTATION_LONGITUDE \\cos $VN_ELLIPTIC_OBLIQUITY
```
where:
- ``$VN_MEAN_SIDEREAL_GREENWICH_TIME`` corresponds to [`mean_sidereal_greenwich_time`](@ref) 
- ``$VN_NUTATION_LONGITUDE ``corresponds to [`nutation_longitude`](@ref) 
- ``$VN_ELLIPTIC_OBLIQUITY ``corresponds to [`elliptic_obliquity`](@ref) 
"""
function apparent_sidereal_greenwich_time(
    mean_sidereal_greenwich_time,
    nutation_longitude, 
    elliptic_obliquity
)
    return mean_sidereal_greenwich_time + nutation_longitude * cosd(elliptic_obliquity)
end

"""
    observer_local_hour(
        observer_longitude, 
        apparent_sidereal_greenwich_time, 
        geocentric_sun_ascension
    )

The angle between the observer’s local meridian and a celestial object measured westward on 
the celestial sphere. 

```math
$VN_OBSERVER_LOCAL_HOUR = $VN_APPARENT_SIDEREAL_GREENWICH_TIME + $VN_OBSERVER_LONGITUDE - 
$VN_GEOCENTRIC_SUN_ASCENSION \\mod 360
```
where:
- ``$VN_OBSERVER_LONGITUDE`` corresponds to the `observer_longitude`
- ``$VN_APPARENT_SIDEREAL_GREENWICH_TIME`` corresponds to 
    [`apparent_sidereal_greenwich_time`](@ref)
- ``$VN_GEOCENTRIC_SUN_ASCENSION`` corresponds to [`geocentric_sun_ascension`](@ref)
The ``\\mod(\\cdot)`` operator ensures that the results remains in the 
`[0, 360]` range.
"""
function observer_local_hour(
    observer_longitude, 
    apparent_sidereal_greenwich_time, 
    geocentric_sun_ascension
)
    return mod360(
        apparent_sidereal_greenwich_time + 
        observer_longitude - 
        geocentric_sun_ascension
    )
end

"""
    reduced_observer_latitude(observer_latitude)

The angle whose tangent equals the tangent of the geodetic latitude scaled by Earth’s 
polar-to-equatorial radius ratio, representing the point’s projection onto the surface of 
the reference ellipsoid. 

```math
$VN_REDUCED_OBSERVER_LATITUDE = \\text{atan}(0.99664719\\tan($VN_OBSERVER_LATITUDE))
```
where ``$VN_OBSERVER_LATITUDE`` corresponds to the `observer_latitude`.
"""
function reduced_observer_latitude(observer_latitude)
    return atand(0.99664719 * tand(observer_latitude))
end

"""
    radial_distance_equatorial_plane(
        observer_latitude,
        observer_altitude,
        reduced_observer_latitude
    )

The perpendicular distance of a point from Earth’s equatorial plane. 

```math
$VN_RADIAL_DISTANCE_FROM_EQUATORIAL_PLANE = \\cos $VN_REDUCED_OBSERVER_LATITUDE + 
$VN_OBSERVER_ALTITUDE \\frac{\\cos $VN_OBSERVER_LATITUDE}{6378140}
```
where:
- ``$VN_OBSERVER_LATITUDE`` corresponds to `observer_latitude`
- ``$VN_OBSERVER_ALTITUDE`` corresponds to `observer_altitude`
- ``$VN_REDUCED_OBSERVER_LATITUDE`` corresponds to [`reduced_observer_latitude`](@ref)
"""
function radial_distance_equatorial_plane(
    observer_latitude,
    observer_altitude,
    reduced_observer_latitude
)
    return cosd(reduced_observer_latitude) + 
           observer_altitude * cosd(observer_latitude) / 6378140.0
end

"""
    radial_distance_rotational_axis(
        observer_latitude, 
        observer_altitude,
        reduced_observer_latitude
    )

The perpendicular distance of a point from Earth’s spin axis. 

```math
$VN_RADIAL_DISTANCE_FROM_ROTATIONAL_AXIS = 0.99664719\\sin $VN_REDUCED_OBSERVER_LATITUDE + 
$VN_OBSERVER_ALTITUDE \\frac{\\sin $VN_OBSERVER_LATITUDE}{6378140}
``` 
where:
- ``$VN_OBSERVER_LATITUDE`` corresponds to `observer_latitude`
- ``$VN_OBSERVER_ALTITUDE`` corresponds to `observer_altitude`
- ``$VN_REDUCED_OBSERVER_LATITUDE`` corresponds to [`reduced_observer_latitude`](@ref)
"""
function radial_distance_rotational_axis(
    observer_latitude, 
    observer_altitude,
    reduced_observer_latitude
)
    return 0.99664719sind(reduced_observer_latitude) + 
           observer_altitude * sind(observer_latitude) / 6378140.0
end

"""
    sun_equatorial_horizontal_parallax(heliocentric_radius)

The angle between the Sun’s direction as seen from Earth’s center and from a point on the 
equator at sea level when the Sun is on the horizon. 
```math
$VN_SUN_EQUATORIAL_HORIZONTAL_PARALLAX = \\frac{8.794}{3600 $VN_HELIOCENTRIC_RADIUS}
```
where ``$VN_HELIOCENTRIC_RADIUS`` corresponds to the [`heliocentric_radius`](@ref), in 
astronomical units.
"""
function sun_equatorial_horizontal_parallax(heliocentric_radius)
    return 8.794 / (3600.0 * heliocentric_radius)
end

"""
    sun_ascension_parallax(
        radial_distance_equatorial_plane, 
        sun_equatorial_horizontal_parallax, 
        geocentric_sun_declination, 
        observer_local_hour
    )

The correction to the Sun’s right ascension due to the difference between geocentric and 
topocentric perspectives. 

```math
    $VN_SUN_ASCENSION_PARALLAX = \\text{atan2}(-$VN_RADIAL_DISTANCE_FROM_EQUATORIAL_PLANE 
    \\sin $VN_SUN_EQUATORIAL_HORIZONTAL_PARALLAX \\sin $VN_OBSERVER_LOCAL_HOUR, 
    \\cos $VN_GEOCENTRIC_SUN_DECLINATION - $VN_RADIAL_DISTANCE_FROM_EQUATORIAL_PLANE 
    \\sin $VN_SUN_EQUATORIAL_HORIZONTAL_PARALLAX \\cos $VN_OBSERVER_LOCAL_HOUR)
```
where:
- ``$VN_RADIAL_DISTANCE_FROM_EQUATORIAL_PLANE`` corresponds to 
    [`radial_distance_equatorial_plane`](@ref)) [m]
- ``$VN_SUN_EQUATORIAL_HORIZONTAL_PARALLAX`` corresponds to 
    [`sun_equatorial_horizontal_parallax`](@ref)
- ``$VN_GEOCENTRIC_SUN_DECLINATION`` corresponds to 
    [`geocentric_sun_declination`](@ref) 
- ``$VN_OBSERVER_LOCAL_HOUR`` corresponds to [`observer_local_hour`](@ref)
"""
function sun_ascension_parallax(
    radial_distance_equatorial_plane, 
    sun_equatorial_horizontal_parallax, 
    geocentric_sun_declination, 
    observer_local_hour
)
    s = -radial_distance_equatorial_plane * 
         sind(sun_equatorial_horizontal_parallax) * 
         sind(observer_local_hour)
    c = cosd(geocentric_sun_declination) - 
        radial_distance_equatorial_plane * 
        sind(sun_equatorial_horizontal_parallax) * 
        cosd(observer_local_hour)
    return atand(s, c)
end

"""
    topocentric_sun_ascension(geocentric_sun_ascension, sun_ascension_parallax)

The Sun’s right ascension as viewed from the observer’s actual location on Earth’s surface. 
 
```math
$VN_TOPOCENTRIC_SUN_ASCENSION = $VN_GEOCENTRIC_SUN_ASCENSION + $VN_SUN_ASCENSION_PARALLAX
```
where:
- ``$VN_GEOCENTRIC_SUN_ASCENSION`` corresponds to [`geocentric_sun_ascension`](@ref)
- ``$VN_SUN_ASCENSION_PARALLAX`` corresponds to [`sun_ascension_parallax`](@ref)
"""
function topocentric_sun_ascension(
    geocentric_sun_ascension, 
    sun_ascension_parallax
)
    return geocentric_sun_ascension + sun_ascension_parallax
end

"""
    topocentric_local_hour(observer_local_hour, sun_ascension_parallax)

The hour angle of a celestial object as seen from the observer’s exact location on Earth’s 
surface rather than its center. 

```math
$VN_TOPOCENTRIC_LOCAL_HOUR = $VN_OBSERVER_LOCAL_HOUR - $VN_SUN_ASCENSION_PARALLAX
``` 
where:
- ``$VN_OBSERVER_LOCAL_HOUR`` corresponds to the [`observer_local_hour`](@ref)
- ``$VN_SUN_ASCENSION_PARALLAX`` corresponds to the [`sun_ascension_parallax`](@ref)
"""
function topocentric_local_hour(observer_local_hour, sun_ascension_parallax)
    return observer_local_hour - sun_ascension_parallax
end

"""
    topocentric_sun_declination(
        radial_distance_equatorial_plane, 
        radial_distance_rotational_axis, 
        sun_equatorial_horizontal_parallax, 
        geocentric_sun_declination, 
        observer_local_hour, 
        sun_ascension_parallax
    )

The Sun’s declination as viewed from the observer’s actual location on Earth’s surface. 

```math
    $VN_TOPOCENTRIC_SUN_DECLINATION = \\text{atan2}((\\sin $VN_GEOCENTRIC_SUN_DECLINATION - 
    $VN_RADIAL_DISTANCE_FROM_ROTATIONAL_AXIS\\sin $VN_SUN_EQUATORIAL_HORIZONTAL_PARALLAX)
    \\cos $VN_SUN_ASCENSION_PARALLAX, \\cos $VN_GEOCENTRIC_SUN_DECLINATION - 
    $VN_RADIAL_DISTANCE_FROM_EQUATORIAL_PLANE\\sin $VN_SUN_EQUATORIAL_HORIZONTAL_PARALLAX 
    \\cos $VN_OBSERVER_LOCAL_HOUR)
```
where:
- ``$VN_RADIAL_DISTANCE_FROM_EQUATORIAL_PLANE`` corresponds to 
    [`radial_distance_equatorial_plane`](@ref) [m]
- ``$VN_RADIAL_DISTANCE_FROM_ROTATIONAL_AXIS`` corresponds to 
    [`radial_distance_rotational_axis`](@ref) [m]
- ``$VN_SUN_EQUATORIAL_HORIZONTAL_PARALLAX`` corresponds to 
    [`sun_equatorial_horizontal_parallax`](@ref)
- ``$VN_GEOCENTRIC_SUN_DECLINATION`` corresponds to [`geocentric_sun_declination`](@ref)
- ``$VN_OBSERVER_LOCAL_HOUR`` corresponds to [`observer_local_hour`](@ref)
- ``$VN_SUN_ASCENSION_PARALLAX`` corresponds to [`sun_ascension_parallax`](@ref)
"""
function topocentric_sun_declination(
    radial_distance_equatorial_plane, 
    radial_distance_rotational_axis, 
    sun_equatorial_horizontal_parallax, 
    geocentric_sun_declination, 
    observer_local_hour, 
    sun_ascension_parallax
)
    s = cosd(sun_ascension_parallax) * (
            sind(geocentric_sun_declination) - 
            radial_distance_rotational_axis * sind(sun_equatorial_horizontal_parallax)
        )
    c = cosd(geocentric_sun_declination) - 
        radial_distance_equatorial_plane * 
        sind(sun_equatorial_horizontal_parallax) * 
        cosd(observer_local_hour)
    return atand(s, c)
end

"""
    topocentric_apparent_elevation(
        observer_latitude,
        topocentric_sun_declination,
        topocentric_local_hour
    )

Topocentric elevation angle without atmospheric refraction correction. 

```math 
    $VN_TOPOCENTRIC_APPARENT_ELEVATION = \\text{asin}(\\sin $VN_OBSERVER_LATITUDE 
    \\sin $VN_TOPOCENTRIC_SUN_DECLINATION + \\cos $VN_OBSERVER_LATITUDE 
    \\cos $VN_TOPOCENTRIC_SUN_DECLINATION \\cos $VN_TOPOCENTRIC_LOCAL_HOUR)
```
where:
- ``$VN_OBSERVER_LATITUDE`` corresponds to `observer_latitude`
- ``$VN_TOPOCENTRIC_SUN_DECLINATION`` corresponds to [`topocentric_sun_declination`](@ref)
- ``$VN_TOPOCENTRIC_LOCAL_HOUR`` corresponds to [`topocentric_local_hour`](@ref)
"""
function topocentric_apparent_elevation(
    observer_latitude,
    topocentric_sun_declination,
    topocentric_local_hour
)
    s = sind(observer_latitude) * sind(topocentric_sun_declination)
    c = cosd(observer_latitude) * cosd(topocentric_sun_declination) * 
        cosd(topocentric_local_hour)
    return asind(s + c)
end

"""
    topocentric_elevation_correction(
        temperature,
        pressure,
        topocentric_apparent_elevation
    )

The adjustment applied to an object’s observed altitude to account for the bending of its 
light as it passes through Earth’s atmosphere. 
```math
    $VN_TOPOCENTRIC_ELEVATION_CORRECTION = \\frac{$VN_PRESSURE}{1010} \\frac{283}{273 + 
    $VN_TEMPERATURE} \\frac{1.02}{60\\tan($VN_TOPOCENTRIC_APPARENT_ELEVATION + 
    \\frac{10.3}{$VN_TOPOCENTRIC_APPARENT_ELEVATION + 5.11})}
```
where:
- ``$VN_TEMPERATURE`` corresponds to the local temperature [C]
- ``$VN_PRESSURE`` corresponds to the local pressure [mbar]
- ``$VN_TOPOCENTRIC_APPARENT_ELEVATION`` corresponds to the 
    [`topocentric_apparent_elevation`](@ref)
"""
function topocentric_elevation_correction(
    temperature,
    pressure,
    topocentric_apparent_elevation
)
    tae = topocentric_apparent_elevation
    return 283.0 * 1.02pressure / (
                1010.0(273 + temperature) * (60tand(tae + 10.3 / (tae + 5.11)))
           )
end

"""
    topocentric_elevation(
        topocentric_apparent_elevation, 
        topocentric_elevation_correction
    )

The angle between the Sun and the observer’s local horizon, measured at the observer’s 
location. 
```math
$VN_TOPOCENTRIC_ELEVATION = $VN_TOPOCENTRIC_APPARENT_ELEVATION + 
$VN_TOPOCENTRIC_ELEVATION_CORRECTION
```
where:
- ``$VN_TOPOCENTRIC_ELEVATION_CORRECTION`` corresponds to the 
    [`topocentric_elevation_correction`](@ref)
- ``$VN_TOPOCENTRIC_APPARENT_ELEVATION`` corresponds to the 
    [`topocentric_apparent_elevation`](@ref)
"""
function topocentric_elevation(
    topocentric_apparent_elevation, 
    topocentric_elevation_correction
)
    return topocentric_elevation_correction + topocentric_apparent_elevation
end

"""
    topocentric_astronomical_azimuth(
        observer_latitude, 
        topocentric_sun_declination, 
        topocentric_local_hour
    )

The Sun’s azimuth measured from true north eastward as seen from the observer’s location, 
based on astronomical (not navigational) convention. 

```math
    $VN_TOPOCENTRIC_ASTRONOMICAL_AZIMUTH = \\text{atan2}(\\sin $VN_TOPOCENTRIC_LOCAL_HOUR, 
    \\cos $VN_OBSERVER_LOCAL_HOUR \\sin $VN_OBSERVER_LATITUDE - 
    \\tan $VN_TOPOCENTRIC_SUN_DECLINATION \\cos $VN_OBSERVER_LATITUDE) \\mod 360
```
where:
- ``$VN_OBSERVER_LATITUDE`` corresponds to `observer_latitude`
- ``$VN_TOPOCENTRIC_SUN_DECLINATION`` corresponds to [`topocentric_sun_declination`](@ref)
- ``$VN_TOPOCENTRIC_LOCAL_HOUR`` corresponds to [`topocentric_local_hour`](@ref)

The ``\\mod(\\cdot)`` operator ensures that the results remains in the `[0, 360]` range.
"""
function topocentric_astronomical_azimuth(
    observer_latitude, 
    topocentric_sun_declination, 
    topocentric_local_hour
)
    s = sind(topocentric_local_hour)
    c = cosd(topocentric_local_hour) * sind(observer_latitude) - 
        tand(topocentric_sun_declination) * cosd(observer_latitude)
    return mod360(atand(s, c))
end

""" 
    topocentric_azimuth(topocentric_astronomical_azimuth)

The Sun’s azimuth from the observer’s location relative to a defined reference direction, 
typically true north, on the horizon. 

```math
$VN_TOPOCENTRIC_AZIMUTH = $VN_TOPOCENTRIC_ASTRONOMICAL_AZIMUTH + 180 \\mod 360
``` 
where ``$VN_TOPOCENTRIC_ASTRONOMICAL_AZIMUTH`` is the 
[`topocentric_astronomical_azimuth`](@ref).
"""
function topocentric_azimuth(topocentric_astronomical_azimuth)
    return mod360(topocentric_astronomical_azimuth + 180)
end
