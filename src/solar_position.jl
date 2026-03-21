using Pkg.Artifacts, Serialization, Dates

const EARTH_PERIODIC_TERMS = open(artifact"data/earth-periodic-terms.jld", "r") do io
    deserialize(io)
end

"""
    SolarPosition(azimuth, elevation, apparent_elevation)

Represents the position of the Sun relative to an observer.

# Fields
- `azimuth::Real`: solar azimuth angle (degrees, measured clockwise from north)
- `elevation::Real`: solar elevation angle above the horizon (degrees)
- `apparent_elevation::Real`: elevation corrected for atmospheric refraction (degrees)
"""
struct SolarPosition{T<:Real}
    azimuth::T
    elevation::T
    apparent_elevation::T
end

"""
    apparent_zenith(solpos::SolarPosition)

The complementary to 90° (degrees) of the Sun's apparent elevation.
"""
apparent_zenith(solpos::SolarPosition) = 90.0 - solpos.apparent_elevation

"""
    zenith(solpos::SolarPosition)

The complementary to 90° (degrees) of the Sun's elevation.
"""
zenith(solpos::SolarPosition) = 90.0 - solpos.elevation

"""
    sunray(solpos::SolarPosition)

The unit vector defining the direction of the Sun's rays, with respect to an observer on 
the Earth's surface. Conventionally, it points from the observer towards the Sun 
[parkin2010solar](@cite). 
"""
sunray(azimuth, apparent_elevation) = [
    sind(solpos.azimuth) * cosd(solpos.apparent_elevation)
    cosd(solpos.azimuth) * cosd(solpos.apparent_elevation)
    sind(solpos.apparent_elevation)
] 

# Limits an angle (in degrees) between 0° and 360°
mod360(angle) = mod(angle, 360)

"""
    delta_T(month::Int, year::Int)

Computes the difference between Terrestrial Dynamical Time (TD) and Universal Time (UT).

It is evaluated as a piecewise polynomial, as per 
[NASA](https://eclipse.gsfc.nasa.gov/SEcat5/deltatpoly.html), which in turn refer to 
[morrison2004historical, stephenson1986atlas](@cite). Note that ``\\Delta T`` is unknown 
for years before `-1999` and after `3000`. Values could be calculated for such intervals, 
although they are not intended to be used for these years.
"""
function delta_T(month::Int, year::Int)
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
""" 
julian_ephemeris_day(julian_day, ΔT) = julian_day + ΔT / 86400.0

"""
    julian_century(julian_day)

A time interval of exactly 36,525 days (365.25 days × 100) used as a standard unit in 
astronomy.
"""
julian_century(julian_day) = (julian_day - 2451545.0) / 36525.0

"""
    julian_ephemeris_century(julian_ephemeris_day)

A 36,525-day interval measured specifically in Ephemeris Time (or its modern dynamical time
scales) for high-precision astronomical modeling. 
"""
julian_ephemeris_century(julian_ephemeris_day) = (julian_ephemeris_day − 2451545.0) / 36525.0

"""
    julian_ephemeris_millenium(julian_ephemeris_century)

A 1,000-year interval equal to 365,250 ephemeris days, defined within Ephemeris Time for 
long-term astronomical calculations. 
"""
julian_ephemeris_millenium(julian_ephemeris_century) = julian_ephemeris_century / 10.0

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
"""
function heliocentric_longitude(julian_ephemeris_millenium)
    μ = heliocentric_polynomial(julian_ephemeris_millenium, EARTH_PERIODIC_TERMS.M)
    return μ |> rad2deg |> mod360
end

"""
    heliocentric_latitude(julian_ephemeris_millenium)

A celestial object's angular position around the Sun measured in the ecliptic plane from a 
defined reference direction.
"""
function heliocentric_latitude(julian_ephemeris_millenium)
    λ = heliocentric_polynomial(julian_ephemeris_millenium, EARTH_PERIODIC_TERMS.L)
    return λ |> rad2deg
end

"""
    heliocentric_radius(julian_ephemeris_millenium)

The distance from the Sun to the celestial object in a chosen heliocentric coordinate 
system. In this case, it corresponds to the Earth-Sun distance.
"""
function heliocentric_radius(julian_ephemeris_millenium)
    return heliocentric_polynomial(julian_ephemeris_millenium, EARTH_PERIODIC_TERMS.R)
end

"""
    geocentric_longitude(heliocentric_longitude)

A celestial object's ecliptic longitude as measured from Earth’s center, referenced to the 
ecliptic plane and the vernal equinox. 
"""
function geocentric_longitude(heliocentric_longitude)
    return mod360(heliocentric_longitude + 180)
end

"""
    geocentric_latitude(heliocentric_latitude)

A celestial object's angular distance north or south of the ecliptic plane as measured from 
Earth’s center. 
"""
geocentric_latitude(heliocentric_latitude) = -heliocentric_latitude

"""
    geocentric_sun_ascension(
        apparent_sun_longitude,
        elliptic_obliquity,
        geocentric_latitude
    )

The angle measured eastward from the geocentric meridian, that locates a direction or point 
relative to Earth’s ellipsoid rather than its rotational axis. 
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
"""
function mean_moon_elongation_from_sun(julian_ephemeris_century)
    jec = julian_ephemeris_century
    return 297.85036 + 445267.111480jec − 0.0019142jec ^ 2 + jec^3 / 189474.0
end

"""
    mean_sun_anomaly(julian_ephemeris_century)

The angular position of the Sun in its elliptical orbit, measured from perihelion and 
increasing uniformly in time. 
"""
function mean_sun_anomaly(julian_ephemeris_century)
    jec = julian_ephemeris_century
    return 357.52772 + 35999.050340jec − 0.0001603jec^2 − jec^3 / 300000.0
end

"""
    mean_moon_anomaly(julian_ephemeris_century)

The uniformly increasing angular position of the Moon in its elliptical orbit measured from 
perigee.
"""
function mean_moon_anomaly(julian_ephemeris_century)
    jec = julian_ephemeris_century
    return 134.96298 + 477198.867398jec + 0.0086972jec^2 + jec^3 / 56250.0
end

"""
    moon_latitude_argument(julian_ephemeris_century)

The angle from the Moon’s ascending node to its position measured along its orbit, using 
mean (unperturbed) orbital elements. 
"""
function moon_latitude_argument(julian_ephemeris_century)
    jec = julian_ephemeris_century
    return 93.27191 + 483202.017538jec - 0.0036825jec^2 + jec^3 / 327270.0
end

"""
    ascending_moon_longitude(julian_ephemeris_century)

The ecliptic longitude referenced to the mean equinox of the date, of the point where the 
Moon’s mean orbit crosses northward through the ecliptic. 
"""
function ascending_moon_longitude(julian_ephemeris_century)
    jec = julian_ephemeris_century
    return 125.04452 - 1934.136261jec + 0.0020708jec^2 + jec^3 / 450000.0
end

"""
    nutation_coefficients(julian_ephemeris_century)

Constructs the vector of weights for evaluating the [`nutation_longitude`](@ref) and 
[`nutation_obliquity`](@ref). 
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
    nutation_longitude(julian_ephemeris_century, nutation_coefficients)

The small periodic variation in Earth’s ecliptic longitude caused by gravitational torques 
from the Moon and Sun. 
"""
function nutation_longitude(julian_ephemeris_century, nutation_coefficients)
    jec = julian_ephemeris_century
    Ψ = EARTH_PERIODIC_TERMS.Ψ
    W = nutation_coefficients
    return sum((Ψ[:, 1] .+ Ψ[:, 2] * jec) .* sind.(W)) / 3.6e7
end

"""
    nutation_obliquity(julian_ephemeris_century, nutation_coefficients)

The small periodic variation in Earth’s axial tilt (obliquity) resulting from 
gravitational perturbations. 
"""
function nutation_obliquity(julian_ephemeris_century, nutation_coefficients)
    jec = julian_ephemeris_century
    E = EARTH_PERIODIC_TERMS.E
    W = nutation_coefficients
    return sum((E[:, 1] .+ E[:, 2] * jec) .* cosd.(W)) / 3.6e7
end

"""
    mean_elliptic_obliquity(julian_ephemeris_millenium)

The angle between Earth’s equator and the mean (long-term averaged) ecliptic defined by an 
unperturbed elliptical Earth orbit. 
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
"""
function elliptic_obliquity(mean_elliptic_obliquity, nutation_obliquity)
    return mean_elliptic_obliquity / 3600.0 + nutation_obliquity
end

"""
    aberration_correction(heliocentric_radius)

An adjustment applied to celestial coordinates to account for the apparent displacement 
caused by Earth’s motion through space. 
"""
aberration_correction(heliocentric_radius) = -20.4898 / (3600.0heliocentric_radius)

"""
    apparent_sun_longitude(geocentric_longitude, nutation_longitude, aberration_correction)

The Sun’s ecliptic longitude after applying corrections for nutation and aberration.
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
"""
reduced_observer_latitude(observer_latitude) = atand(0.99664719 * tand(observer_latitude))

"""
    radial_distance_equatorial_plane(
        observer_latitude,
        observer_altitude,
        reduced_observer_latitude
    )

The perpendicular distance of a point from Earth’s equatorial plane. 
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
    topocentric_elevation(topocentric_apparent_elevation, topocentric_elevation_correction)

The angle between the Sun and the observer’s local horizon, measured at the observer’s 
location. 
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
"""
function topocentric_azimuth(topocentric_astronomical_azimuth)
    return mod360(topocentric_astronomical_azimuth + 180)
end

"""
    spa(location::Location, datetime::DateTime)

Computes the solar position in terms of Solar azimuth, elevation, and apparent elevation, 
following the solar position algorithm [reda2004solar](@cite).
"""
function spa(location::Location, datetime::DateTime)
    ΔT = delta_T(Dates.month(datetime), Dates.year(datetime))

    jd = Dates.datetime2julian(datetime)
    jed = julian_ephemeris_day(jd, ΔT)
    jc = julian_century(jd)
    jec = julian_ephemeris_century(jed)
    jem = julian_ephemeris_millenium(jec)

    hel_lat = heliocentric_latitude(jem)
    hel_lon = heliocentric_longitude(jem)
    hel_radius = heliocentric_radius(jem)

    geo_lon = geocentric_longitude(hel_lon)
    geo_lat = geocentric_latitude(hel_lat)

    nut_coeff = nutation_coefficients(jec)
    nut_lon = nutation_longitude(jec, nut_coeff)
    nut_obl = nutation_obliquity(jec, nut_coeff)

    mean_ell_obl = mean_elliptic_obliquity(jem)
    ell_obl = elliptic_obliquity(mean_ell_obl, nut_obl)
    aberr_corr = aberration_correction(hel_radius)

    app_sun_lon = apparent_sun_longitude(geo_lon, nut_lon, aberr_corr)
    geo_sun_asc = geocentric_sun_ascension(app_sun_lon, ell_obl, geo_lat)
    geo_sun_dec = geocentric_sun_declination(app_sun_lon, ell_obl, geo_lat)

    mean_sun_lon = mean_sun_longitude(jem)
    mean_sid_gw_time = mean_sidereal_greenwich_time(jd, jc)
    app_sid_gw_time = apparent_sidereal_greenwich_time(mean_sid_gw_time, nut_lon, ell_obl)

    obs_local_hr = observer_local_hour(location.longitude, app_sid_gw_time, geo_sun_asc)
    red_obs_lat = reduced_observer_latitude(location.latitude)
    rad_dist_eq_plane = radial_distance_equatorial_plane(
        location.latitude,
        location.altitude,
        red_obs_lat,
    )

    rad_dist_rot_axis = radial_distance_rotational_axis(
        location.latitude,
        location.altitude,
        red_obs_lat,
    )

    sun_eq_horiz_px = sun_equatorial_horizontal_parallax(hel_radius)
    sun_asc_px = sun_ascension_parallax(
        rad_dist_eq_plane, 
        sun_eq_horiz_px,
        geo_sun_dec,
        obs_local_hr,
    )

    top_sun_asc = topocentric_sun_ascension(geo_sun_asc, sun_asc_px)
    top_local_hr = topocentric_sun_ascension(obs_local_hr, sun_asc_px)
    top_sun_dec = topocentric_sun_declination(
        rad_dist_eq_plane,
        rad_dist_rot_axis,
        sun_eq_horiz_px,
        geo_sun_dec,
        obs_local_hr,
        sun_asc_px,
    )

    top_app_elev = topocentric_apparent_elevation(
        location.latitude, 
        top_sun_dec, 
        top_local_hr,
    )

    top_elev_corr = topocentric_elevation_correction(
        location.temperature,
        location.pressure,
        top_app_elev,
    )

    top_elev = topocentric_elevation(top_app_elev, top_elev_corr)
    top_astr_azimuth = topocentric_astronomical_azimuth(
        location.latitude,
        top_sun_dec,
        top_local_hr,
    )
    top_azimuth = topocentric_azimuth(top_astr_azimuth)

    top_azimuth, top_elev, top_app_elev = promote(top_azimuth, top_elev, top_app_elev)
    return SolarPosition(top_azimuth, top_elev, top_app_elev)
end
