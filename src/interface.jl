using Dates

"""
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

function Base.getproperty(solpos::SolarPosition, s::Symbol)
    if s === :zenith
        return 90 - getfield(solpos, :elevation)
    elseif s === :apparent_zenith 
        return 90 - getfield(solpos, :apparent_elevation)
    else
        return getfield(solpos, s)
    end
end

Base.propertynames(solpos::SolarPosition, private=false) =
    (fieldnames(typeof(solpos))..., :zenith, :apparent_zenith)

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
