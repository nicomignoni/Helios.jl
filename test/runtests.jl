using Test, Dates, Helios

function print_error(name, test, correct)
    println(
        "$(name):
        - test: $(test)
        - correct: $(correct)
        - error: $(100(test - correct)/correct |> abs) %\n"
    )
end

@testset "solar position" begin
    # Test data
    latitude = 39.742476
    longitude = -105.1786
    altitude = 1830.14
    temperature = 11.0
    pressure = 820
    tz = -7
    datetime = DateTime(2003, 10, 17, 12 - tz, 30, 30)
    ΔT = 67

    correct = (
        julian_day=2452930.312847,
        julian_century=0.037928,
        heliocentric_longitude=24.0182616917,
        heliocentric_latitude=-0.0001011219,
        heliocentric_radius=0.9965422974,
        geocentric_longitude=204.0182616917,
        geocentric_latitude=0.0001011219,
        geocentric_sun_ascension=202.22741,
        geocentric_sun_declination=-9.31434,
        nutation_longitude=-0.00399840,
        nutation_obliquity=0.00166657,
        elliptic_obliquity=23.440465,
        apparent_sun_longitude=204.0085519281,
        observer_local_hour=11.105900,
        topocentric_sun_ascension=202.22704,
        topicentric_local_hour=11.10629,
        topocentric_sun_declination=-9.316179,
        topocentric_azimuth=194.34024,
        topocentric_elevation=39.88838
    )

    jd = Helios.Dates.datetime2julian(datetime)
    jed = Helios.julian_ephemeris_day(jd, ΔT)
    jc = Helios.julian_century(jd)
    jec = Helios.julian_ephemeris_century(jed)
    jem = Helios.julian_ephemeris_millenium(jec)

    hel_lat = Helios.heliocentric_latitude(jem)
    hel_lon = Helios.heliocentric_longitude(jem)
    hel_radius = Helios.heliocentric_radius(jem)

    geo_lon = Helios.geocentric_longitude(hel_lon)
    geo_lat = Helios.geocentric_latitude(hel_lat)

    nut_coeff = Helios.nutation_coefficients(jec)
    nut_lon = Helios.nutation_longitude(jec, nut_coeff)
    nut_obl = Helios.nutation_obliquity(jec, nut_coeff)

    mean_ell_obl = Helios.mean_elliptic_obliquity(jem)
    ell_obl = Helios.elliptic_obliquity(mean_ell_obl, nut_obl)
    aberr_corr = Helios.aberration_correction(hel_radius)

    app_sun_lon = Helios.apparent_sun_longitude(geo_lon, nut_lon, aberr_corr)
    geo_sun_asc = Helios.geocentric_sun_ascension(app_sun_lon, ell_obl, geo_lat)
    geo_sun_dec = Helios.geocentric_sun_declination(app_sun_lon, ell_obl, geo_lat)

    mean_sun_lon = Helios.mean_sun_longitude(jem)
    mean_sid_gw_time = Helios.mean_sidereal_greenwich_time(jd, jc)
    app_sid_gw_time = Helios.apparent_sidereal_greenwich_time(
        mean_sid_gw_time, nut_lon, ell_obl
    )

    obs_local_hr = Helios.observer_local_hour(
        longitude, app_sid_gw_time, geo_sun_asc
    )
    red_obs_lat = Helios.reduced_observer_latitude(latitude)
    rad_dist_eq_plane = Helios.radial_distance_equatorial_plane(
        latitude, altitude, red_obs_lat,
    )

    rad_dist_rot_axis = Helios.radial_distance_rotational_axis(
        latitude, altitude, red_obs_lat
    )

    sun_eq_horiz_px = Helios.sun_equatorial_horizontal_parallax(hel_radius)
    sun_asc_px = Helios.sun_ascension_parallax(
        rad_dist_eq_plane, sun_eq_horiz_px, geo_sun_dec, obs_local_hr
    )

    top_sun_asc = Helios.topocentric_sun_ascension(geo_sun_asc, sun_asc_px)
    top_local_hr = Helios.topocentric_sun_ascension(obs_local_hr, sun_asc_px)
    top_sun_dec = Helios.topocentric_sun_declination(
        rad_dist_eq_plane,
        rad_dist_rot_axis,
        sun_eq_horiz_px,
        geo_sun_dec,
        obs_local_hr,
        sun_asc_px,
    )

    top_app_elev = Helios.topocentric_apparent_elevation(
        latitude, top_sun_dec, top_local_hr
    )

    top_elev_corr = Helios.topocentric_elevation_correction(
        temperature, pressure, top_app_elev
    )

    top_elev = Helios.topocentric_elevation(top_app_elev, top_elev_corr)
    top_astr_azimuth = Helios.topocentric_astronomical_azimuth(
        latitude, top_sun_dec, top_local_hr
    )
    top_azimuth = Helios.topocentric_azimuth(top_astr_azimuth)

    print_error("Julian day", jd, correct.julian_day)
    print_error("Julian century", jc, correct.julian_century)
    print_error("Heliocentric longitude", hel_lon, correct.heliocentric_longitude)
    print_error("Heliocentric latitude", hel_lat, correct.heliocentric_latitude)
    print_error("Heliocentric radius", hel_radius, correct.heliocentric_radius)
    print_error("Geocentric longitude", geo_lon, correct.geocentric_longitude)
    print_error("Geocentric latitude", geo_lat, correct.geocentric_latitude)
    print_error("Nutation longitude", nut_lon, correct.nutation_longitude)
    print_error("Nutation obliquity", nut_obl, correct.nutation_obliquity)
    print_error("Elliptic obliquity", ell_obl, correct.elliptic_obliquity)
    print_error("Apparent Sun longitude", app_sun_lon, correct.apparent_sun_longitude)
    print_error("Geocentric Sun ascension", geo_sun_asc, correct.geocentric_sun_ascension)
    print_error("Geocentric Sun declination", geo_sun_dec, correct.geocentric_sun_declination)
    print_error("Obsever hour angle", obs_local_hr, correct.observer_local_hour)
    print_error("Topocentric Sun declination", top_sun_dec, correct.topocentric_sun_declination)
    print_error("Topocentric Sun ascension", top_sun_asc, correct.topocentric_sun_ascension)
    print_error("Topocentric local hour", top_local_hr, correct.topicentric_local_hour)
    print_error("Topocentric elevation", top_elev, correct.topocentric_elevation)
    print_error("Topocentric azimuth", top_azimuth, correct.topocentric_azimuth)
end

@testset "clearsky" begin
    location = Location(0.0, 0.0, 0.0)
    solpos = SolarPosition(0.0, 0.0, 30.4)

    correct = (
        ineichen=(dni=1007.606890209482, dhi=42.47213240572813, ghi=552.3552398128525),
        haurwitz=(ghi=494.477044,)
    )

    # Ineichen
    irradiance = clearsky_ineichen(
        location,
        now(); # not really needed, just a slack argument for this test
        solpos=solpos,
        relative_airmass = Helios.ATMOSPHERIC_PRESSURE / location.pressure,
        linke_turbidity = 2.1,
        extraterrestial_radiation = 1364.0,
        perez_enhancement = false
    )
    print_error("Ineichen (DNI)", irradiance.dni, correct.ineichen.dni)
    print_error("Ineichen (DHI)", irradiance.dhi, correct.ineichen.dhi)
    print_error("Ineichen (GHI)", irradiance.ghi, correct.ineichen.ghi)

    # Haurwitz
    ghi = clearsky_haurwitz(solpos)
    print_error("Haurwitz (GHI)", ghi, correct.haurwitz.ghi)
end
