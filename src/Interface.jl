module Interface

using Dates, SolarFunctions.Model, FunctionFlow

# Solar position algorithm
delta_T = Node(Model.delta_T, :datetime)

julian_day = Node(datetime2julian, :datetime)

julian_ephemeris_day = Node(
    Model.julian_ephemeris_day,
    julian_day,
    Ambiguity(:delta_T, delta_T; default=delta_T)
)
julian_century = Node(Model.julian_century, julian_day)

julian_ephemeris_century = Node(Model.julian_ephemeris_century, julian_ephemeris_day)

julian_ephemeris_millenium = Node(Model.julian_ephemeris_millenium, julian_ephemeris_century)

heliocentric_latitude = Node(Model.heliocentric_latitude, julian_ephemeris_millenium)

heliocentric_longitude = Node(Model.heliocentric_longitude, julian_ephemeris_millenium)

heliocentric_radius = Node(Model.heliocentric_radius, julian_ephemeris_millenium)

geocentric_longitude = Node(Model.geocentric_longitude, heliocentric_longitude)

geocentric_latitude = Node(Model.geocentric_latitude, heliocentric_latitude)

nutation_coefficients = Node(Model.nutation_coefficients, julian_ephemeris_century)

nutation_longitude = Node(Model.nutation_longitude, julian_ephemeris_century, nutation_coefficients)

nutation_obliquity = Node(Model.nutation_obliquity, julian_ephemeris_century, nutation_coefficients)

mean_elliptic_obliquity = Node(Model.mean_elliptic_obliquity, julian_ephemeris_millenium)

elliptic_obliquity = Node(Model.elliptic_obliquity, mean_elliptic_obliquity, nutation_obliquity)

aberration_correction = Node(Model.aberration_correction, heliocentric_radius)

apparent_sun_longitude = Node(
    Model.apparent_sun_longitude,
    geocentric_longitude,
    nutation_longitude,
    aberration_correction
)

geocentric_sun_ascension = Node(
    Model.geocentric_sun_ascension,
    apparent_sun_longitude,
    elliptic_obliquity,
    geocentric_latitude,
)

geocentric_sun_declination = Node(
    Model.geocentric_sun_declination,
    apparent_sun_longitude,
    elliptic_obliquity,
    geocentric_latitude,
)

mean_sun_longitude = Node(Model.mean_sun_longitude, julian_ephemeris_millenium)

mean_sidereal_greenwich_time = Node(
    Model.mean_sidereal_greenwich_time, 
    julian_day, 
    julian_century
)

apparent_sidereal_greenwich_time = Node(
    Model.apparent_sidereal_greenwich_time,
    mean_sidereal_greenwich_time,
    nutation_longitude,
    elliptic_obliquity
)

observer_local_hour = Node(
    Model.observer_local_hour,
    :observer_longitude,
    apparent_sidereal_greenwich_time,
    geocentric_sun_ascension
)

reduced_observer_latitude = Node(Model.reduced_observer_latitude, :observer_latitude)

radial_distance_equatorial_plane = Node(
    Model.radial_distance_equatorial_plane,
    :observer_latitude,
    :observer_altitude,
    reduced_observer_latitude
)

radial_distance_rotational_axis = Node(
    Model.radial_distance_rotational_axis,
    :observer_latitude,
    :observer_altitude,
    reduced_observer_latitude
)

sun_equatorial_horizontal_parallax = Node(
    Model.sun_equatorial_horizontal_parallax, 
    heliocentric_radius
)

sun_ascension_parallax = Node(
    Model.sun_ascension_parallax,
    radial_distance_equatorial_plane,
    sun_equatorial_horizontal_parallax,
    geocentric_sun_declination,
    observer_local_hour
)

topocentric_sun_ascension = Node(
    Model.topocentric_sun_ascension,
    geocentric_sun_ascension,
    sun_ascension_parallax
)

topocentric_local_hour = Node(
    Model.topocentric_sun_ascension,
    observer_local_hour,
    sun_ascension_parallax
)

topocentric_sun_declination = Node(
    Model.topocentric_sun_declination,
    radial_distance_equatorial_plane,
    radial_distance_rotational_axis,
    sun_equatorial_horizontal_parallax,
    geocentric_sun_declination,
    observer_local_hour,
    sun_ascension_parallax
)

topocentric_apparent_elevation = Node(
    Model.topocentric_apparent_elevation,
    :observer_latitude,
    topocentric_sun_declination,
    topocentric_local_hour
)

topocentric_elevation_correction = Node(
    Model.topocentric_elevation_correction,
    :temperature,
    :pressure,
    topocentric_apparent_elevation
)

topocentric_elevation = Node(
    Model.topocentric_elevation,
    topocentric_apparent_elevation,
    topocentric_elevation_correction
)

topocentric_astronomical_azimuth = Node(
    Model.topocentric_astronomical_azimuth,
    :observer_latitude,
    topocentric_sun_declination,
    topocentric_local_hour
)

topocentric_azimuth = Node(
    Model.topocentric_azimuth,
    topocentric_astronomical_azimuth
)

# Atmosphere
relative_airmass_simple = Node(Model.relative_airmass_simple, topocentric_elevation)

relative_airmass_kasten1966 = Node(
    Model.relative_airmass_kasten1966,
    topocentric_apparent_elevation
)

relative_airmass_youngirvine1967 = Node(
    Model.relative_airmass_youngirvine1967,
    topocentric_elevation
)

relative_arimass_kastenyoung1989 = Node(
    Model.relative_arimass_kastenyoung1989,
    topocentric_apparent_elevation
)

relative_airmass_gueymard1993 = Node(
    Model.relative_airmass_gueymard1993,
    topocentric_apparent_elevation
)

relative_airmass_young1994 = Node(
    Model.relative_airmass_young1994,
    topocentric_elevation
)

relative_airmass_pickering2002 = Node(
    Moddel.relative_airmass_pickering2002,
    topocentric_apparent_elevation
)

relative_airmass_gueymard2003 = Node(
    Model.relative_airmass_gueymard2003,
    topocentric_apparent_elevation
)

absolute_airmass = Node(
    Ambiguity(
        relative_airmass_simple,
        relative_airmass_kasten1966,
        relative_airmass_youngirvine1967,
        relative_arimass_kastenyoung1989,
        relative_airmass_gueymard1993,
        relative_airmass_young1994,
        relative_airmass_pickering2002,
        relative_airmass_gueymard2003;
        default=relative_arimass_kastenyoung1989
    ),
    :pressure
)

# Clearsky
clearsky_ineichen = Node(
    Model.clearsky_ineichen,
    topocentric_apparent_elevation,
    :observer_altitude,
    absolute_airmass,
    linke_turbidity,
    extraterrestial_radiation,
    :perez_enhancement
)

end # module Interface
