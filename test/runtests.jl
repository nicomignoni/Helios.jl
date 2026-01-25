using Test, Dates, SolarFunctions.Interface

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

    print_error(
        "Julian day",
        Interface.julian_day(datetime=datetime),
        correct.julian_day
    )
    print_error(
        "Julian century",
        Interface.julian_ephemeris_day(:delta_T, datetime=datetime, delta_T=ΔT),
        correct.julian_century
    )

    print_error(
        "Heliocentric longitude",
        Interface.heliocentric_longitude(:delta_T, datetime=datetime, delta_T=ΔT),
        correct.heliocentric_longitude
    )
    print_error(
        "Heliocentric latitude",
        Interface.heliocentric_latitude(:delta_T, datetime=datetime, delta_T=ΔT),
        correct.heliocentric_latitude
    )
    print_error(
        "Heliocentric radius",
        Interface.heliocentric_radius(:delta_T, datetime=datetime, delta_T=ΔT),
        correct.heliocentric_radius
    )
    print_error(
        "Geocentric longitude",
        Interface.geocentric_longitude(:delta_T, datetime=datetime, delta_T=ΔT),
        correct.geocentric_longitude
    )
    print_error(
        "Geocentric latitude",
        Interface.geocentric_latitude(:delta_T, datetime=datetime, delta_T=ΔT),
        correct.geocentric_latitude
    )
    print_error(
        "Nutation longitude",
        Interface.nutation_longitude(:delta_T, datetime=datetime, delta_T=ΔT),
        correct.nutation_longitude
    )
    print_error(
        "Nutation obliquity",
        Interface.nutation_obliquity(:delta_T, datetime=datetime, delta_T=ΔT),
        correct.nutation_obliquity
    )
    print_error(
        "Elliptic obliquity",
        Interface.elliptic_obliquity(:delta_T, datetime=datetime, delta_T=ΔT),
        correct.elliptic_obliquity
    )
    print_error(
        "Apparent Sun longitude",
        Interface.apparent_sun_longitude(:delta_T, datetime=datetime, delta_T=ΔT),
        correct.apparent_sun_longitude
    )
    print_error(
        "Geocentric Sun ascension",
        Interface.geocentric_sun_ascension(:delta_T, datetime=datetime, delta_T=ΔT),
        correct.geocentric_sun_ascension
    )
    print_error(
        "Geocentric Sun declination",
        Interface.geocentric_sun_declination(:delta_T, datetime=datetime, delta_T=ΔT),
        correct.geocentric_sun_declination
    )
    print_error(
        "Obsever hour angle",
        Interface.observer_local_hour(:delta_T, observer_longitude=longitude, datetime=datetime, delta_T=ΔT),
        correct.observer_local_hour
    )
    print_error(
        "Topocentric Sun declination",
        Interface.topocentric_sun_declination(
            :delta_T,
            observer_latitude=latitude,
            observer_altitude=altitude,
            observer_longitude=longitude,
            datetime=datetime,
            delta_T=ΔT
        ),
        correct.topocentric_sun_declination
    )
    print_error(
        "Topocentric Sun ascension",
        Interface.topocentric_sun_ascension(
            :delta_T,
            observer_latitude=latitude,
            observer_longitude=longitude,
            observer_altitude=altitude,
            datetime=datetime,
            delta_T=ΔT
        ),
        correct.topocentric_sun_ascension
    )
    print_error(
        "Topocentric hour angle",
        Interface.topocentric_local_hour(
            :delta_T,
            observer_longitude=longitude,
            observer_latitude=latitude,
            observer_altitude=altitude,
            datetime=datetime,
            delta_T=ΔT
        ),
        correct.topicentric_local_hour
    )
    print_error(
        "Topocentric elevation",
        Interface.topocentric_elevation(
            :delta_T,
            observer_longitude=longitude,
            observer_latitude=latitude,
            observer_altitude=altitude,
            temperature=temperature,
            pressure=pressure,
            datetime=datetime,
            delta_T=ΔT
        ),
        correct.topocentric_elevation
    )
    print_error(
        "Topocentric azimuth",
        Interface.topocentric_azimuth(
            :delta_T,
            observer_latitude=latitude,
            observer_longitude=longitude,
            observer_altitude=altitude,
            datetime=datetime,
            delta_T=ΔT
        ),
        correct.topocentric_azimuth
    )
end
#
# @testset "clearsky" begin
#     sun_apparent_elevation = 30.4
#     observer_altitude = 0.0
#     absolute_airmass = 1.0
#     linke_turbidity = 2.1
#     extraterrestial_radiation = 1364.0
#     perez_enanchement = false
#
#     correct = (
#         ineichen=(
#             dni=1007.606890209482,
#             dhi=42.47213240572813,
#             ghi=552.3552398128525,
#         ),
#         haurwitz=(
#             ghi=494.477044,
#         )
#     )
#
#     # Ineichen
#     dni, dhi, ghi = Model.clearsky_ineichen(
#         sun_apparent_elevation,
#         observer_altitude,
#         absolute_airmass,
#         linke_turbidity,
#         extraterrestial_radiation,
#         perez_enanchement
#     )
#     print_error("Ineichen (DNI)", dni, correct.ineichen.dni)
#     print_error("Ineichen (DHI)", dhi, correct.ineichen.dhi)
#     print_error("Ineichen (GHI)", ghi, correct.ineichen.ghi)
#
#     # Haurwitz
#     ghi = Model.clearsky_haurwitz(sun_apparent_elevation)
#     print_error("Haurwitz (GHI)", ghi, correct.haurwitz.ghi)
# end
