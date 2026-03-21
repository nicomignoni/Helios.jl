using Dates

const SOLAR_CONSTANT=1366.1 # W/m²

"""
    Irradiance(dni, dhi, ghi)

Container for solar irradiance components.

# Fields
- `dni`: Direct Normal Irradiance [W/m²].
- `dhi`: Diffuse Horizontal Irradiance [W/m²].
- `ghi`: Global Horizontal Irradiance [W/m²].
"""
struct Irradiance{T<:Real}
    dni::T
    dhi::T
    ghi::T
end

"""
    day_angle(datetime::DateTime, offset::Int)

Calculates the day angle [degrees] for the Earth's orbit around the Sun. 

For the Spencer method, `offset=1`; for the ASCE method, `offset=0`.
"""
function day_angle(datetime::DateTime, offset)
    return 360.0(dayofyear(datetime) - offset) / 365.0
end

"""
    extraterrestrial_irradiance_spencer1971(datetime::DateTime)

Calculates the extraterrestrial irradiance [W/m²], as proposed in 
[spencer1971fourier](@cite).
"""
function extraterrestrial_irradiance_spencer1971(datetime::DateTime)
    da = day_angle(datetime, 1.0)
    R = 1.00011 + 0.034221cosd(da) + 0.00128sind(da) + 0.000719cosd(2da) + 7.7e-05sind(2da)
    return R * SOLAR_CONSTANT
end

"""
    extraterrestrial_irradiance_asce(datetime::DateTime)

Calculates the extraterrestrial irradiance [W/m²] using the ASCE evapotranspiration 
equation [walter2000asce](@cite).
"""
function extraterrestrial_irradiance_asce(datetime::DateTime)
    da = day_angle(datetime, 0.0)
    R = 1 + 0.033cosd(da) 
    return R * SOLAR_CONSTANT
end

"""
    extraterrestrial_irradiance_nrel(heliocentric_radius)

Calculates the extraterrestrial irradiance [W/m²] as implemented by the NREL SPA 
[reda2004solar](@cite). 
"""
function extraterrestrial_irradiance_nrel(heliocentric_radius)
    return SOLAR_CONSTANT / heliocentric_radius^2
end

"""
    angle_of_incidence_projection(
        solpos::SolarPosition, surface_tilt, surface_roll, surface_azimuth
    )

Calculates the dot product of the sun position unit vector and the surface normal unit 
vector, i.e., the cosine of the angle of incidence.
"""
function angle_of_incidence_projection(
    solpos::SolarPosition, surface_tilt, surface_roll, surface_azimuth
)
    return cosd(surface_roll) * cosd(surface_tilt) * sind(solpos.apparent_elevation) +
           cosd(solpos.apparent_elevation) * sind(solpos.azimuth) * (
                sind(surface_tilt) * sind(surface_azimuth) + 
                cosd(surface_azimuth) * cosd(surface_tilt) * sind(surface_roll)
           ) +
           cosd(solpos.azimuth) * cosd(solpos.apparent_elevation) * (
                cosd(surface_azimuth) * sind(surface_tilt) - 
                sind(surface_azimuth) * cosd(surface_tilt) * sind(surface_roll)
           )
end

"""
    beam_component(solpos::SolarPosition, dni, surface_tilt, surface_roll, surface_azimuth)

Calculates the beam component of the plane of array irradiance.
"""
function beam_component(solpos::SolarPosition, dni, surface_tilt, surface_roll, surface_azimuth)
    aoi_proj = angle_of_incidence_projection(solpos, surface_tilt, surface_roll, surface_azimuth) 
    return max(0.0, aoi_proj) * dni
end

"""
    sky_diffuse_component_isotropic(dhi, surface_tilt, surface_roll)

Determine diffuse irradiance from the sky on a tilted surface using the isotropic sky 
model [loutzenhiser2007empirical, kamphuis2020perspectives](@cite).
"""
function sky_diffuse_component_isotropic(dhi, surface_tilt, surface_roll)
    return 0.5 * dhi * (1 + cosd(surface_tilt) * cosd(surface_roll)) 
end

"""
    sky_diffuse_component_klucher(
        solpos::SolarPosition, dhi, ghi, 
        surface_tilt, surface_roll, surface_azimuth
    )

Determine diffuse irradiance from the sky on a tilted surface using the Klucher (1979) 
model [klucher1979evaluation, loutzenhiser2007empirical](@cite).
"""
function sky_diffuse_component_klucher(
    solpos::SolarPosition, dhi, ghi, 
    surface_tilt, surface_roll, surface_azimuth
)
    F = 1 - (dhi / ghi)^2
    aoi_proj = angle_of_incidence_projection(solpos, surface_tilt, surface_roll, surface_azimuth)
    return 0.5(1 + cosd(surface_tilt) * cosd(surface_roll)) *
           (1 + F * (1 - cosd(surface_tilt) * cosd(surface_roll))^1.5) *
           (1 + F * max(0.0, aoi_proj)^2 * sind(solpos.zenith)^3)
end

"""
    sky_diffuse_component_haydavies(
        solpos::SolarPosition, dni, dni_extra, dhi, 
        surface_tilt, surface_roll, surface_azimuth
    )

Determine diffuse irradiance from the sky on a tilted surface using the Hay and Davies 
(1980) model [hay5972calculations, loutzenhiser2007empirical](@cite).
"""
function sky_diffuse_component_haydavies(
    solpos::SolarPosition, dni, dni_extra, dhi, 
    surface_tilt, surface_roll, surface_azimuth
)
    aoi_proj = angle_of_incidence_projection(solpos, surface_tilt, surface_roll, surface_azimuth)
    projection_ratio = max(0.0, aoi_proj) / max(sind(solpos.elevation), 0.01745)
    anisotropy_index = dni / dni_extra
    poa_isotropic = dhi * (1 - anisotropy_index) * 0.5(1 + cosd(surface_tilt) * cosd(surface_roll))
    poa_circumsolar = dhi * (anisotropy_index * projection_ratio)
    return max(0.0, poa_isotropic) + max(0.0, poa_circumsolar)
end

"""
    sky_diffuse_component_reindl(
        solpos::SolarPosition, dni, dni_extra, dhi, ghi,
        surface_tilt, surface_roll, surface_azimuth
    )

Determine the diffuse irradiance from the sky on a tilted surface using the Reindl (1990) 
model [loutzenhiser2007empirical](@cite).
"""
function sky_diffuse_component_reindl(
    solpos::SolarPosition, dni, dni_extra, dhi, ghi,
    surface_tilt, surface_roll, surface_azimuth
)
    aoi_proj = angle_of_incidence_projection(solpos, surface_tilt, surface_roll, surface_azimuth)
    sin_solar_elevation = sind(solpos.elevation)
    anisotropy_index = dni / dni_extra
    rb = max(0.0, aoi_proj) / max(sin_solar_elevation, 0.01745)
    hb = max(0.0, dni * sin_solar_elevation)
    hb_to_ghi = ghi == 0.0 ? 0.0 : hb / ghi

    t1 = 1.0 - anisotropy_index
    t2 = 0.5(1.0 + cosd(surface_tilt) * cosd(surface_roll))
    t3 = 1.0 + sqrt(hb_to_ghi) * (1.0 - cosd(surface_tilt) * cosd(surface_roll))^1.5
    return max(0.0, dhi * (anisotropy_index * rb + t1 * t2 * t3))
end

"""
    sky_diffuse_component_king(
        solpos::SolarPosition, dhi, ghi, surface_tilt, surface_roll, surface_azimuth
    )

Determine diffuse irradiance from the sky on a tilted surface using the King model.
"""
function sky_diffuse_component_king(
    solpos::SolarPosition, dhi, ghi, surface_tilt, surface_roll, surface_azimuth
)
    sky_diffuse = 0.5dhi * (1.0 + cosd(surface_tilt) * cosd(surface_roll)) + 
                  0.5ghi * (0.012 * (90.0 - solpos.elevation - 0.04)) * 
                  (1.0 - cosd(surface_tilt) * cosd(surface_roll))
    return max(0.0, sky_diffuse)
end

function sky_diffuse_component_perez()

end

function sky_diffuse_component_perez_driesse()

end

"""
    ground_diffuse_component(ghi, albedo, surface_tilt, surface_roll)

Estimate diffuse irradiance on a tilted surface from ground reflections 
[loutzenhiser2007empirical](@cite).
"""
function ground_diffuse_component(ghi, albedo, surface_tilt, surface_roll)
    return 0.5 * albedo * ghi * (1.0 - cosd(surface_tilt) * cosd(surface_roll))
end

"""
    total_irradiance( 
        solpos::SolarPosition, irrad::Irradiance,
        surface_tilt, surface_roll, surface_azimuth;
        albedo, sky_diffuse_component
    )

Determine total in-plane irradiance and its beam, sky diffuse and ground
reflected components, using the specified sky diffuse irradiance model.
"""
function total_irradiance(
    solpos::SolarPosition, irrad::Irradiance,
    surface_tilt, surface_roll, surface_azimuth;
    albedo=ALBEDO.urban,
    sky_diffuse_component=sky_diffuse_component_isotropic(irrad.dhi, surface_tilt, surface_roll)
)
    return beam_component(solpos, irrad.dni, surface_tilt, surface_roll, surface_azimuth) + 
           sky_diffuse_component +  
           ground_diffuse_component(irrad.ghi, albedo, surface_tilt, surface_roll)
end

