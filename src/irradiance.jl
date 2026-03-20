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
    angle_of_incidence_projection(p::Panel, solpos::SolarPosition)

Calculates the dot product of the sun position unit vector and the surface normal unit 
vector, i.e., the cosine of the angle of incidence.
"""
function angle_of_incidence_projection(p::Panel, solpos::SolarPosition)
    return cosd(p.roll) * cosd(p.tilt) * sind(solpos.apparent_elevation) +
           cosd(solpos.apparent_elevation) * sind(solpos.azimuth) * (
                sind(p.tilt) * sind(p.azimuth) + 
                cosd(p.azimuth) * cosd(p.tilt) * sind(p.roll)
           ) +
           cosd(solpos.azimuth) * cos(solpos.apparent_elevation) * (
                cosd(p.azimuth) * sin(p.tilt) - 
                sind(p.azimuth) * cosd(p.tilt) * sind(p.roll)
           )
end

"""
    beam_component(p::Panel, solpos::SolarPosition, dni)

Calculates the beam component of the plane of array irradiance.
"""
function beam_component(p::Panel, solpos::SolarPosition, dni)
    aoi_proj = angle_of_incidence_projection(p, solpos) 
    return aoi_proj * dni
end

"""
    sky_diffuse_component_isotropic(p::Panel, dhi)

Determine diffuse irradiance from the sky on a tilted surface using the isotropic sky 
model [loutzenhiser2007empirical, kamphuis2020perspectives](@cite).
"""
function sky_diffuse_component_isotropic(p::Panel, dhi)
    return 0.5 * dhi * (1 + cosd(p.tilt) * cosd(p.roll)) 
end

"""
    sky_diffuse_component_klucher(p::Panel, solpos::SolarPosition, dhi, ghi)

Determine diffuse irradiance from the sky on a tilted surface using the Klucher (1979) 
model [klucher1979evaluation, loutzenhiser2007empirical](@cite).
"""
function sky_diffuse_component_klucher(p::Panel, solpos::SolarPosition, dhi, ghi)
    F = 1 - (dhi / ghi)^2
    aoi_proj = max(0.0, angle_of_incidence_projection(p, solpos))
    return 0.5(1 + cosd(p.tilt) * cosd(p.roll)) *
           (1 + F * (1 - cosd(p.tilt) * cosd(p.roll))^1.5) *
           (1 + F * aoi_proj^2 * sind(solpos.zenith)^3)
end

"""
    sky_diffuse_component_haydavies(p::Panel, solpos::SolarPosition, dhi, dhi_extra, dhi)

Determine diffuse irradiance from the sky on a tilted surface using the Hay and Davies 
(1980) model [hay5972calculations, loutzenhiser2007empirical](@cite).
"""
function sky_diffuse_component_haydavies(p::Panel, solpos::SolarPosition, dni, dni_extra, dhi)
    aoi_proj = max(0, angle_of_incidence_projection(p, solpos))
    projection_ratio = aoi_proj / max(cos(solpos.zenith), 0.01745)
    anisotropy_index = dni / dni_extra
    poa_isotropic = max(0.0, dhi * (1 - anisotropy_index) * 0.5(1 + cosd(p.tilt) * cosd(p.roll)))
    poa_circumsolar = max(0.0, dhi * (anisotropy_index * projection_ratio))
    return poa_isotropic + poa_circumsolar
end

"""
    sky_diffuse_component_reindl(p::Panel, solpos::SolarPosition, dni, dhi_extra, dhi, ghi)

Determine the diffuse irradiance from the sky on a tilted surface using the Reindl (1990) 
model [loutzenhiser2007empirical](@cite).
"""
function sky_diffuse_component_reindl(p::Panel, solpos::SolarPosition, dni, dhi_extra, dhi, ghi)
    aoi_proj = max(0.0, angle_of_incidence_projection(p, solpos))
    cos_solar_zenith = cosd(solpos.zenith)
    anisotropy_index = dni / dni_extra
    rb = cos_solar_zenith / max(cos_solar_zenith, 0.01745)
    hb = np.maximum(0.0, dni * cos_solar_zenith)
    hb_to_ghi = ghi == 0.0 ? 0.0 : hb / ghi

    t1 = 1 - anisotropy_index
    t2 = 0.5(1 + cosd(p.tilt) * cos(p.roll))
    t3 = 1 + sqrt(hb_to_ghi) * (1 - cosd(p.tilt) * cosd(p.roll))^1.5
    return max(0.0, dhi * (anisotropy_index * rb + t1 * t2 * t3))
end

"""
    sky_diffuse_component_king(p::Panel, solpos::SolarPosition, dhi, ghi)

Determine diffuse irradiance from the sky on a tilted surface using the King 
model.
"""
function sky_diffuse_component_king(p::Panel, solpos::SolarPosition, dhi, ghi)
    sky_diffuse = 0.5dhi * (1 + cosd(p.tilt) * cosd(p.roll)) + 
                  0.5ghi * (0.012 * solpos.zenith - 0.04) * (1 - cosd(p.tilt) * cos(p.roll))
    return max(0.0, sky_diffuse)
end

function sky_diffuse_component_perez()

end

function sky_diffuse_component_perez_driesse()

end

"""
    ground_diffuse_component(p::Panel, albedo, ghi)

Estimate diffuse irradiance on a tilted surface from ground reflections 
[loutzenhiser2007empirical](@cite).
"""
function ground_diffuse_component(p::Panel, albedo, ghi)
    return 0.5 * albedo * ghi * (1 - cosd(p.tilt) * cosd(p.roll))
end

"""
total_irradiance(
    p::Panel, 
    solpos::SolarPosition, 
    irrad::Irradiance;
    albedo,
    sky_diffuse_component
)
"""
function total_irradiance(
    p::Panel, 
    solpos::SolarPosition, 
    irrad::Irradiance;
    albedo=ALBEDO.urban,
    sky_diffuse_component=sky_diffuse_component_isotropic(p, irrad.dhi)
)
    return beam_component(p, solpos, irrad.dni) + 
           sky_diffuse_component +  
           ground_diffuse_component(p, albedo, irrad.ghi)
end

