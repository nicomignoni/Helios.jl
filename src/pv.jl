"""
    surface_normal(surface_tilt, surface_roll, surface_azimuth)

Computes the normal (unit) vector of the oriented surface.
"""
function surface_normal(surface_tilt, surface_roll, surface_azimuth)
    [
        sind(surface_tilt) * sind(surface_azimuth) + 
            cosd(surface_tilt) * cosd(surface_azimuth) * sind(surface_roll),
        sind(surface_tilt) * cosd(surface_azimuth) - 
            cosd(surface_tilt) * sind(surface_azimuth) * sind(surface_roll),
        cosd(surface_tilt) * cosd(surface_roll)
    ]
end
