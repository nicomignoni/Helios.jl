using StaticArrays

mutable struct Panel{T}
    size::SVector{2, T}
    pos::SVector{3, T}
    tilt::T
    roll::T
    azimuth::T
end

function Panel(size; pos=zeros(SVector{3}), tilt=0.0, roll=0.0, azimuth=180.0)
    Panel(size, pos, tilt, roll, azimuth)
end

function area(p::Panel)
    prod(p.size) 
end

function normal(p::Panel)
    [
        sind(p.tilt) * sind(p.azimuth) + cosd(p.tilt) * cosd(p.azimuth) * sind(p.roll),
        sind(p.tilt) * cosd(p.azimuth) - cosd(p.tilt) * sind(p.azimuth) * sind(p.roll),
        cosd(p.tilt) * cosd(p.roll)
    ]
end
