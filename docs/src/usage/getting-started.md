# Getting started

```@example QUICKSTART
using Dates, Plots, Helios

# Create a location and a year-long time range
location = Location(41.11148, 16.8554; altitude=6) # latitude, longitude, altitude [meters]
time_range = let n = now(); n:Day(1):(n + Year(1)); end

# Compute azimuth and apparent elevation for such location and time range
solpos = Vector{SolarPosition}(undef, length(time_range))
for (i, dt) in enumerate(time_range)
    solpos[i] = spa(location, dt)
end

# Let's create an analemma
analemma = scatter(
    getfield.(solpos, :azimuth),
    getfield.(solpos, :apparent_elevation),
    xlabel = "Solar azimuth [deg]",
    ylabel = "Solar elevation [deg]",
    legend = false,
    grid = true,
)
```

```@example QUICKSTART
time_range = let today = DateTime(today()); today:Minute(1):(today + Day(1)); end

# Compute the irradiance for daylight hours
daylight_times = DateTime[] 
irradiance = (
    ineichen = Irradiance[],
    simplified_solis = Irradiance[]
)
for dt in time_range
    solpos = spa(location, dt)
    if solpos.apparent_elevation > 0
        push!(daylight_times, dt)
        # using the Ineichen model
        push!(irradiance.ineichen, clearsky_ineichen(location, dt; solpos=solpos))
        # and the simplified_solis model
        push!(irradiance.simplified_solis, clearsky_simplified_solis(location, dt, solpos=solpos))
    end
end

# and plot its components
plt = plot(xlabel = "Time", ylabel = "Irradiance [W/m^2]", grid=true, layout=(2,1), linx=:x)
for (c, l) in ((:dni, "DNI"), (:dhi, "DHI"), (:ghi, "GHI"))
    plot!(plt, Time.(daylight_times), getfield.(irradiance.ineichen, c), label=l, lw=3, subplot=1)
    plot!(plt, Time.(daylight_times), getfield.(irradiance.simplified_solis, c), label=l, lw=3, subplot=2)
end

plt
```
