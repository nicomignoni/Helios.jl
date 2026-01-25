# API

## Solar position algorithm
The following functions implement the Solar Position Algorithm, as described in 
[reda2004solar](@cite). 

```@autodocs
Modules = [SolarFunctions.Model]
Pages = ["src/spa.jl"]
```

## Atmosphere

```@autodocs
Modules = [SolarFunctions.Model]
Pages = ["src/atmosphere.jl"]
```

## Clear-sky
The following functions implement different models to calculate the irradiance components 
under the clear-sky assumption.

```@autodocs
Modules = [SolarFunctions.Model]
Pages = ["src/clearsky.jl"]
```

## Irradiance

```@autodocs
Modules = [SolarFunctions.Model]
Pages = ["src/irradiance.jl"]
```
