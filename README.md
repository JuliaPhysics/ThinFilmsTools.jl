# ThinFilmsTools.jl

[![The MIT License](https://img.shields.io/badge/license-MIT-orange.svg?style=flat-square)](http://opensource.org/licenses/MIT)
[![Build Status](https://travis-ci.com/lnacquaroli/ThinFilmsTools.jl.svg?branch=master)](https://travis-ci.com/lnacquaroli/ThinFilmsTools.jl)

[ThinFilmsTools.jl](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/Home) provides tools for the design and characterisation of thin films written in Julia. More documentation, examples and details can be found in the [Wiki](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/Home) page.

The are two major methods working for this already: [TMMO1DIsotropic.jl](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/TMMO1DIsotropic.jl) (old [TMMOptics.jl](https://github.com/lnacquaroli/TMMOptics.jl)) and [ThreeOmegaMethod.jl](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/ThreeOmegaMethod.jl). Looking forward to expand it with other techniques.

## Installation

This package is not yet registered. It can be installed in Julia with the following ([see further](https://docs.julialang.org/en/v1/stdlib/Pkg/index.html#Adding-unregistered-packages-1)):
```julia
julia> ]
(v1.1) pkg> add https://github.com/lnacquaroli/ThinFilmsTools.jl
```

[ThinFilmsTools.jl](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/Home) is compatible with Julia version 1.1 or later, and includes two main modules for the calculation of parameters of thin films: [TMMO1DIsotropic.jl](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/TMMO1DIsotropic.jl) (replacing [TMMOptics.jl](https://github.com/lnacquaroli/TMMOptics.jl)) for the simulation of optical properties of thin films, [FitReflectance.jl](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/FitReflectance.jl) to fit thin films reflectance spectrum, [FitTransmittance.jl](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/FitTransmittance.jl) to fit thin films transmittance spectrum and [ThreeOmegaMethod.jl](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/ThreeOmegaMethod.jl) to model the thermal properties of thin films based on the 3ω method.

The package also contains a number of indices of refraction for different materials in a database [RefractiveIndicesDB.jl](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/RefractiveIndicesDB.jl) ready to use.

For the simulation of index of refraction mixtures, there is also available a database [MixingRules.jl](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/MixingRules.jl) that contains several mixing rules for this.

A bunch of functions and recipes are included ([PlottingTools.jl](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/PlottingTools.jl)) for convenience to plot results from TMMO1DIsotropic.jl, FitReflectance.jl FitTransmittance.jl and ThreeOmegaMethod.jl.

