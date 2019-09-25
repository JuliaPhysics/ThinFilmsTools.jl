# ThinFilmsTools.jl

[![The MIT License](https://img.shields.io/badge/license-MIT-orange.svg?style=flat-square)](http://opensource.org/licenses/MIT)
[![Build Status](https://travis-ci.com/lnacquaroli/ThinFilmsTools.jl.svg?branch=master)](https://travis-ci.com/lnacquaroli/ThinFilmsTools.jl)
[![DOI](https://zenodo.org/badge/200238502.svg)](https://zenodo.org/badge/latestdoi/200238502)


[ThinFilmsTools.jl](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/Home) provides tools for the design and characterisation of thin films written in Julia. More documentation, examples and details can be found in the [Wiki](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/Home) page.

## Installation

This package is not yet registered. It can be installed in Julia with the following ([see further](https://docs.julialang.org/en/v1/stdlib/Pkg/index.html#Adding-unregistered-packages-1)):
```julia
julia> ]
pkg> add https://github.com/lnacquaroli/ThinFilmsTools.jl
```

[ThinFilmsTools.jl](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/Home) is compatible with Julia version 1.1 or later.

The package includes modules for the calculation of parameters of thin films: [TMMO1DIsotropic](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/TMMO1DIsotropic) for the simulation of optical properties of thin films, [FitTMMO1DIsotropic](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/FitTMMO1DIsotropic) to fit thin films spectrum and [ThreeOmegaMethod](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/ThreeOmegaMethod) to model the thermal properties of thin films based on the 3ω method.

The package also contains a number of indices of refraction for different materials in a database [RefractiveIndicesDB.jl](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/RefractiveIndicesDB.jl) ready to use.

For the simulation of index of refraction mixtures, there is also available a database [RefractiveIndicesModels.jl](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/RefractiveIndicesModels.jl) that contains several mixing rules for this.

A bunch of functions and recipes are included ([PlottingTools.jl](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/PlottingTools.jl)) for convenience to plot results from FitTMMO1DIsotropic, TMMO1DIsotropic and ThreeOmegaMethod.
