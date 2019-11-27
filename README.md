# ThinFilmsTools.jl

[![The MIT License](https://img.shields.io/badge/license-MIT-orange.svg?style=flat-square)](http://opensource.org/licenses/MIT)
[![Build Status](https://travis-ci.com/lnacquaroli/ThinFilmsTools.jl.svg?branch=master)](https://travis-ci.com/lnacquaroli/ThinFilmsTools.jl)
[![DOI](https://zenodo.org/badge/200238502.svg)](https://zenodo.org/badge/latestdoi/200238502)


[ThinFilmsTools.jl](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/Home) provides tools for the design and characterisation of thin films written in Julia. Documentation, examples and details can be found in the [Wiki](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/Home) page.

Check [julia language](https://julialang.org/) for how to install and setup the environment.

## Installation

This package is not yet registered. It can be installed in Julia with the following ([see further](https://docs.julialang.org/en/v1/stdlib/Pkg/index.html#Adding-unregistered-packages-1)):
```julia
julia> ]
pkg> add https://github.com/lnacquaroli/ThinFilmsTools.jl
```

[ThinFilmsTools.jl](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/Home) is compatible with Julia version 1.1 or later.

The package includes modules for the calculation of parameters of thin films: [tmm_optics](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/TMMOptics) for the simulation of optical properties of thin films, [fit_tmm_optics](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/FitTMMOptics) to fit thin films spectrum and [three_omega](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/ThreeOmega) to model the thermal properties of thin films based on the 3ω method.

The package also contains a number of indices of refraction for different materials in a database [RefractiveIndicesDB.jl](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/RefractiveIndicesDB.jl) ready to use.

For the simulation of index of refraction mixtures, there is also available a database [RefractiveIndicesModels.jl](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/RefractiveIndicesModels.jl) that contains several mixing rules for this.

The [fit_curve_model](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/FitCurveModel) allows fitting the spectral line shapes of different processes using a simple function. It supports three models for an arbitrary number of peaks: Gaussian, Lorentzian, and Voigtian.

A bunch of functions and recipes are included ([PlottingTools.jl](https://github.com/lnacquaroli/ThinFilmsTools.jl/wiki/PlottingTools.jl)) for convenience to plot results from fit_tmm_optics, tmm_optics and three_omega.
