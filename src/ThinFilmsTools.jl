module ThinFilmsTools

using Statistics
using LinearAlgebra
using Printf
using QuadGK
using Plots # using RecipesBase
using LaTeXStrings
using Interpolations
using Optim

include("CommonStructures.jl")
export PlaneWave,
       LayerTMMO1DIso,
       ModelFit,
       BoundariesFit

include("CommonUtils.jl")
export Utils

include("TMMO1DIsotropic.jl")
export TMMO1DIsotropic

include("ThreeOmegaMethod.jl")
export ThreeOmegaMethod,
       LayerTOM,
       HeaterGeometry,
       Source

include("RefractiveIndicesDB.jl")
export RIdb

include("RefractiveIndicesModels.jl")
export RI

include("SpectraDBExamples.jl")
export SpectraDB

include("FitThinFilmSpectrum.jl")
export SpaceSolution2D,
       TheoreticalSpectrum,
       NormalizeReflectance,
       FitTMMO1DIsotropic,
       Reflectance,
       Transmittance,
       Generic,
       DBR,
       DBRAlpha,
       MC1d,
       MC1dAlpha

include("PlottingTools.jl")
export TOMplot,
       TMMOPlotSpectra1D,
       TMMOPlotSpectraAngle1D,
       TMMOPlotSpectra2D,
       TMMOPlotEMF2D,
       TMMOPlotEMFAngle2D,
       TMMOPlotNprofile,
       TMMOPlotDispersion,
       PlotFitSpectrum,
       SpaceSolutionOFplot,
       PlotFitSpectrumEllip

end # ThinFilmsTools
