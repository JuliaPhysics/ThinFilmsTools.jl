module ThinFilmsTools

using Printf: @sprintf
using QuadGK

include("TMMO1DIsotropic.jl")
export TMMO1DIsotropic, PlaneWave, LayerTMMO1DIso

include("ThreeOmegaMethod.jl")
export LayerTOM, HeaterGeometry, Source, ThreeOmegaMethod

using Plots, LaTeXStrings # To use with PlottingTools.jl
pyplot(grid=false)
#gr(size=(600,450), grid=false, format=:svg)
include("PlottingTools.jl")
export TOMplot, TMMOplotSpectra1D, TMMOplotSpectraAngle1D, TMMOplotSpectra2D, TMMOplotEMF2D, TMMOplotEMFAngle2D, TMMOplotNprofile, TMMOplotdispersion

# datapath = joinpath(@__DIR__, "data/")
# include(datapath * "RIdb.jl")
include("RefractiveIndicesDB.jl")
export RIdb

include("MixingRules.jl")
export DF

end # ThinFilmsTools
