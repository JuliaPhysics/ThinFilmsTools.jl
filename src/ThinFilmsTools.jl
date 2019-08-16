module ThinFilmsTools

using Statistics
using Printf: @sprintf
using QuadGK
using Plots
using LaTeXStrings # To use with PlottingTools.jl
pyplot(grid=false)

include("TMMO1DIsotropic.jl")
export TMMO1DIsotropic, PlaneWave, LayerTMMO1DIso, TransferMatrix

include("ThreeOmegaMethod.jl")
export LayerTOM, HeaterGeometry, Source, ThreeOmegaMethod

include("PlottingTools.jl")
export TOMplot, TMMOplotSpectra1D, TMMOplotSpectraAngle1D
export TMMOplotSpectra2D, TMMOplotEMF2D, TMMOplotEMFAngle2D
export TMMOplotNprofile, TMMOplotdispersion

include("RefractiveIndicesDB.jl")
export RIdb

include("MixingRules.jl")
export DF

end # ThinFilmsTools
