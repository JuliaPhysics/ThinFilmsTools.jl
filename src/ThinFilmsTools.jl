module ThinFilmsTools

using Printf: @sprintf
using QuadGK

include("TMMO1DIsotropic.jl")
export TMMO1DIsotropic, PlaneWave, LayerTMMO1DIso

include("ThreeOmegaMethod.jl")
export LayerTOM, HeaterGeometry, Source, ThreeOmegaMethod

using Plots, LaTeXStrings # To use with PlottingTools.jl
pyplot(grid=false)
include("PlottingTools.jl")
export TOMplot, TMMOplotSpectra1D, TMMOplotSpectraAngle1D, TMMOplotSpectra2D, TMMOplotEMF2D, TMMOplotEMFAngle2D, TMMOplotNprofile, TMMOplotdispersion

include("RefractiveIndicesDB.jl")
export RIdb

include("MixingRules.jl")
export DF

end # ThinFilmsTools
