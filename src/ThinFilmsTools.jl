module ThinFilmsTools

using Printf: @sprintf
using QuadGK
using Plots
using LaTeXStrings
pyplot(grid=false)

include("TMMOptics.jl")
export TMMOptics, PlaneWave, LayerTMMO

include("ThreeOmegaMethod.jl")
export LayerTOM, HeaterGeometry, Source, ThreeOmegaMethod

include("PlottingTools.jl")
export TOMplot, TMMOplotSpectra1D, TMMOplotSpectraAngle1D, TMMOplotSpectra2D, TMMOplotEMF2D, TMMOplotEMFAngle2D, TMMOplotNprofile, TMMOplotdispersion

end # ThinFilmsTools
