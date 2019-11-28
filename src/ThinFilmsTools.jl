module ThinFilmsTools

include("Utils.jl")
export Utils

include("CommonStructures.jl")
using .CommonStructures
export PlaneWave,
       LayerTMMO,
       ModelFit,
       BoundariesFit

include("RefractiveIndicesDB.jl")
export RIdb

include("RefractiveIndicesModels.jl")
export RI

include("SpectraDBExamples.jl")
export SpectraDB

include("TransferMatrixMethod.jl")
using .TransferMatrixMethod
export tmm_optics, TMMOptics

include("FTFSpectrum.jl")
using .FTFSpectrum
export FitTMMOptics,
       Reflectance,
       Transmittance,
       Ellipsometry,
       NoAlpha,
       UseAlpha,
       MeanAbs,
       SumAbs,
       SumMeanAbs,
       space_solution_ema,
       theoretical_spectrum,
       normalize_reflectance,
       fit_tmm_optics

include("ThreeOmegaMethod.jl")
using .ThreeOmegaMethod
export ThreeOmega,
       LayerTOM,
       HeaterGeometry,
       Source,
       three_omega

include("PlottingTools.jl")
using .PlottingTools
export TOMPlot,
       Spectrum1D,
       SpectrumAngle1D,
       Spectrum2D,
       EMF2D,
       EMFAngle2D,
       FitSpectrum,
       FitSpectrumEllip,
       SpaceSolution,
       RIprofile,
       PBGDispersion1D,
       PBGDispersion1Dalt,
       PBGDispersion1Dimre,
       PBGDispersion2D,
       PBGDispersion2Dalt

include("FitBellCurves.jl")
using .FitBellCurves
export FitCurveModel, fit_curve_model

end # module
