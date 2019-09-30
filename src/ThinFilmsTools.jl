module ThinFilmsTools

include("Utils.jl")
export Utils

include("RefractiveIndicesDB.jl")
export RIdb

include("RefractiveIndicesModels.jl")
export RI

include("CommonStructures.jl")
using .CommonStructures
export PlaneWave,
       LayerTMMO1DIso,
       ModelFit,
       BoundariesFit

include("SpectraDBExamples.jl")
export SpectraDB

include("MatrixMethod1DIsotropic.jl")
using .MatrixMethod1DIsotropic
export TMMO1DIsotropic

include("FitThinFilmSpectrum1DIsotropic.jl")
using .FitThinFilmPattern
export SpaceSolutionEMA,
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

include("ThreeOmegaMethod.jl")
using .ThreeΩMethod
export ThreeOmegaMethod,
       LayerTOM,
       HeaterGeometry,
       Source

include("PlottingTools.jl")
using .PlottingTools
export TOMPlot,
       TMMOPlotSpectra1D,
       TMMOPlotSpectraAngle1D,
       TMMOPlotSpectra2D,
       TMMOPlotEMF2D,
       TMMOPlotEMFAngle2D,
       TMMOPlotNprofile,
       TMMOPlotDispersion,
       PlotFitSpectrum,
       SpaceSolutionOFplot

include("FitBellShapeCurves.jl")
using .FitBellShapeCurves
export CurveModelSolution,
       FitCurveModel

end # ThinFilmsTools
