using Optim
using Plots
pyplot(reuse=false, size=(640, 480), grid=false)
closeall()
using ThinFilmsTools

# Wavelength range [nm]
λ = 400:1000
# Angle of incidence [degrees]
θ = [5.]
# Polarisation (1.0 = p, 0.0 = s, between 0.0 and 1.0 = average)
pol = 0.5
beam = PlaneWave(λ, θ; p=pol)

# Refractive indices of incident (0) and substrate (2)
incident = RIdb.air(beam.λ)
emergent = RIdb.silicon(beam.λ)

# Define the RI model to use
layers = [ LayerTMMO1DIso(incident),
           ModelFit(:looyengaspheresbin; N=(incident, emergent)),
           ModelFit(:looyengaspheresbin; N=(incident, emergent)),
           LayerTMMO1DIso(emergent) ]

# Absolute reflectance spectrum to fit stored in Utils
Rexp = Utils.getBraggSpectrum(beam.λ)
plot(TMMOPlotSpectra1D(), beam.λ, Rexp)
gui()

options = Optim.Options(g_abstol=1e-8, g_reltol=1e-8, iterations=10^5, show_trace=true, store_trace=true);

seed = [[69, 0.54], [119.0, 0.77]]
solOptim = FitTMMO1DIsotropic(Reflectance(), seed, beam, Rexp, layers; arrange=DBR(), L=[9], options=options, alg=SAMIN(), lb=0.5.*seed, ub=1.5.*seed)
plot(PlotFitSpectrum(), solOptim.beam.λ, solOptim.spectrumExp, solOptim.spectrumFit)
gui()

# Optimisation using alpha parameter
seed2 = [[69, 0.54], [119.0, 0.77], [0.995]]
solOptim2 = FitTMMO1DIsotropic(Reflectance(), seed2, beam, Rexp, layers; arrange=DBRAlpha(), L=[9], options=options, alg=SAMIN(), lb=0.5.*seed2, ub=1.5.*seed2)
plot(PlotFitSpectrum(), solOptim2.beam.λ, solOptim2.spectrumExp, solOptim2.spectrumFit)
gui()

