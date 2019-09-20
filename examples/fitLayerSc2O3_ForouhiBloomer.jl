using Optim
using Plots
pyplot(reuse=false, size=(640, 480), grid=false)
closeall()
using ThinFilmsTools

# Wavelength range [nm]
λ = 200:1100
# Angle of incidence [degrees]
θ = [0.]
# Polarisation (1.0 = p, 0.0 = s, between 0.0 and 1.0 = average)
pol = 0.5
beam = PlaneWave(λ, θ; p=pol)

# Refractive indices of incident (0) and substrate (2)
incident = RIdb.air(beam.λ)
emergent = RIdb.fusedsilicauv2(beam.λ)

# Define the RI model to use
layers = [ LayerTMMO1DIso(incident),
           ModelFit(:forouhibloomer),
           LayerTMMO1DIso(emergent) ]

# Measured absolute transmittance
Texp = Utils.getScandiaSpectrum(beam.λ)./100
plot(TMMOPlotSpectra1D(), beam.λ, Texp)
gui()

## Optimization using Optim

options = Optim.Options(g_abstol=1e-8, g_reltol=1e-8, iterations=10^5, show_trace=true, store_trace=true);
seed = [vcat(95.0, [0.13, 15.0, 63.0, 4.0, 2.0])]
solOptim = FitTMMO1DIsotropic(Transmittance(), seed, beam, Texp, layers; options=options, alg=SAMIN(), lb=0.1.*seed, ub=2.0.*seed)

plot(PlotFitSpectrum(), solOptim.beam.λ, solOptim.spectrumExp, solOptim.spectrumFit, xaxis=("Wavelength [nm]"), yaxis=("Transmittance"))
gui()

