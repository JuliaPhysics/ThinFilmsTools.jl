using Optim
using Plots
pyplot(reuse=false, size=(640, 480), grid=false)
closeall()
using ThinFilmsTools

##
function glassTransmittance(beam, incident, emergent)
    # Create a measured spectrum of reflection
    layers2 = [ LayerTMMO1DIso(incident),
                LayerTMMO1DIso(RIdb.glass(beam.λ./1e3); d=400.),
                LayerTMMO1DIso(emergent) ]
    sol = TMMO1DIsotropic(beam, layers2)
    return beam.p.*sol.Spectra.Tp .+ (1.0 - beam.p).*sol.Spectra.Ts
end
##

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
           ModelFit(:cauchyurbach),
           LayerTMMO1DIso(emergent) ]

# Create transmittance spectrum to fit
Texp = glassTransmittance(beam, incident, emergent)
plot(TMMOPlotSpectra1D(), beam.λ, Texp)
gui()

options = Optim.Options(g_abstol=1e-8, g_reltol=1e-8, iterations=10^5, show_trace=true, store_trace=true);
seed = [vcat(420.0, [1.4, 0.23, 1.0, 0.0, 1.0, 1.0])]
solOptim = FitTMMO1DIsotropic(Transmittance(), seed, beam, Texp, layers; options=options, alg=SAMIN(), lb=0.5.*seed, ub=1.5.*seed)

plot(PlotFitSpectrum(), solOptim.beam.λ, solOptim.spectrumExp, solOptim.spectrumFit, xaxis=("Wavelength [nm]"), yaxis=("Transmittance"), tickfont=font(12), legendfont=font(10))
gui()

