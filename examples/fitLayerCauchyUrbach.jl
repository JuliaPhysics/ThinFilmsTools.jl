# Load modules
using Plots, LaTeXStrings
pyplot()
using ThinFilmsTools
using Optim

##
function glass_transmittance(beam, incident, emergent)
    # Create a measured spectrum of reflection
    layers2 = [ LayerTMMO(incident),
                LayerTMMO(RIdb.glass(beam.λ./1e3); d=400.),
                LayerTMMO(emergent) ]
    sol = tmm_optics(beam, layers2)
    return vec(beam.p.*sol.Spectra.Tp .+ (1.0 - beam.p).*sol.Spectra.Ts)
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
layers = [
    LayerTMMO(incident), # 1
    ModelFit(:cauchyurbach), # 2
    LayerTMMO(emergent), # 3
]

# Create transmittance spectrum to fit
Texp = glass_transmittance(beam, incident, emergent)
# plot(beam.λ, Texp)
# gui()

options = Optim.Options(
    g_abstol=1e-8, g_reltol=1e-8, iterations=10^5, show_trace=true, store_trace=true,
);

seed = [vcat(420.0, [1.4, 0.23, 1.0, 0.0, 1.0, 1.0])]

solOptim = fit_tmm_optics(
    Transmittance(Texp), seed, beam, layers, order;
    options=options, alg=SAMIN(), lb=0.5.*seed, ub=1.5.*seed,
)

plot(FitSpectrum(),
    solOptim.beam.λ, solOptim.spectrumExp, solOptim.spectrumFit,
    xaxis=("Wavelength [nm]"),
    yaxis=("Transmittance"),
)
gui()
