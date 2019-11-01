# Load modules
using Plots, LaTeXStrings
pyplot()
using ThinFilmsTools
using Optim

# Wavelength range [nm]
λ = 200:1100
# Angle of incidence [degrees]
θ = [0.]
# Polarisation (1.0 = p, 0.0 = s, between 0.0 and 1.0 = average)
pol = 0.5
beam = PlaneWave(λ, θ; p=pol)

# Refractive indices of incident (0) and substrate (2)
incident = RIdb.air(beam.λ)
emergent = RIdb.fusedsilicauv(beam.λ)

# Define the RI model to use
layers = [
    LayerTMMO(incident),
    ModelFit(:forouhibloomer),
    LayerTMMO(emergent),
]

# Measured absolute transmittance
Texp = SpectraDB.ScandiaSpectrum(beam.λ)
plot(beam.λ, Texp)
gui()

options = Optim.Options(
    g_abstol=1e-8, g_reltol=1e-8, iterations=10^3, show_trace=true, store_trace=false,
);

seed = [vcat(90.0, # thickness
        [1.6, 5.0, # ϵinf, Eg
         0.15, 16, 65, # A1, B1, C1
         ],
    ),
]

solOptim = FitTMMOptics(
    Transmittance(Texp), seed, beam, layers;
    options=options, alg=LBFGS(),
)

plot(FitSpectrum(), solOptim.Beam.λ, solOptim.spectrumExp, solOptim.spectrumFit)
gui()

