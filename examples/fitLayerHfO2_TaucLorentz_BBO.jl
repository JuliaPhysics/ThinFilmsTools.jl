# Load modules
using Plots, LaTeXStrings
pyplot()
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
emergent = RIdb.fusedsilicauv(beam.λ)

# Define the RI model to use
layers = [
    LayerTMMO(incident),
    ModelFit(:tauclorentz),
    LayerTMMO(emergent),
]

# Measured absolute transmittance
Texp = SpectraDB.HafniaSpectrum(beam.λ)
plot(beam.λ, Texp)
gui()

## Optimization using BBO

# One oscillator
seed = [vcat(175.0, # thickness
             [3.0, 5.5,
              545, 6.2, 13.2], # osc 1
        ),
]
lb = 0.01.*seed
ub = 2.0.*seed
lb[1][1] = 150
ub[1][1] = 200
lb[1][3] = 4.0
ub[1][3] = 6.0

solOptim = FitTMMOptics(
    Transmittance(Texp), seed, beam, layers;
    alg=:BBO,
    SearchRange=Utils.unfoldbnd(lb,ub),
    Method=:de_rand_1_bin,
)

plot(FitSpectrum(),
    solOptim.Beam.λ, solOptim.spectrumExp, solOptim.spectrumFit,
    xaxis=("Wavelength [nm]"), yaxis=("Transmittance"),
)
gui()
