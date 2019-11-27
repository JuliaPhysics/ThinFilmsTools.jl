# Load modules
using Plots, LaTeXStrings
pyplot()
using ThinFilmsTools
using Optim

# Wavelength range [nm]
λ = 200:2100
# Angle of incidence [degrees]
θ = [60.]
# Polarisation (1.0 = p, 0.0 = s, between 0.0 and 1.0 = average)
pol = 0.5
beam = PlaneWave(λ, θ; p=pol)

# Refractive indices of incident (0) and substrate (2)
incident = RIdb.air(beam.λ)
emergent = RIdb.fused_silica_uv(beam.λ)

# Define the RI model to use
layers = [
    LayerTMMO(incident),
    ModelFit(:tauclorentz),
    LayerTMMO(emergent),
]

# Measured absolute transmittance
Ψexp, Δexp, = SpectraDB.hafnia_ellips(beam.λ)
plot(beam.λ, Ψexp)
plot!(beam.λ, Δexp)
gui()

# Spectrum type: ellipsometry with the input experimental spectra
spectype = Ellipsometry([Ψexp Δexp])

## Optimization using BBO

# One oscillator
seed = [vcat(100.0, # thickness
             [3.0, 5.5,
              545, 6.2, 13.2], # osc 1
        ),
]
lb = 0.01.*seed
ub = 2.0.*seed
lb[1][1] = 80
ub[1][1] = 200
lb[1][3] = 4.0
ub[1][3] = 6.0

solOptim = fit_tmm_optics(
    spectype, seed, beam, layers;
    alg=:BBO,
    SearchRange=Utils.unfoldbnd(lb,ub),
)

plot(FitSpectrumEllip(),
    solOptim.beam.λ, solOptim.spectrumExp, solOptim.spectrumFit,
    xaxis=("Wavelength [nm]"), yaxis=(),
)
gui()
