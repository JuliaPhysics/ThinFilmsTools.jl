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
emergent = RIdb.fusedsilicauv2(beam.λ)

# Define the RI model to use
layers = [
    LayerTMMO(incident),
    ModelFit(:tauclorentz),
    LayerTMMO(emergent),
]

# Measured absolute transmittance
Texp = SpectraDB.HafniaSpectrum(beam.λ)
plot(Spectrum1D(), beam.λ, Texp)
gui()

## Optimization using Optim

options = Optim.Options(
    g_abstol=1e-8, g_reltol=1e-8, iterations=10^5, show_trace=true, store_trace=true,
);

# One oscillator
seed = [vcat(175.0, # thickness
             [3.0, # high-frequency index
               2.5, 2.5, 2.37, 5.0], # osc 1
        ),
]

solOptim = FitTMMOptics(
    Transmittance(), seed, beam, Texp, layers;
    options=options, alg=SAMIN(), lb=0.2.*seed, ub=10.0.*seed,
)

plot(FitSpectrum(),
    solOptim.Beam.λ, solOptim.spectrumExp, solOptim.spectrumFit,
    xaxis=("Wavelength [nm]"), yaxis=("Transmittance"),
)
gui()

# julia> solOptim.optParams
# 1-element Array{Array{Float64,1},1}:
#  [185.57495457185868, 1.4538249665571539, 24.99997160640067, 6.801747171200166, 0.47400001243358125, 4.579667310612569]
#
# julia> solOptim.objfunMin
# 0.0004221112894584093

# Two oscillators
seed2 = [vcat(185.0, # thickness
              [2.0, # high-frequency permittivity
                25.0, 7.0, 0.5, 5.0, # osc 1
                10.0, 7.0, 0.5, 5.0, # osc 2
              ]
        ),
]

solOptim2 = FitTMMOptics(
    Transmittance(), seed2, beam, Texp, layers;
    options=options, alg=SAMIN(), lb=0.1.*seed2, ub=10.0.*seed2,
)

plot(FitSpectrum(),
    solOptim2.Beam.λ, solOptim2.spectrumExp, solOptim2.spectrumFit,
    xaxis=("Wavelength [nm]"), yaxis=("Transmittance"),
)
gui()

# julia> solOptim2.optParams
# 1-element Array{Array{Float64,1},1}:
#  [189.09445508077323, 0.39609259415335923, 158.43421436032554, 8.610532933615291, 3.4744910306832337, 5.015700019476466, 15.249954668687487, 0.7000170470355342, 0.06250726388751021, 1.0341375083174082]
#
# julia> solOptim2.objfunMin
# 0.0001506570743847816


