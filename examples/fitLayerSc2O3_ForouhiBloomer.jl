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
    ModelFit(:forouhibloomer),
    LayerTMMO(emergent),
]

# Measured absolute transmittance
Texp = SpectraDB.ScandiaSpectrum(beam.λ)
plot(Spectrum1D(), beam.λ, Texp)
gui()

## Optimization using Optim

options = Optim.Options(
    g_abstol=1e-8, g_reltol=1e-8, iterations=10^5, show_trace=true, store_trace=true,
);

# One oscillator
seed = [vcat(95.0, # thickness
              [2.0, # high-frequency index
               0.2, 11.0, 69.0, 3.0], # osc 1
        ),
]
solOptim = FitTMMOptics(
    Transmittance(), seed, beam, Texp, layers;
    options=options, alg=SAMIN(), lb=0.2.*seed, ub=2.0.*seed,
)
plot(FitSpectrum(),
    solOptim.Beam.λ, solOptim.spectrumExp, solOptim.spectrumFit,
    xaxis=("Wavelength [nm]"), yaxis=("Transmittance"),
)
gui()

# julia> solOptim.optParams
# 1-element Array{Array{Float64,1},1}:
 # [91.1219205143226, 1.5418143348321705, 0.2043718097297613, 17.43135747991159, 77.71747716440301, 4.35561678410554]
#
# julia> solOptim.objfunMin
# 6.394348409151524e-5

# Two oscillators
seed2 = [vcat(95.0, # thickness
              [2.0, # high-frequency index
               0.2, 11.0, 69.0, 3.0, # osc 1
               0.2, 11.0, 69.0, 3.0, # osc 2
              ]
        ),
]
solOptim2 = FitTMMOptics(
    Transmittance(), seed2, beam, Texp, layers;
    options=options, alg=SAMIN(), lb=0.1.*seed2, ub=2.0.*seed2,
)
plot(FitSpectrum(),
    solOptim2.Beam.λ, solOptim2.spectrumExp, solOptim2.spectrumFit,
    xaxis=("Wavelength [nm]"), yaxis=("Transmittance"),
)
gui()

# julia> solOptim2.objfunMin
# 6.153001209410525e-5
#
# julia> solOptim2.optParams
# 1-element Array{Array{Float64,1},1}:
#  [91.79294036278226, 1.502110799503048, 0.04461148596620742, 16.929726884406634, 72.18163776104026, 3.4149721274723626, 0.22343536771187994, 19.796636189142653, 108.40252015354704, 4.806383876703072]
