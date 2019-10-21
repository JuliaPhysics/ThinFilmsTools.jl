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
    ModelFit(:forouhibloomermodified),
    LayerTMMO(emergent),
]

# Measured absolute transmittance
Texp = SpectraDB.HafniaSpectrum(beam.λ)
# plot(Spectrum1D(), beam.λ, Texp)
# gui()

## Optimization using Optim

options = Optim.Options(
    g_abstol=1e-8, g_reltol=1e-8, iterations=10^5, show_trace=true, store_trace=true,
);

# One oscillator
seed = [vcat(175.0, # thickness
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
#  [187.183792816788, 4.252808817196439, 4.999999998836596, 4.99999998125792, 3.123370291788414, 3.626953034098114]
#
# julia> solOptim.objfunMin
# 0.003391222097106318

# Two oscillators
seed2 = [vcat(185.0, # thickness
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

# julia> solOptim2.optParams
# 1-element Array{Array{Float64,1},1}:
#  [196.9483935634827, 3.4299090967378145, 9.999954258245525, 6.766755402479663, 0.30000000518861997, 2.4550653312529804, 1.5700141364784361, 0.4619717620787732, 0.39815650814881204, 0.5000174726770397]
#
# julia> solOptim2.objfunMin
# 0.00042755324878245287

# Three oscillators
seed3 = [vcat(195.0, # thickness
              [4.0, # high-frequency permittivity
               10.0, 7.0, 0.3, 3.0, # osc 1
               1.5, 0.4, 0.3, 0.5, # osc 2
               0.5, 0.1, 0.1, 0.1, # osc 3
              ]
        ),
]

solOptim3 = FitTMMOptics(
    Transmittance(), seed3, beam, Texp, layers;
    options=options, alg=SAMIN(), lb=0.1.*seed3, ub=2.0.*seed3,
)

plot(FitSpectrum(),
    solOptim3.Beam.λ, solOptim3.spectrumExp, solOptim3.spectrumFit,
    xaxis=("Wavelength [nm]"), yaxis=("Transmittance"),
)
gui()

# julia> solOptim3.objfunMin
# 0.00014874006735026815
#
# julia> solOptim3.optParams
# 1-element Array{Array{Float64,1},1}:
#  [192.75792016187577, 3.2411769190166146, 17.46281049095821, 6.567840560816103, 0.13772655091842556, 2.9673002141683855, 2.983451857126771, 0.3972297899281654, 0.5995812652718431, 0.9438184481128262, 0.9999875548110491, 0.19967018840936243, 0.1999535561385171, 0.010014114569596897]
