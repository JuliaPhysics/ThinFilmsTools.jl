# Load modules
using Plots, LaTeXStrings
pyplot()
using ThinFilmsTools
using Optim

# Wavelength range [nm] and angle of incidence [degrees]
beam = PlaneWave(200:2100, [60.])

# Refractive indices of media
incident = RIdb.air(beam.λ)
emergent = RIdb.fused_silica_uv(beam.λ)

# Define the RI model to use
layers = [
    LayerTMMO(incident),
    ModelFit(:tauclorentz),
    LayerTMMO(emergent),
]

# Data
Ψexp, Δexp = SpectraDB.tantala_spectrum(beam.λ)
plot(1240.0./beam.λ, Ψexp)
plot!(1240.0./beam.λ, Δexp)
gui()

## Search
seed = [vcat(92.0, # thickness
             [1.5, 4.3, # ϵinf, Eg
              340, 5.0, 2.8, # osc 1, A1, E01, Γ1
              715, 2, 0.5, # osc 2
              41, 5, 0.8, # osc 3
             ],
        ),
]

# Spectrum type: ellipsometry with the input experimental spectra
spectype = Ellipsometry([Ψexp Δexp])

## Optimization using Optim
options = Optim.Options(
    g_abstol=1e-8, g_reltol=1e-8, iterations=10^5, show_trace=true, store_trace=true,
);

solOptim0 = fit_tmm_optics(
    spectype, seed, beam, layers;
    options=options, alg = LBFGS(), oftype=SumMeanAbs(),
)

plot(FitSpectrumEllip(),
    solOptim0.beam.λ, solOptim0.spectrumExp, solOptim0.spectrumFit,
    xaxis=("Wavelength [nm]"),
)
gui()
# plot(FitSpectrum(),
#     solOptim0.Beam.λ, solOptim0.spectrumExp, solOptim0.spectrumFit,
#     xaxis=("Wavelength [nm]"),
# )
# gui()
# julia> solOptim0.optParams
# 1-element Array{Array{Float64,1},1}:
#  [91.9939987559028, 1.4024927693924976, 4.359976200010805, 339.99868736615565, 5.093422704869424, 2.642028933308104, 715.0005849269817, 1.9891814175404072, 0.7171551391581684, 41.000188071713396, 4.94789281453286, 0.7311345356462031]
#
# julia> solOptim0.objfunMin
# 0.015678263617751792

lb = 0.5.*solOptim0.optParams
ub = 2.0.*solOptim0.optParams
# Thickness
lb[1][1] = 0.9.*solOptim0.optParams[1][1]
ub[1][1] = 1.2.*solOptim0.optParams[1][1]
# Egap
lb[1][3] = 3.0
ub[1][3] = 6.0
# solOptim = FitTMMOptics(
    # Ellipsometry(), solOptim0.optParams, beam, hcat(Ψexp, Δexp), layers;
    # options=options, alg=SAMIN(), lb=lb, ub=ub,
# )
solOptim = fit_tmm_optics(
    spectype, solOptim0.optParams, beam, layers;
    options=options, alg=SAMIN(), lb=lb, ub=ub, oftype=SumMeanAbs(),
)
plot(FitSpectrumEllip(),
    solOptim.beam.λ, solOptim.spectrumExp, solOptim.spectrumFit,
    xaxis=("Wavelength [nm]"),
)
gui()
# plot(FitSpectrum(),
#     solOptim.Beam.λ, solOptim.spectrumExp, solOptim.spectrumFit,
#     xaxis=("Wavelength [nm]"),
# )
# gui()

# julia> solOptim.objfunMin
# 0.01011540926111779
#
# julia> solOptim.optParams
# 1-element Array{Array{Float64,1},1}:
#  [92.41268518127337, 1.294018564396567, 4.258287813517591, 348.0769327939244, 3.7171269504011133, 4.97982821966445, 574.0431000488865, 3.7406423898551893, 0.580618765881699, 81.73269057321525, 5.273106528658877, 1.0466795973678233]
