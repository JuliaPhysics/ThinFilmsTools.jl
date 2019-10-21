# Load modules
using Plots, LaTeXStrings
pyplot()
using ThinFilmsTools
using Optim

##
function glass_transmittance(beam, incident, emergent)
    layers = [
        LayerTMMO(incident),
        LayerTMMO(RI.looyenga([0.5 0.5],[RIdb.air(beam.λ) RIdb.silicon(beam.λ)]); d=100.),
        LayerTMMO(RIdb.glass(beam.λ./1e3); d=300.),
        LayerTMMO(RIdb.glass(beam.λ./1e3); d=100.),
        LayerTMMO(emergent),
    ]
    sol = TMMOptics(beam, layers)[1]
    return beam.p.*sol.Tp .+ (1.0 - beam.p).*sol.Ts
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
    ModelFit(:looyenga; N=(incident, emergent)), # 2
    LayerTMMO(RIdb.glass(beam.λ./1e3); d=300.), # 3
    ModelFit(:cauchyurbach), # 4
    LayerTMMO(emergent), # 5
]

# Create transmittance spectrum to fit
Texp = glass_transmittance(beam, incident, emergent)
# plot(Spectrum1D(), beam.λ, Texp)
# gui()

options = Optim.Options(
    g_abstol=1e-8, g_reltol=1e-8, iterations=10^5, show_trace=true, store_trace=true,
);

# Seeds of the ModelFit that appears in layers. The order here matters and the elements in each seed must match the number of elements accepted in the ModelFit layer used
seed = [
    [100.0, 0.5], # intial guess for :looyenga model, layers[2]
    [100.0, 1.5, 0.23, 1.0, 0.0, 1.0, 1.0], # intial guess for :cauchyurbach model, layers[4]
]

solOptim = FitTMMOptics(
    Transmittance(), seed, beam, Texp, layers;
    options=options, alg=SAMIN(), lb=0.5.*seed, ub=1.5.*seed,
)

plot(FitSpectrum(),
    solOptim.Beam.λ, solOptim.spectrumExp, solOptim.spectrumFit,
    xaxis=("Wavelength [nm]"),
    yaxis=("Transmittance"),
)
gui()

# julia> solOptim.objfunMin
# 2.0852427777269855e-7
#
# julia> solOptim.optParams
# 2-element Array{Array{Float64,1},1}:
#  [100.20261188423837, 0.5009380743144798]
#  [100.01313444247826, 1.5059167041688846, 0.1150000000000001, 0.5000000000000016, 0.0, 0.6175516492383015, 1.23581135410979]
