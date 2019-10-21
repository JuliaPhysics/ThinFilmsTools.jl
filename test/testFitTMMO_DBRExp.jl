using Optim

λ = 400:1000
θ = [5.]
pol = 0.5
beam = PlaneWave(λ, θ; p=pol)

# Refractive indices of incident (0) and substrate (2)
incident = RIdb.air(beam.λ)
emergent = RIdb.silicon(beam.λ)

# Define the RI model to use
layers = [
    LayerTMMO(incident),
    ModelFit(:looyenga; N=(incident, emergent)),
    ModelFit(:looyenga; N=(incident, emergent)),
    LayerTMMO(emergent),
]

# Set the order of the layers (build the system) to fit
order = [
    1 # incident medium
    repeat([2, 3], 9) # top DBR of the MC
    4 # substrate
]

# Absolute reflectance spectrum to fit stored in Utils
Rexp = SpectraDB.BraggSpectrum(beam.λ)

options = Optim.Options(
    g_abstol=1e-8, g_reltol=1e-8, iterations=10^5, show_trace=true, store_trace=true,
);

# Seeds for each ModelFit layer defined above (without alpha)
seed = [
    [69, 0.54], # layers[2]
    [119.0, 0.77], # layers[3]
]

solOptim = FitTMMOptics(
    Reflectance(), seed, beam, Rexp, layers;
    order=order, options=options, alg=SAMIN(), lb=0.5.*seed, ub=1.5.*seed,
)

# Seeds for each ModelFit layer defined above plus the alpha (with alpha)
seed2 = [
    [69, 0.54], # layers[2]
    [119.0, 0.77], # layers[3]
    [0.995], # alpha
]

solOptim2 = FitTMMOptics(
    Reflectance(), seed2, beam, Rexp, layers;
    order=order, options=options, alg=SAMIN(), lb=0.5.*seed2, ub=1.5.*seed2, alpha=true,
)

