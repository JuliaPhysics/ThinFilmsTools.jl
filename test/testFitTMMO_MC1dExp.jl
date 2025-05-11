using Optim

# Wavelength range [nm]
λ = 400:1000
# Angle of incidence [degrees]
θ = [5.]
# Polarisation (p->1.0, s->0.0, weighted->between 0.0 and 1.0)
pol = 0.5
beam = PlaneWave(λ, θ; p=pol)

# Components of EMA
n1 = RIdb.air(beam.λ)
n2 = RIdb.silicon(beam.λ)
n3 = RIdb.glass(beam.λ)

# Define the RI model to use
layers = [
    LayerTMMO(n1), # 1
    ModelFit(:looyenga; N=(ninc=n1, nhost=n2)), # 2
    ModelFit(:looyenga; N=(ninc=n1, nhost=n2)), # 3
    LayerTMMO(n3), # 4
]

# Set the order of the layers (build the system) with the layers
order = [1, # incident medium
         2, 3, 2, 3, 2, 3, 2, 3, # top DBR of the MC
         2, 2, # defect of the MC
         3, 2, 3, 2, 3, 2, 3, 2, # bottom DBR of the MC
         4, # substrate
]

# Absolute reflectance spectrum to fit stored in Utils
Rexp = SpectraDB.fp_spectrum(beam.λ)

options = Optim.Options(
    g_abstol=1e-8, iterations=10^5, show_trace=true, store_trace=false,
    #g_reltol=1e-8,    # g_reltol no longer supported in Optim.jl
);

# Seeds for each ModelFit layer defined above (No alpha)
seed = [
    [119.0, 0.8], # layers[2]
    [76.0, 0.5], # layers[3]
]

solOptim = fit_tmm_optics(
    Reflectance(Rexp), seed, beam, layers;
    order=order, options=options, alg=SAMIN(), lb=0.5.*seed, ub=1.5.*seed,
)

# Seeds for each ModelFit layer defined above plus the alpha (with alpha)
seed2 = [
    [120, 0.8],
    [83.0, 0.5],
    [0.995], # alpha
]

solOptim2 = fit_tmm_optics(
    Reflectance(Rexp), seed2, beam, layers;
    order=order, options=options, alg=SAMIN(), lb=0.5.*seed2, ub=1.5.*seed2, alpha=true,
)

