beam = PlaneWave(200:1100, [0.0])

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
Texp = SpectraDB.hafnia_spectrum(beam.λ)

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
solOptim = fit_tmm_optics(
    Transmittance(Texp), seed, beam, layers;
    alg=:BBO,
    SearchRange=Utils.unfoldbnd(lb,ub),
    Method=:de_rand_1_bin,
)

