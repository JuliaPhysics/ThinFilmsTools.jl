using Optim

λ = 400:1000
θ = [5.]
pol = 0.5
beam = PlaneWave(λ, θ; p=pol)

incident = RIdb.air(beam.λ)
emergent = RIdb.silicon(beam.λ)

layers = [ LayerTMMO1DIso(incident),
           ModelFit(:looyengaspheresbin; N=(incident, emergent)),
           ModelFit(:looyengaspheresbin; N=(incident, emergent)),
           LayerTMMO1DIso(emergent) ]

Rexp = SpectraDB.BraggSpectrum(beam.λ)

options = Optim.Options(g_abstol=1e-8, g_reltol=1e-8, iterations=10^5, show_trace=true, store_trace=true);

seed = [[69., 0.54], [119.0, 0.77]]
solOptim = FitTMMO1DIsotropic(Reflectance(), seed, beam, Rexp, layers; arrange=DBR(), L=[9], options=options, alg=SAMIN(), lb=0.5.*seed, ub=1.5.*seed)

seed2 = [[69, 0.54], [119.0, 0.77], [0.995]]
solOptim2 = FitTMMO1DIsotropic(Reflectance(), seed2, beam, Rexp, layers; arrange=DBRAlpha(), L=[9], options=options, alg=SAMIN(), lb=0.5.*seed2, ub=1.5.*seed2)

