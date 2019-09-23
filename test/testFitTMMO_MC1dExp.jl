using Optim

λ = 400:1000
θ = [5.]
pol = 0.5
beam = PlaneWave(λ, θ; p=pol)

n1 = RIdb.air(beam.λ)
n2 = RIdb.silicon(beam.λ)
n3 = RIdb.glass(beam.λ)

layers = [ LayerTMMO1DIso(n1),
           ModelFit(:looyengaspheresbin; N=(n1, n2)),
           ModelFit(:looyengaspheresbin; N=(n1, n2)),
           LayerTMMO1DIso(n3) ]

Rexp = SpectraDB.FPSpectrum(beam.λ)

options = Optim.Options(g_abstol=1e-8, g_reltol=1e-8, iterations=10^6, show_trace=true, store_trace=true);

seed = [[119.0, 0.8], [76.0, 0.5]]
solOptim = FitTMMO1DIsotropic(Reflectance(), seed, beam, Rexp, layers; arrange=MC1d(), L=[4, 4], Ld=[2], options=options, alg=SAMIN(), lb=0.5.*seed, ub=1.5.*seed)
solOptim2 = FitTMMO1DIsotropic(Reflectance(), seed, beam, Rexp, layers; arrange=MC1d(), L=[4, 4], Ld=[2], options=options, alg=NelderMead())

seed2 = [[120, 0.8], [83.0, 0.5], [0.995]]
solOptim3 = FitTMMO1DIsotropic(Reflectance(), seed2, beam, Rexp, layers; arrange=MC1dAlpha(), L=[4, 4], Ld=[2], options=options, alg=SAMIN(), lb=0.5.*seed2, ub=1.5.*seed2)
solOptim4 = FitTMMO1DIsotropic(Reflectance(), seed2, beam, Rexp, layers; arrange=MC1dAlpha(), L=[4, 4], Ld=[2], options=options, alg=NelderMead())

