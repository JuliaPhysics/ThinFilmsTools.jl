using Optim

ftype = Reflectance()

λ = 250:900
θ = [5.]
pol = 0.5
beam = PlaneWave(λ, θ; p=pol)

incident = RIdb.air(beam.λ)
emergent = RIdb.silicon(beam.λ)

layers = [ LayerTMMO1DIso(incident),
           ModelFit(:looyengaspheresbin; N=(incident, emergent)),
           LayerTMMO1DIso(emergent) ]

Rexp = SpectraDB.SL1ExpSpectrum(beam.λ)
Rref = SpectraDB.SL1RefSpectrum(beam.λ)
Rthe = TheoreticalSpectrum(ftype, beam, incident, emergent)
Rexp_norm = NormalizeReflectance(beam.λ, [beam.λ Rexp], [beam.λ Rthe], [beam.λ Rref])

seed = vcat(3318, 0.85)
options = Optim.Options(g_abstol=1e-8, g_reltol=1e-8, iterations=10^5, store_trace=true, show_trace=true)
solOptim = FitTMMO1DIsotropic(ftype, [seed], beam, Rexp_norm, layers; alg=SAMIN(), options=options, lb=[0.5.*seed], ub=[1.5.*seed])
solOptim2 = FitTMMO1DIsotropic(ftype, [seed], beam, Rexp_norm, layers; alg=NelderMead(), options=options)
