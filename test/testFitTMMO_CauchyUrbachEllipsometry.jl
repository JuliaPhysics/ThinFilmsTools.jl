using Optim

function glassDeltaPsi(beam)
    layers2 = [ LayerTMMO1DIso(incident),
                LayerTMMO1DIso(RIdb.glass(beam.λ./1e3); d=400.),
                LayerTMMO1DIso(emergent) ]
    sol = TMMO1DIsotropic(beam, layers2)
    ratio = sol.Spectra.ρp./sol.Spectra.ρs
    return atan.(abs.(ratio)), angle.(ratio)
end

λ = 400:1000
θ = [60.]
pol = 0.5
beam = PlaneWave(λ, θ; p=pol)

incident = RIdb.air(beam.λ)
emergent = RIdb.silicon(beam.λ)

layers = [ LayerTMMO1DIso(incident),
           ModelFit(:cauchyurbach),
           LayerTMMO1DIso(emergent) ]

Ψ, Δ = glassDeltaPsi(beam)

options = Optim.Options(g_abstol=1e-8, g_reltol=1e-8, iterations=10^5, show_trace=true, store_trace=true);
seed = [vcat(420.0, [1.4, 0.23, 1.0, 0.0, 1.0, 1.0])]
solOptim = FitTMMO1DIsotropic(Ellipsometry(), seed, beam, hcat(Ψ, Δ), layers; options=options, alg=SAMIN(), lb=0.5.*seed, ub=1.5.*seed)
solOptim2 = FitTMMO1DIsotropic(Ellipsometry(), seed, beam, hcat(Ψ, Δ), layers; options=options, alg=NelderMead())
