# Load modules
#using Plots, LaTeXStrings
#pyplot()
#using ThinFilmsTools

# Load modules
using Plots, LaTeXStrings
pyplot()
using PyCall; pygui(true)
closeall()
include("/home/leniac/JuliaLangDev/ThinFilmsTools/src14/ThinFilmsTools.jl")
using Main.ThinFilmsTools

let
    λ = [632.8]
    n = [RIdb.air(λ) RIdb.silicon(λ)]
    p = LinRange(0,1,100)
    neff_looyenga = [RI.looyenga([p[i] 1.0-p[i]], n)[1] for i in eachindex(p)]
    neff_bruggeman = [RI.bruggeman(p[i], n)[1] for i in eachindex(p)]
    neff_maxwell = [RI.maxwell([p[i] 1.0-p[i]], n)[1] for i in eachindex(p)]
    neff_monecke = [RI.monecke(p[i], n)[1] for i in eachindex(p)]
    plot(p, real.(neff_looyenga), label="Looyenga", line=(:solid))
    plot!(p, real.(neff_bruggeman), label="Bruggeman", line=(:dash))
    plot!(p, real.(neff_maxwell), label="Maxwell", line=(:dashdot))
    plot!(p, real.(neff_monecke), label="Monecke", line=(3.0, :solid))
    xaxis!("Porosity, %")
    yaxis!("Effective dielectric permittivity")
    gui()
end


## lorentzlorenz EMA
# D. Schwarz et al. / Thin Solid Films 519 (2011) 2994-2997
let
    plt = plot();
    f = [0.0,0.25,0.5,0.75,1.0]
    ϵg = complex.(1.0:0.005:1.8)
    ϵₘ = 2.1*ones(length(ϵg))
    n = sqrt.([ϵg ϵₘ])
    for i in eachindex(f)
        neff = RI.lorentzlorenz(f[i], n)
        plot!(plt, real.(ϵg), real.(neff.^2), label="f = $(f[i])")
    end
    xaxis!("Guest dielectric constant, ϵg")
    yaxis!("Lorentz-Lorenz effective permittivity")
    gui()
end
