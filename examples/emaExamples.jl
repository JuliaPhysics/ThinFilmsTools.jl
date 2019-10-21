# Load modules
using Plots, LaTeXStrings
pyplot()
using ThinFilmsTools

let
    λ = [632.8]
    n = [RIdb.air(λ) RIdb.silicon(λ)]
    p = LinRange(0,1,100)
    neff_looyenga = [RI.looyenga([p[i] 1.0.-p[i]], n)[1] for i in eachindex(p)]
    neff_bruggeman = [RI.bruggeman(p[i], n)[1] for i in eachindex(p)]
    neff_maxwell = [RI.maxwell([p[i] 1.0.-p[i]], n)[1] for i in eachindex(p)]
    neff_monecke = [RI.monecke(p[i], n)[1] for i in eachindex(p)]
    plot(p, real.(neff_looyenga), label="Looyenga", line=(:solid))
    plot!(p, real.(neff_bruggeman), label="Bruggeman", line=(:dash))
    plot!(p, real.(neff_maxwell), label="Maxwell", line=(:dashdot))
    plot!(p, real.(neff_monecke), label="Monecke", line=(2.5, :solid))
    xaxis!("Porosity, %")
    yaxis!("Effective dielectric permittivity")
    gui()
end
