# Load modules
using Plots, LaTeXStrings
pyplot(reuse=false, grid=false)
closeall()
using ThinFilmsTools

function main()
    # Define beam
    λ = LinRange(400, 700, 500) # wavelength range [nm]
    θ = [0.]
    beam = PlaneWave(λ, θ)
    # Define layers
    layers1 = [ LayerTMMO1DIso(RIdb.air(beam.λ)),
                LayerTMMO1DIso(RIdb.dummy(beam.λ, 1.38, 0.); type=:OT, d=1/4.),
                LayerTMMO1DIso(RIdb.dummy(beam.λ, 2.15, 0.); type=:OT, d=1/2.),
                LayerTMMO1DIso(RIdb.bk7(beam.λ)) ]
    layers2 = [ LayerTMMO1DIso(RIdb.air(beam.λ)),
                LayerTMMO1DIso(RIdb.dummy(beam.λ, 1.38, 0.); type=:OT, d=1/4.),
                LayerTMMO1DIso(RIdb.dummy(beam.λ, 1.9, 0.); type=:OT, d=1/2.),
                LayerTMMO1DIso(RIdb.bk7(beam.λ)) ]
    # Reference wavelenth [nm]
    λ0 = 509.8
    # call main script
    sol1 = TMMO1DIsotropic(beam, layers1; λ0=λ0)
    sol2 = TMMO1DIsotropic(beam, layers2; λ0=λ0)
    return sol1, sol2
end

sol1, sol2 = main()

# plot the R, T and A spectra
plot(TMMOPlotSpectra1D(), sol1.Beam.λ, sol1.Spectra.Rs,
     label=L"Ta$_2$O$_5$", line=(:solid),
     ylims=(0.0, 0.025), xlims=(sol1.Beam.λ[1], sol1.Beam.λ[end]),
     yaxis=("Reflectance"))
plot!(TMMOPlotSpectra1D(),
      sol2.Beam.λ, sol2.Spectra.Rs, label=L"Y$_2$O$_3$", line=(:dashdot))
gui()

# plot the refractive index profile
TMMOPlotNprofile(sol1)
### Note that if you choose θ outside the range it will show the EMF for one extrema. Same for λ.
# title!(L"Ta$_2$O$_5$")

TMMOPlotNprofile(sol2)
### Note that if you choose θ outside the range it will show the EMF for one extrema. Same for λ.
# title!(L"Y$_2$O$_3$")
