# Load modules
using Plots, LaTeXStrings
pyplot()
using ThinFilmsTools

function main()
    # Define beam
    λ = LinRange(400,700,500) # wavelength range [nm]
    θ = [0.] # Angle of incidence
    beam = PlaneWave(λ, θ)
    # Define layers for two systems
    layers1 = [
        LayerTMMO(RIdb.air(beam.λ)),
        LayerTMMO(RIdb.dummy(beam.λ,1.38,0.0); type=:OT, d=1/4.),
        LayerTMMO(RIdb.dummy(beam.λ,2.15,0.0); type=:OT, d=1/2.),
        LayerTMMO(RIdb.bk7(beam.λ)),
    ]
    layers2 = [
        LayerTMMO(RIdb.air(beam.λ)),
        LayerTMMO(RIdb.dummy(beam.λ,1.38,0.0); type=:OT, d=1/4.),
        LayerTMMO(RIdb.dummy(beam.λ,1.90,0.0); type=:OT, d=1/2.),
        LayerTMMO(RIdb.bk7(beam.λ)),
    ]
    # Reference wavelenth [nm]
    λ0 = 509.8
    # call main script
    sol1 = tmm_optics(beam, layers1; λ0=λ0)
    sol2 = tmm_optics(beam, layers2; λ0=λ0)
    return sol1, sol2
end

sol1, sol2 = main()

# plot the R, T and A spectra
plot(Spectrum1D(),
     sol1.beam.λ, sol1.Spectra.Rs,
     label=L"Ta$_2$O$_5$", line=(:solid),
     ylims=(0.0, 0.025),
     xlims=(sol1.Beam.λ[1], sol1.Beam.λ[end]),
     yaxis=("Reflectance"),
);
plot!(Spectrum1D(),
      sol2.beam.λ, sol2.Spectra.Rs,
      label=L"Y$_2$O$_3$",
      line=(:dashdot),
);
gui()

# plot the refractive index profile
plot(RIprofile(), sol1, title=L"Ta$_2$O$_5$")
gui()
### Note that if you choose θ outside the range it will show the EMF for one extrema. Same for λ.

plot(RIprofile(), sol2, title=L"Y$_2$O$_3$")
gui()
### Note that if you choose θ outside the range it will show the EMF for one extrema. Same for λ.
