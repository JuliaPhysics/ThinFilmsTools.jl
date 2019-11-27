# Load modules
using Plots, LaTeXStrings
pyplot()
using ThinFilmsTools

function main()
    # Define beam
    λ = LinRange(200,5000,1000) # wavelength range [nm]
    θ = LinRange(0,90,90)
    beam = PlaneWave(λ,θ)
    # Define layers
    l0 = LayerTMMO(RIdb.air(beam.λ))
    l1 = LayerTMMO(RIdb.dummy(beam.λ,1.45,0.0); type=:OT, d=1/4.)
    l2 = LayerTMMO(RIdb.dummy(beam.λ,3.45,0.0); type=:OT, d=1/4.)
    layers = [l0, l1, l2, l1, l2, l1, l2, l1, l2, l0]
    # Reference wavelenth
    λ0 = 700.0
    # call main script
    sol = tmm_optics(beam, layers; λ0=λ0, emfflag=true, h=10, pbgflag=true)
    return sol
end

function plotEMF(sol, ϕ)
    ϕ_ = findmin(abs.(sol.beam.θ .- ϕ))[2][1]
    plot(EMF2D(),
              sol.beam.λ, sol.Misc.ℓ, log10.(sol.Field.emfp[:,ϕ_,:]),
              title="Log EMF intesnsity, p-TM",
    )
    gui()
    plot(EMF2D(),
              sol.beam.λ, sol.Misc.ℓ, log10.(sol.Field.emfs[:,ϕ_,:]),
              title="Log EMF intesnsity, s-TE",
    )
    gui()
    return nothing
end

sol = main()

# plot the R, T and A spectra
plot(Spectrum2D(),
        sol.beam.λ, sol.beam.θ, sol.Spectra.Rp,
        title="Reflectance, p-TM", clims=(0.0, 1.0),
)
gui()
plot(Spectrum2D(),
        sol.beam.λ, sol.beam.θ, sol.Spectra.Rs,
        title="Reflectance, s-TE", clims=(0.0, 1.0),
)
gui()

# Plot the EMF at normal incidence
plotEMF(sol, 0.0)

# Plot the EMF at 15 degrees
plotEMF(sol, 15.0)

# plot the refractive index profile
plot(RIprofile(), sol; plotemf=true)
gui()
### Note that if you choose θ outside the range it will show the EMF for one extrema. Same for λ.

# plot the photonic dispersion with custom function
plot(PBGDispersion2D(), sol.Bloch; wave=:p)
gui()

plot(PBGDispersion2D(), sol.Bloch; wave=:s)
gui()

plot(PBGDispersion2Dalt(), sol.Bloch)
gui()
