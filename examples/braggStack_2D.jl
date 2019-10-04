# Load modules
using Plots, LaTeXStrings
pyplot(reuse=false, grid=false)
closeall()
using ThinFilmsTools

function main()
    # Define beam
    λ = LinRange(200, 5000, 1000) # wavelength range [nm]
    θ = LinRange(0, 90, 90)
    beam = PlaneWave(λ, θ)
    # Define layers
    l0 = LayerTMMO1DIso(RIdb.air(beam.λ))
    l1 = LayerTMMO1DIso(RIdb.dummy(beam.λ, 1.45, 0.); type=:OT, d=1/4.)
    l2 = LayerTMMO1DIso(RIdb.dummy(beam.λ, 3.45, 0.); type=:OT, d=1/4.)
    layers = [l0, l1, l2, l1, l2, l1, l2, l1, l2, l0]
    # Reference wavelenth
    λ0 = 700.
    # call main script
    return TMMO1DIsotropic(beam, layers; λ0=λ0, emfflag=true, h=10, pbgflag=true)
end

function plotEMF(sol, ϕ)
    ϕ_ = findmin(abs.(sol.Beam.θ .- ϕ))[2][1]
    plot(
         plot(TMMOPlotEMF2D(),
              sol.Beam.λ, sol.Misc.ℓ, log10.(sol.Field.emfp[:,ϕ_,:]),
              title="Log EMF intesnsity, p-TM"),
         plot(TMMOPlotEMF2D(),
              sol.Beam.λ, sol.Misc.ℓ, log10.(sol.Field.emfs[:,ϕ_,:]),
              title="Log EMF intesnsity, s-TE"),
         layout=(2,1))
    gui()
    return nothing
end

sol = main()

# plot the R, T and A spectra
plt1 = plot(TMMOPlotSpectra2D(),
            sol.Beam.λ, sol.Beam.θ, sol.Spectra.Rp,
            title="Reflectance, p-TM", clims=(0.0, 1.0));
# gui()

plt2 = plot(TMMOPlotSpectra2D(),
            sol.Beam.λ, sol.Beam.θ, sol.Spectra.Rs,
            title="Reflectance, s-TE", clims=(0.0, 1.0));
# gui()
plot(plt1, plt2, layout=(2,1))
gui()

# Plot the EMF at normal incidence
plotEMF(sol, 0.0)

# Plot the EMF at 15 degrees
plotEMF(sol, 15.0)

# plot the refractive index profile
plot(TMMOPlotNprofile(), sol; plotemf=true)
gui()
### Note that if you choose θ outside the range it will show the EMF for one extrema. Same for λ.

# plot the photonic dispersion with custom function
plot(TMMOPlotDispersion2D(), sol.Bloch; wave=:p)
gui()

plot(TMMOPlotDispersion2D(), sol.Bloch; wave=:s)
gui()

plot(TMMOPlotDispersion2Dalt(), sol.Bloch)
gui()
