# Load modules
using Plots, LaTeXStrings
pyplot(reuse=false, grid=false)
closeall()
using ThinFilmsTools

function main()
    # Define beam
    λ = LinRange(200, 1000, 600) # wavelength range [nm]
    θ = [0.] # angle of incidence [degrees]
    beam = PlaneWave(λ, θ)
    # Define layers
    l0 = LayerTMMO1DIso(RIdb.air(beam.λ))
    l1 = LayerTMMO1DIso(RIdb.dummy(beam.λ, 1.45, 0.); type=:OT, d=1/4.)
    l2 = LayerTMMO1DIso(RIdb.dummy(beam.λ, 3.45, 0.); type=:OT, d=1/4.)
    layers = [l0, l1, l2, l1, l2, l1, l2, l1, l2, l0]
    # Reference wavelenth [nm]
    λ0 = 700.
    # call main script
    return TMMO1DIsotropic(beam, layers; λ0=λ0, emfflag=true, h=10, pbgflag=true)
end

sol = main()

# plot the R, T and A spectra
plt1 = plot(TMMOPlotSpectra1D(),
            sol.Beam.λ, [sol.Spectra.Rp, sol.Spectra.Tp, 1.0.-(sol.Spectra.Rp.+sol.Spectra.Tp)], label=["Reflectance" "Transmittance" "Absorbance"], line=([:solid :dash :dashdot]), ylims=(0.0,1.0), xlims=(sol.Beam.λ[1], sol.Beam.λ[end]));
# gui()

# plot the EMF pattern for normal incidence
plt2 = plot(TMMOPlotEMF2D(),
            sol.Beam.λ, sol.Misc.ℓ, log10.(sol.Field.emfp[:,1,:]),
            title=("Log of EMF intensity"));
# gui()
plot(plt1, plt2, layout=(2,1))
gui()

# plot the refractive index profile
plot(TMMOPlotNprofile(), sol; plotemf=true)
gui()
### Note that if you choose θ outside the range it will show the EMF for one extrema. Same for λ.

# plot the photonic dispersion with custom function
plot(TMMOPlotDispersion1D(), sol.Bloch)
gui()

plot(TMMOPlotDispersion1Dalt(), sol.Bloch)
gui()

plot(TMMOPlotDispersion1Dimre(), sol.Bloch)
gui()
