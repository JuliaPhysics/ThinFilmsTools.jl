# Load modules
using Plots, LaTeXStrings
pyplot(reuse=false, grid=false)
closeall()
using ThinFilmsTools

# Define beam
λi = 200 # intial wavelength [nm]
λf = 1000 # final wavelength [nm]
λ = LinRange(λi, λf, λf-λi) # wavelength range [nm]
θ = [0.] # angle of incidence [degrees]
beam = PlaneWave(λ, θ)

# Define layers
l0 = LayerTMMO1DIso(RIdb.air(beam.λ))
l1 = LayerTMMO1DIso(RIdb.dummy(beam.λ, 1.45, 0.); type=:OT, d=1/4.)
l2 = LayerTMMO1DIso(RIdb.dummy(beam.λ, 3.45, 0.); type=:OT, d=1/4.)
layers = [l0, l1, l2, l1, l2, l1, l2, l1, l2, l0]

# Reference wavelenth
λ0 = 700.

# call main script
sol = TMMO1DIsotropic(beam, layers; λ0=λ0, emfflag=true, h=10, pbgflag=true)

# plot the R, T and A spectra
plot(TMMOPlotSpectra1D(), sol.Beam.λ, [sol.Spectra.Rp, sol.Spectra.Tp, 1.0.-(sol.Spectra.Rp.+sol.Spectra.Tp)], label=["Reflectance" "Transmittance" "Absorbance"], line=([:solid :dash :dashdot]), ylims=(0.0,1.0), xlims=(sol.Beam.λ[1], sol.Beam.λ[end]));
gui()

# plot the EMF pattern for normal incidence
plot(TMMOPlotEMF2D(), sol.Beam.λ, sol.Misc.ℓ, log10.(sol.Field.emfp[:,1,:]), title=("Log of EMF intensity"))
gui()

# plot the refractive index profile
TMMOPlotNprofile(sol;)

# plot the photonic dispersion with custom function
TMMOPlotDispersion(sol)
