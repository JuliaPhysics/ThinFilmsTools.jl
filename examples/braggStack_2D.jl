# Load modules
using Plots, LaTeXStrings
pyplot(reuse=false, grid=false)
closeall()
using ThinFilmsTools

# Define beam
λi = 200 # intial wavelength [nm]
λf = 5000 # final wavelength [nm]
λ = LinRange(λi, λf, 1000) # wavelength range [nm]
λ0 = 700. # reference wavelength
θi = 0
θf = 90
θ = LinRange(0, 90, θf-θi)
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
plot(TMMOPlotSpectra2D(), sol.Beam.λ, sol.Beam.θ, sol.Spectra.Rp, title="Reflectance, p-TM", clims=(0.0, 1.0))
gui()

plot(TMMOPlotSpectra2D(), sol.Beam.λ, sol.Beam.θ, sol.Spectra.Rs, title="Reflectance, s-TE", clims=(0.0, 1.0))
gui()

# plot the EMF pattern for normal incidence and 30 degrees
plot(TMMOPlotEMF2D(), sol.Beam.λ, sol.Misc.ℓ, log10.(sol.Field.emfp[:,1,:]), title="Log EMF intesnsity, p-TM")
gui()

plot(TMMOPlotEMF2D(), sol.Beam.λ, sol.Misc.ℓ, log10.(sol.Field.emfs[:,1,:]), title="Log EMF intesnsity, s-TE")
gui()

aux1 = findmin(abs.(sol.Beam.θ .- 15.0))[2][1]
plot(TMMOPlotEMF2D(), sol.Beam.λ, sol.Misc.ℓ, log10.(sol.Field.emfp[:,aux1,:]), title="Log EMF intesnsity, p-TM, θ = 15 [°]")
gui()

plot(TMMOPlotEMF2D(), sol.Beam.λ, sol.Misc.ℓ, log10.(sol.Field.emfs[:,aux1,:]), title="Log EMF intesnsity, s-TE, θ = 15 [°]")
gui()

# plot the refractive index profile
TMMOPlotNprofile(sol)

# plot the photonic dispersion with custom function
TMMOPlotDispersion(sol)