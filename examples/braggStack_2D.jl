# Load modules
using Plots, LaTeXStrings
pyplot(reuse=false, grid=false)
closeall()
using ThinFilmsTools

# Define beam
λi = 200 # intial wavelength [nm]
λf = 5000 # final wavelength [nm]
λ = LinRange(λi, λf, λf-λi) # wavelength range [nm]
λ0 = 700. # reference wavelength
θi = 0
θf = 90
θ = LinRange(0, 90, θf-θi)
beam = PlaneWave(λ, λ0, θ)

# Define layers
l0 = LayerTMMO1DIso(type=:GT, n=RIdb.air(beam.λ), d=0.)
l1 = LayerTMMO1DIso(type=:OT, n=RIdb.dummy(beam.λ, 1.45, 0.), d=1/4.)
l2 = LayerTMMO1DIso(type=:OT, n=RIdb.dummy(beam.λ, 3.45, 0.), d=1/4.)

layers = [l0, l1, l2, l1, l2, l1, l2, l1, l2, l0]

# call main script
sol = TMMO1DIsotropic(beam, layers; emfflag=true, h=10, pbgflag=true)

### Optional examples to plot results

# plot the R, T and A spectra
# plot(TMMOplotSpectra2D(), sol.Beam.λ, sol.Beam.θ, sol.Spectra.Rp)
# gui()
#
# plot(TMMOplotSpectra2D(), sol.Beam.λ, sol.Beam.θ, sol.Spectra.Tp)
# gui()
#
# # plot the EMF pattern for normal incidence
# plot(TMMOplotEMF2D(), sol.Beam.λ, sol.Misc.ℓ, log10.(sol.Field.emfp[:,1,:]))
# gui()
#
# # plot the refractive index profile
# TMMOplotNprofile(sol)

# plot the photonic dispersion with custom function
TMMOplotdispersion(sol)
