# Load modules
using Plots, LaTeXStrings
pyplot(reuse=false, grid=false)
closeall()
include("ThinFilmsTools.jl")
using Main.ThinFilmsTools
include("RIdb.jl") # collection of refractive indexes data
using Main.RIdb: air, silicon, dummy

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
l0 = LayerTMMO(type=:GT, n=air(beam.λ), d=0.)
l1 = LayerTMMO(type=:OT, n=dummy(beam.λ, 1.45, 0.), d=1/4., nλ0=dummy([beam.λ0], 1.45, 0.))
l2 = LayerTMMO(type=:OT, n=dummy(beam.λ, 3.45, 0.), d=1/4., nλ0=dummy([beam.λ0], 3.45, 0.))

layers = [l0, l1, l2, l1, l2, l1, l2, l1, l2, l0]

# call main script
sol = TMMOptics(beam, layers; emfflag=true, h=10, pbgflag=true)

### Optional examples to plot results

# plot the R, T and A spectra
plot(TMMOplotSpectra2D(), sol.Beam.λ, sol.Beam.θ, sol.Spectra.Rp)
gui()

plot(TMMOplotSpectra2D(), sol.Beam.λ, sol.Beam.θ, sol.Spectra.Tp)
gui()

# plot the EMF pattern for normal incidence
plot(TMMOplotEMF2D(), sol.Beam.λ, sol.Misc.ℓ, log10.(sol.Field.emfp[:,1,:]))
gui()

# plot the refractive index profile
TMMOplotNprofile(sol)

# plot the photonic dispersion with custom function
TMMOplotdispersion(sol)
