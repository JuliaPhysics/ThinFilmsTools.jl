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
λf = 1000 # final wavelength [nm]
λ = LinRange(λi, λf, λf-λi) # wavelength range [nm]
λ0 = 700. # reference wavelength
θ = [0.] # angle of incidence [degrees]
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
p1 = plot(TMMOplotSpectra1D(), sol.Beam.λ, [sol.Spectra.Rp, sol.Spectra.Tp, 1.0.-(sol.Spectra.Rp.+sol.Spectra.Tp)], label=["Reflectance" "Transmittance" "Absorbance"], line=([:solid :dash :dashdot]), ylims=(0.0,1.0), xlims=(sol.Beam.λ[1], sol.Beam.λ[end]));
plot(p1)
gui()

# plot the EMF pattern for normal incidence
plot(TMMOplotEMF2D(), sol.Beam.λ, sol.Misc.ℓ, log10.(sol.Field.emfp[:,1,:]), title=("Log of EMF intensity"))
gui()

# plot the refractive index profile
TMMOplotNprofile(sol;)
### Note that if you choose θ outside the range it will show the EMF for one extrema. Same for λ.

# plot the photonic dispersion with custom function
TMMOplotdispersion(sol)
