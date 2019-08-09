# Load modules
using Plots, LaTeXStrings
pyplot(reuse=false, grid=false)
closeall()
include("ThinFilmsTools.jl")
using Main.ThinFilmsTools
include("RIdb.jl") # collection of refractive indexes data
using Main.RIdb: air, silicon
include("MixingRules.jl") # collection of mixing rules for dielectric functions
using Main.MixingRules: looyengaspheres

# Define beam
λi = 400 # intial wavelength [nm]
λf = 1000 # final wavelength [nm]
λ = LinRange(λi, λf, λf-λi+1) # wavelength range [nm]
λ0 = 700. # reference wavelength
θ = [0.] # angle of incidence [degrees]
beam = PlaneWave(λ, λ0, θ)

# Define layers
l0 = LayerTMMO(type=:GT, n=air(beam.λ), d=0.)
l1 = LayerTMMO(type=:OT, n=looyengaspheres(air(beam.λ),silicon(beam.λ),0.54), d=1/4.)
l2 = LayerTMMO(type=:OT, n=looyengaspheres(air(beam.λ),silicon(beam.λ),0.86), d=1/4.)
l3 = LayerTMMO(type=:GT, n=silicon(beam.λ), d=0.)
layers = vec([l0 l1 l2 l1 l2 l1 l2 l1 l2 l2 l2 l2 l1 l2 l1 l2 l1 l2 l1 l3])

# call main script
sol = TMMOptics(beam, layers; emfflag=true, h=10)

### Optional examples to plot results

# plot the R, T and A spectra
plot(TMMOplotSpectra1D(), sol.Beam.λ, [sol.Spectra.Rp, sol.Spectra.Tp, 1.0.-(sol.Spectra.Rp.+sol.Spectra.Tp)], label=["Reflectance" "Transmittance" "Absorbance"], line=([:solid :dash :dashdot]), xlims=(sol.Beam.λ[1], sol.Beam.λ[end]), yaxis=("Transmittance", (0.,1.0)));
gui()

plot(TMMOplotEMF2D(), sol.Beam.λ, sol.Misc.ℓ, sol.Field.emfp[:,1,:])
gui()

# plot the refractive index profile
TMMOplotNprofile(sol)
