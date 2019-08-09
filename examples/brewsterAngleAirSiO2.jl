# https://www.rp-photonics.com/brewsters_angle.html

# Load modules
using Plots, LaTeXStrings
pyplot(reuse=false, grid=false)
closeall()
include("ThinFilmsTools.jl")
using Main.ThinFilmsTools
include("RIdb.jl")
using Main.RIdb: air, dummy, bk7

# Define beam
λ = [1064.]
λ0 = 510.
θi = 0.01
θf = 90
θ = LinRange(θi, θf, 901)
beam = PlaneWave(λ, λ0, θ)

# Define layers: notice that lambda is outside the range
layers = [ LayerTMMO(type=:GT, n=dummy(beam.λ, 1., 0.), d=0., nλ0=dummy([beam.λ0], 1., 0.)),
           LayerTMMO(type=:GT, n=dummy(beam.λ, 1.5, 0.), d=10., nλ0=dummy([beam.λ0], 1.5, 0.)),
           LayerTMMO(type=:GT, n=dummy(beam.λ, 1.5, 0.), d=0., nλ0=dummy([beam.λ0], 1.5, 0.)) ]

# call main script
sol = TMMOptics(beam, layers)

### Optional examples to plot results

# plot the R, T and A spectra
p1 = plot(TMMOplotSpectraAngle1D(), sol.Beam.θ, [sol.Spectra.Rp[1,:], sol.Spectra.Rs[1,:], sol.Spectra.Rs[1,:]./sol.Spectra.Rp[1,:]./1000.], label=["p-wave" "s-wave" "(Rs/Rp)/1000"], line=([:solid :dash :dashdot]), xlims=(sol.Beam.θ[1], sol.Beam.θ[end]), yaxis=("Transmittance", (0.,0.2)));
plot(p1)
gui()

# plot the refractive index profile
TMMOplotNprofile(sol)
