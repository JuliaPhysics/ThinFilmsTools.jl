# Load modules
using Plots, LaTeXStrings
pyplot(reuse=false, grid=false, size=(640,480))
closeall()
include("ThinFilmsTools.jl")
using Main.ThinFilmsTools

function nacl(x)
    x = x/1000
    n=@. sqrt(complex.(1+0.00055+0.19800./(1-(0.050./x).^2)+0.48398./(1-(0.100./x).^2)+0.38696./(1-(0.128./x).^2)+0.25998./(1-(0.158./x).^2)+0.08796./(1-(40.50./x).^2)+3.17064./(1-(60.98./x).^2)+0.30038./(1-(120.34./x).^2)))
    return real.(n)
end

# Define beam
λi = 4900 # intial wavelength [nm]
λf = 15000 # final wavelength [nm]
λ = LinRange(λi, λf, Integer(round((λf-λi+1)/7))) # wavelength range [nm]
λ0 = 12000. # reference wavelength
θi = 0
θf = 80
θ = LinRange(θi, θf, θf-θi+1) # angle of incidence [degrees]
beam = PlaneWave(λ, λ0, θ)

# Define layers
l0 = LayerTMMO1DIso(type=:GT, n=RIdb.air(beam.λ), d=0.)
l1 = LayerTMMO1DIso(type=:GT, n=RIdb.dummy(beam.λ, 4.6, 0.), d=800.)
l2 = LayerTMMO1DIso(type=:GT, n=RIdb.dummy(beam.λ, 1.6, 0.), d=1650.)
l3 = LayerTMMO1DIso(type=:GT, n=nacl(beam.λ), d=0.)
layers = [l0, l1, l2, l1, l2, l1, l2, l1, l2, l1, l3]

# call main script
sol = TMMO1DIsotropic(beam, layers; pbgflag=true)

### Optional examples to plot results

# plot the R, T and A spectra
t45 = findmin(abs.(sol.Beam.θ .- 45.0))[2][1]
t80 = findmin(abs.(sol.Beam.θ .- 80.0))[2][1]
p1 = plot(TMMOplotSpectra1D(), sol.Beam.λ, [sol.Spectra.Rp[:,1,:], sol.Spectra.Rs[:,1,:]], label=["p-wave" "s-wave"], yaxis=("Reflectance", (0., 1.)), xlims=(sol.Beam.λ[1], sol.Beam.λ[end]), title=(L"$\theta$ = 0 [$\degree$]"));
p2 = plot(TMMOplotSpectra1D(), sol.Beam.λ, [sol.Spectra.Rp[:,t45,:], sol.Spectra.Rs[:,t45,:]], label=["p-wave" "s-wave"], yaxis=("Reflectance", (0., 1.)), xlims=(sol.Beam.λ[1], sol.Beam.λ[end]), title=(L"$\theta$ = 45 [$\degree$]"));
p3 = plot(TMMOplotSpectra1D(), sol.Beam.λ, [sol.Spectra.Rp[:,t80,:], sol.Spectra.Rs[:,t80,:]], label=["p-wave" "s-wave"], yaxis=("Reflectance", (0., 1.)), xlims=(sol.Beam.λ[1], sol.Beam.λ[end]), title=(L"$\theta$ = 80 [$\degree$]"));
plot(p1, p2, p3, layout=(3,1))
gui()

plot(TMMOplotSpectra2D(), sol.Beam.λ, sol.Beam.θ, sol.Spectra.Rs, title=("s-wave Reflectance"))
gui()

# plot the refractive index profile
TMMOplotNprofile(sol)

# plot the photonic dispersion with custom function
TMMOplotdispersion(sol)
