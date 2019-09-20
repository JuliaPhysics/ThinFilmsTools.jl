# Load modules
using Plots, LaTeXStrings
pyplot(reuse=false, grid=false, size=(640,480))
closeall()
using ThinFilmsTools

# Substrate index of refraction taken from www.refractiveindex.info
function nacl(x)
    x = x/1000
    n=@. sqrt(complex.(1+0.00055+0.19800./(1-(0.050./x).^2)+0.48398./(1-(0.100./x).^2)+0.38696./(1-(0.128./x).^2)+0.25998./(1-(0.158./x).^2)+0.08796./(1-(40.50./x).^2)+3.17064./(1-(60.98./x).^2)+0.30038./(1-(120.34./x).^2)))
    return real.(n)
end

# Define beam
λi = 4900 # intial wavelength [nm]
λf = 15000 # final wavelength [nm]
λ = LinRange(λi, λf, 5000) # wavelength range [nm]
θi = 0
θf = 80
θ = LinRange(θi, θf, θf-θi+1) # angle of incidence [degrees]
beam = PlaneWave(λ, θ)

# Define layers
l0 = LayerTMMO1DIso(RIdb.air(beam.λ))
l1 = LayerTMMO1DIso(RIdb.dummy(beam.λ, 4.6, 0.); d=800.)
l2 = LayerTMMO1DIso(RIdb.dummy(beam.λ, 1.6, 0.); d=1650.)
l3 = LayerTMMO1DIso(nacl(beam.λ))
layers = [l0, l1, l2, l1, l2, l1, l2, l1, l2, l1, l3]

# Reference wavelength
λ0 = 12000.

# call main script
sol = TMMO1DIsotropic(beam, layers; λ0=λ0, pbgflag=true)

# plot the R, T and A spectra
t45 = findmin(abs.(sol.Beam.θ .- 45.0))[2][1]
t80 = findmin(abs.(sol.Beam.θ .- 80.0))[2][1]
p1 = plot(TMMOPlotSpectra1D(), sol.Beam.λ, [sol.Spectra.Rp[:,1,:], sol.Spectra.Rs[:,1,:]], label=["p-wave" "s-wave"], yaxis=("Reflectance", (0., 1.)), xlims=(sol.Beam.λ[1], sol.Beam.λ[end]), title=(L"$\theta$ = 0 [$\degree$]"), line=([:solid :dashdot]));
p2 = plot(TMMOPlotSpectra1D(), sol.Beam.λ, [sol.Spectra.Rp[:,t45,:], sol.Spectra.Rs[:,t45,:]], label=["p-wave" "s-wave"], yaxis=("Reflectance", (0., 1.)), xlims=(sol.Beam.λ[1], sol.Beam.λ[end]), title=(L"$\theta$ = 45 [$\degree$]"), line=([:solid :dashdot]));
p3 = plot(TMMOPlotSpectra1D(), sol.Beam.λ, [sol.Spectra.Rp[:,t80,:], sol.Spectra.Rs[:,t80,:]], label=["p-wave" "s-wave"], yaxis=("Reflectance", (0., 1.)), xlims=(sol.Beam.λ[1], sol.Beam.λ[end]), title=(L"$\theta$ = 80 [$\degree$]"), line=([:solid :dashdot]));
plot(p1, p2, p3, layout=(3,1))
gui()

plot(TMMOPlotSpectra2D(), sol.Beam.λ, sol.Beam.θ, sol.Spectra.Rs, title=("s-wave Reflectance"))
gui()

plot(TMMOPlotSpectra2D(), sol.Beam.λ, sol.Beam.θ, sol.Spectra.Rp, title=("p-wave Reflectance"))
gui()

# plot the refractive index profile
TMMOPlotNprofile(sol)

# plot the photonic dispersion with custom function
TMMOPlotDispersion(sol)
