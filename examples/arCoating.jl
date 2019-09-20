# Load modules
using Plots, LaTeXStrings
pyplot(reuse=false, grid=false)
using ThinFilmsTools

# Define beam
λi = 400 # intial wavelength [nm]
λf = 1000 # final wavelength [nm]
λ = LinRange(λi, λf, λf-λi+1) # wavelength range [nm]
θ = [0.] # angle of incidence [degrees]
beam = PlaneWave(λ, θ)

# Define layers with their parameters
layers = [ LayerTMMO1DIso(RIdb.air(beam.λ)),
           LayerTMMO1DIso(RI.looyengaspheresbin(0.89,RIdb.air(beam.λ),RIdb.silicon(beam.λ)); d=77.),
           LayerTMMO1DIso(RI.looyengaspheresbin(0.70,RIdb.air(beam.λ),RIdb.silicon(beam.λ)); d=56.),
           LayerTMMO1DIso(RI.looyengaspheresbin(0.41,RIdb.air(beam.λ),RIdb.silicon(beam.λ)), d=39.),
           LayerTMMO1DIso(RIdb.silicon(beam.λ)) ]

# Reference wavelength [nm]
λ0 = 700.0

# call main script
sol = TMMO1DIsotropic(beam, layers)

# plot the R, T and A spectra
plot(TMMOPlotSpectra1D(), sol.Beam.λ, [sol.Spectra.Rp, sol.Spectra.Tp, 1.0 .- (sol.Spectra.Rp .+ sol.Spectra.Tp)], label=["Reflectance" "Transmittance" "Absorbance"], line=([:solid :dash :dashdot]))
gui()

# plot the refractive index profile
TMMOPlotNprofile(sol;)

