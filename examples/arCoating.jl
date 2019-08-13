# Load modules
using Plots, LaTeXStrings
pyplot(reuse=false, grid=false)
closeall()
using ThinFilmsTools

# Define beam
λi = 400 # intial wavelength [nm]
λf = 1000 # final wavelength [nm]
λ = LinRange(λi, λf, λf-λi+1) # wavelength range [nm]
λ0 = 700. # reference wavelength
θ = [0.] # angle of incidence [degrees]
beam = PlaneWave(λ, λ0, θ)

# Define layers with their parameters
layers = [ LayerTMMO1DIso(type=:GT, n=RIdb.air(beam.λ), d=0.),
           LayerTMMO1DIso(type=:GT, n=DF.looyengaspheres(RIdb.air(beam.λ),RIdb.silicon(beam.λ),0.89), d=77.),
           LayerTMMO1DIso(type=:GT, n=DF.looyengaspheres(RIdb.air(beam.λ),RIdb.silicon(beam.λ),0.70), d=56.),
           LayerTMMO1DIso(type=:GT, n=DF.looyengaspheres(RIdb.air(beam.λ),RIdb.silicon(beam.λ),0.41), d=39.),
           LayerTMMO1DIso(type=:GT, n=RIdb.silicon(beam.λ), d=0.) ]

# call main script
sol = TMMO1DIsotropic(beam, layers; emfflag=true, h=10)

### Optional examples to plot results

# plot the R, T and A spectra
plot(TMMOplotSpectra1D(), sol.Beam.λ, [sol.Spectra.Rp, sol.Spectra.Tp, 1.0.-(sol.Spectra.Rp.+sol.Spectra.Tp)], label=["Reflectance" "Transmittance" "Absorbance"], line=([:solid :dash :dashdot]))
gui()

# plot the refractive index profile
TMMOplotNprofile(sol;)
### Note that if you choose θ outside the range it will show the EMF for one extrema. Same for λ.
