using Plots
pyplot(reuse=false, size=(640, 480), grid=false)
closeall()
using Optim
using ThinFilmsTools

# Type of fitting
ftype = Reflectance()

# Wavelength range [nm]
λ = 250:900
# Angle of incidence [degrees]
θ = [5.]
# Polarisation (1.0 = p, 0.0 = s, between 0.0 and 1.0 = average)
pol = 0.5
beam = PlaneWave(λ, θ; p=pol)

# Refractive indices of incident (0) and substrate (2)
incident = RIdb.air(beam.λ)
emergent = RIdb.silicon(beam.λ)

# Define the RI model to use
layers = [ LayerTMMO1DIso(incident),
           ModelFit(:looyengaspheresbin; N=(incident, emergent)),
           LayerTMMO1DIso(emergent) ]

# Raw measured spectrum stored in Utils
Rexp = Utils.getSL2ExpSpectrum(beam.λ)
# Reference measured spectrum stored in Utils
Rref = Utils.getSL2RefSpectrum(beam.λ)
# Theoretical reflectance spectrum for the reference
Rthe = TheoreticalSpectrum(ftype, beam, incident, emergent)
# Calculate the absolute normalised measured spectra to fit 
Rexp_norm = NormalizeReflectance(beam.λ, [beam.λ Rexp], [beam.λ Rthe], [beam.λ Rref])

## Seed for the optimization algorithm

# Boundaries for the space solution and number of grid points
b = BoundariesFit(7000.0, 7300.0, 0.5, 0.65)

# Brute force search
sol = SpaceSolution2D(ftype, b, beam, Rexp_norm, layers)

plot(SpaceSolutionOFplot(), sol.od, sol.x, sol.solSpace, xaxis=("Optical thickness [nm]"), yaxis=("Porosity"); num_levels=50)
gui()
plot(PlotFitSpectrum(), sol.beam.λ, sol.spectrumExp, sol.spectrumFit, xaxis=("Wavelength [nm]"), yaxis=("Reflectance"))
gui()

## Optimization using Optim

# Take the seed as the output from the 2D search
seed = vcat(sol.optThickness, sol.optParams)

options = Optim.Options(g_abstol=1e-8, g_reltol=1e-8, iterations=10^5, store_trace=true, show_trace=true)

solOptim = FitTMMO1DIsotropic(ftype, [seed], beam, Rexp_norm, layers; alg=SAMIN(), options=options, lb=[0.8.*seed], ub=[1.2.*seed])
plot(PlotFitSpectrum(), solOptim.beam.λ, solOptim.spectrumExp, solOptim.spectrumFit)
gui()

