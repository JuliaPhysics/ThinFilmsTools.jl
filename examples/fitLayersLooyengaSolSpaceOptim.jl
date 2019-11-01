# Load modules
using Plots, LaTeXStrings
pyplot()
using ThinFilmsTools
using Optim

##
function get_reflectance(ftype, λ, incident, emergent)
    # Raw measured spectrum stored in Utils
    Rexp = SpectraDB.SL2ExpSpectrum(beam.λ)
    # Reference measured spectrum stored in Utils
    Rref = SpectraDB.SL2RefSpectrum(beam.λ)
    # Theoretical reflectance spectrum for the reference
    Rthe = TheoreticalSpectrum(ftype, beam, incident, emergent)
    # Calculate the absolute normalised measured spectra to fit
    Rexp_norm = NormalizeReflectance(beam.λ, [beam.λ Rexp], [beam.λ Rthe], [beam.λ Rref])
end
##

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
layers = [
    LayerTMMO(incident), # 1
    ModelFit(:looyenga; N=(ninc=incident, nhost=emergent)), # 2
    LayerTMMO(emergent), # 3
]

# Get the spectrum to fit
Rexp = get_reflectance(Reflectance(), beam.λ, incident, emergent)

## Seed for the optimization algorithm

# Boundaries for the space solution and number of grid points
b = BoundariesFit(7000.0, 7300.0, 0.5, 0.65)

# Brute force search
sol = SpaceSolutionEMA(Reflectance(Rexp), b, beam, layers)

plot(SpaceSolution(),
    sol.od, sol.p, sol.solSpace,
    xaxis="Optical thickness [nm]",
    yaxis="Porosity"; num_levels=50)
gui()
plot(FitSpectrum(),
    sol.Beam.λ, sol.spectrumExp, sol.spectrumFit,
    xaxis="Wavelength [nm]",
    yaxis="Reflectance",
)
gui()

## Optimization using Optim

# Take the seed as the output from the 2D search
seed = vcat(sol.optThickness, sol.optParams)

options = Optim.Options(
    g_abstol=1e-8, g_reltol=1e-8, iterations=10^5, store_trace=true, show_trace=true,
)

solOptim = FitTMMOptics(
    ftype, [seed], beam, Rexp, layers;
    alg=SAMIN(), options=options, lb=[0.5.*seed], ub=[1.5.*seed],
)

plot(FitSpectrum(),
    solOptim.Beam.λ, solOptim.spectrumExp, solOptim.spectrumFit,
    xaxis="Wavelength [nm]",
    yaxis="Reflectance",
)
gui()
