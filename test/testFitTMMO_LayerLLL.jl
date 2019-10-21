using Optim

##
function get_reflectance(ftype, λ, incident, emergent)
    # Raw measured spectrum stored in Utils
    Rexp = SpectraDB.SL1ExpSpectrum(beam.λ)
    # Reference measured spectrum stored in Utils
    Rref = SpectraDB.SL1RefSpectrum(beam.λ)
    # Theoretical reflectance spectrum for the reference
    Rthe = TheoreticalSpectrum(ftype, beam, incident, emergent)
    # Calculate the absolute normalised measured spectra to fit
    Rexp_norm = NormalizeReflectance(beam.λ, [beam.λ Rexp], [beam.λ Rthe], [beam.λ Rref])
end
##

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
layers = [
    LayerTMMO(incident), # 1
    ModelFit(:looyenga; N=(incident, emergent)), # 2
    LayerTMMO(emergent), # 3
]

# Get the spectrum to fit
Rexp = get_reflectance(ftype, beam.λ, incident, emergent)

seed = [3300, 0.85]
options = Optim.Options(
    g_abstol=1e-8, g_reltol=1e-8, iterations=10^5, store_trace=true, show_trace=true,
)

solOptim = FitTMMOptics(
    ftype, [seed], beam, Rexp, layers;
    alg=SAMIN(), options=options, lb=[0.5.*seed], ub=[1.5.*seed],
)

solOptim2 = FitTMMOptics(
    ftype, [seed], beam, Rexp, layers;
    alg=NelderMead(), options=options,
)
