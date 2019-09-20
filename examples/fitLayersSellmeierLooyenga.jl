using Plots
pyplot(reuse=false, size=(640, 480), grid=false)
closeall()
using Optim
using ThinFilmsTools

##
function bilayerReflectance(beam, incident, emergent)
    # Create a measured spectrum of reflection
    layers2 = [ LayerTMMO1DIso(incident),
                LayerTMMO1DIso(RI.looyengaspheresbin(0.70,incident,emergent); d=300.),
                LayerTMMO1DIso(RIdb.glass(beam.λ./1e3); d=250.),
                LayerTMMO1DIso(emergent) ]
    sol = TMMO1DIsotropic(beam, layers2)
    return beam.p.*sol.Spectra.Rp .+ (1.0 - beam.p).*sol.Spectra.Rs
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
layers = [ LayerTMMO1DIso(incident),
           ModelFit(:looyengaspheresbin; N=(incident, emergent)),
           ModelFit(:sellmeier),
           LayerTMMO1DIso(emergent) ]

# Create absolute reflectance spectrum to fit
Rexp_norm = bilayerReflectance(beam, incident, emergent)
plot(TMMOPlotSpectra1D(), beam.λ, Rexp_norm)
gui()

options = Optim.Options(g_abstol=1e-8, g_reltol=1e-8, iterations=10^6, show_trace=true, store_trace=true);

seed = [[300, 0.7],
        vcat(250.0, [1.0, 0.23, 1.0, 6e-3, 2e-2, 1e-2])]
solOptim = FitTMMO1DIsotropic(Reflectance(), seed, beam, Rexp_norm, layers; options=options, alg=SAMIN(), lb=0.5.*seed, ub=1.5.*seed)

plot(PlotFitSpectrum(), solOptim.beam.λ, solOptim.spectrumExp, solOptim.spectrumFit)
gui()
