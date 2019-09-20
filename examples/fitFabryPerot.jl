using Optim
pyplot(reuse=false, size=(640, 480), grid=false)
closeall()
using Plots
using ThinFilmsTools

function fpReflectance(beam, incident, emergent)
    l0 = LayerTMMO1DIso(RIdb.air(beam.λ))
    l1 = LayerTMMO1DIso(RI.looyengaspheresbin(0.8, RIdb.air(beam.λ), RIdb.silicon(beam.λ)), type=:OT, d=1/4.)
    l2 = LayerTMMO1DIso(RI.looyengaspheresbin(0.5, RIdb.air(beam.λ), RIdb.silicon(beam.λ)), type=:OT, d=1/4.)
    l3 = LayerTMMO1DIso(RIdb.glass(beam.λ))
    l = vec([l0 repeat([l1 l2], 1, 4) repeat([l1], 1, 2) repeat([l2 l1], 1, 4) l3])
    sol = TMMO1DIsotropic(beam, l)
    return Utils.averagePolarisation(beam.p, sol.Spectra.Rp, sol.Spectra.Rs)
end

# Wavelength range [nm]
λ = 400:1000
# Angle of incidence [degrees]
θ = [5.]
# Polarisation (p->1.0, s->0.0, weighted->between 0.0 and 1.0)
pol = 0.5
beam = PlaneWave(λ, θ; p=pol)

# Components of EMA
n1 = RIdb.air(beam.λ)
n2 = RIdb.silicon(beam.λ)

# Refractive indices of incident and emergent media
incident = deepcopy(n1)
emergent = RIdb.glass(beam.λ)

# Define the RI model to use
layers = [ LayerTMMO1DIso(incident),
           ModelFit(:looyengaspheresbin; N=(n1, n2)),
           ModelFit(:looyengaspheresbin; N=(n1, n2)),
           LayerTMMO1DIso(emergent) ]

# Create transmittance spectrum to fit
Rexp = fpReflectance(beam, incident, emergent)
plot(TMMOPlotSpectra1D(), beam.λ, Rexp)
gui()

options = Optim.Options(g_abstol=1e-8, g_reltol=1e-8, iterations=10^6, show_trace=true, store_trace=true);

seed2 = [[120, 0.75], [80.0, 0.54]]
solOptim2 = FitTMMO1DIsotropic(Reflectance(), seed2, beam, Rexp, layers; arrange=MC1d(), L=[4, 4], Ld=[2], options=options, alg=SAMIN(), lb=0.5.*seed, ub=1.5.*seed)
plot(PlotFitSpectrum(), solOptim2.beam.λ, solOptim2.spectrumExp, solOptim2.spectrumFit)
gui()

seed = [[135, 0.8], [83.0, 0.5], [0.995]]
solOptim = FitTMMO1DIsotropic(Reflectance(), beam, Rexp, layers; arrange=MC1dAlpha(), L=[4, 4], Ld=[2], options=options, alg=SAMIN(), lb=0.5.*seed, ub=1.5.*seed)
plot(PlotFitSpectrum(), solOptim.beam.λ, solOptim.spectrumExp, solOptim.spectrumFit)
gui()

