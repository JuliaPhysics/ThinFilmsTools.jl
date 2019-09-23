using Optim
using Plots
pyplot(reuse=false, size=(640, 480), grid=false)
closeall()
using ThinFilmsTools


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
n3 = RIdb.glass(beam.λ)

# Define the RI model to use
layers = [ LayerTMMO1DIso(n1),
           ModelFit(:looyengaspheresbin; N=(n1, n2)),
           ModelFit(:looyengaspheresbin; N=(n1, n2)),
           LayerTMMO1DIso(n3) ]

# Absolute reflectance spectrum to fit stored in Utils
Rexp = SpectraDB.FPSpectrum(beam.λ)
plot(TMMOPlotSpectra1D(), beam.λ, Rexp)
gui()

options = Optim.Options(g_abstol=1e-8, g_reltol=1e-8, iterations=10^6, show_trace=true, store_trace=true);

seed = [[119.0, 0.8], [76.0, 0.5]]
solOptim = FitTMMO1DIsotropic(Reflectance(), seed, beam, Rexp, layers; arrange=MC1d(), L=[4, 4], Ld=[2], options=options, alg=SAMIN(), lb=0.5.*seed, ub=1.5.*seed)

plot(PlotFitSpectrum(), solOptim.beam.λ, solOptim.spectrumExp, solOptim.spectrumFit)
gui()

seed2 = [[120, 0.8], [83.0, 0.5], [0.995]]
solOptim2 = FitTMMO1DIsotropic(Reflectance(), seed2, beam, Rexp, layers; arrange=MC1dAlpha(), L=[4, 4], Ld=[2], options=options, alg=SAMIN(), lb=0.5.*seed2, ub=1.5.*seed2)
plot(PlotFitSpectrum(), solOptim2.beam.λ, solOptim2.spectrumExp, solOptim2.spectrumFit)
gui()


