# https://www.horiba.com/en_en/spectroscopic-ellipsometry/
# Typically, ellipsometers do not measure ψ and Δ directly. Instead, they measure functions of ψ and Δ, for instance: Is, Ic, and Ic', which are functions of ψ and Δ according to Is = sin2ψ.sin Δ, Ic = sin 2ψ cos Δ, and Ic' = cos 2ψ. When combined, Is and Ic provide an accurate measurement of Δ over the full range from 0° to 360° and Is and Ic' provide an accurate measurement of ψ over the full range from 0° to 90°.
# The fitting procedure works on Ψ and Δ directly, since different ellipsometers may use different setup. Therefore the user must convert the Is, Ic, or whatever parameters into Ψ and Δ, with Δ into radians. Important: the program works for Ψ ∈ [0, π/2] and Δ ∈ [-π, π]. Sometimes equipment gives Δ from 0 to 2π, if that is the case, you need to subtract π from your data.

using Plots
pyplot(reuse=false, size=(640, 480), grid=false)
closeall()
using Optim
using ThinFilmsTools

##
function glassDeltaPsi(beam, incident, emergent)
    # Create a measured spectrum of reflection
    layers2 = [ LayerTMMO1DIso(incident),
                LayerTMMO1DIso(RIdb.glass(beam.λ./1e3); d=400.),
                LayerTMMO1DIso(emergent) ]
    sol = TMMO1DIsotropic(beam, layers2)
    ratio = sol.Spectra.ρp./sol.Spectra.ρs
    return atan.(abs.(ratio)), angle.(ratio)
end
##

# Wavelength range [nm]
λ = 400:1000
# Angle of incidence [degrees]
θ = [60.]
# Polarisation (1.0 = p, 0.0 = s, between 0.0 and 1.0 = average)
pol = 0.5
beam = PlaneWave(λ, θ; p=pol)

# Refractive indices of incident (0) and substrate (2)
incident = RIdb.air(beam.λ)
emergent = RIdb.silicon(beam.λ)

# Define the RI model to use
layers = [ LayerTMMO1DIso(incident),
           ModelFit(:cauchyurbach),
           LayerTMMO1DIso(emergent) ]

# Create transmittance spectrum to fit
Ψ, Δ = glassDeltaPsi(beam, incident, emergent)
plot(TMMOPlotSpectra1D(), beam.λ, Ψ, tickfont=font(12), legendfont=font(10))
plot!(TMMOPlotSpectra1D(), beam.λ, Δ, tickfont=font(12), legendfont=font(10))
gui()

options = Optim.Options(g_abstol=1e-8, g_reltol=1e-8, iterations=10^5, show_trace=true, store_trace=true);

seed = [vcat(420.0, [1.4, 0.23, 1.0, 0.0, 1.0, 1.0])]
solOptim = FitTMMO1DIsotropic(Ellipsometry(), seed, beam, hcat(Ψ, Δ), layers; options=options, alg=SAMIN(), lb=0.5.*seed, ub=1.5.*seed)

plot(PlotFitSpectrumEllip(), beam.λ, solOptim.spectrumExp, solOptim.spectrumFit, xaxis=("Wavelength [nm]"), tickfont=font(12), legendfont=font(10))
gui()

plot(PlotFitSpectrumEllip(), beam.λ, rad2deg.(solOptim.spectrumExp), rad2deg.(solOptim.spectrumFit), xaxis=("Wavelength [nm]"), tickfont=font(12), legendfont=font(10))
gui()

