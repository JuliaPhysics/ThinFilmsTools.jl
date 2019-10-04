# Load modules
using Plots, LaTeXStrings
pyplot(reuse=false, grid=false, size=(640,480))
closeall()
using ThinFilmsTools

function main()
    # Define beam
    λ = LinRange(310, 1000, 500) # wavelength range [nm]
    θ = [0.] # angle of incidence [degrees]
    beam = PlaneWave(λ, θ)
    # Find wavelength closest to 405.0 nm
    aux1 = findmin(abs.(beam.λ .- 405.0))[2][1]
    # Define layers
    layers = [LayerTMMO1DIso(RIdb.air(beam.λ)),
              LayerTMMO1DIso(RIdb.sno2f(beam.λ); d=150., nλ0=RIdb.sno2f([beam.λ[aux1]])),
              LayerTMMO1DIso(RIdb.silicon(beam.λ))]
    # call main script
    return TMMO1DIsotropic(beam, layers)
end

sol = main()

# plot the R, T and A spectra
plot(TMMOPlotSpectra1D(),
     sol.Beam.λ, [sol.Spectra.Rp, sol.Spectra.Tp, 1.0.-(sol.Spectra.Rp.+sol.Spectra.Tp)],
     label=["Reflectance" "Transmittance" "Absorbance"], line=([:solid :dash :dashdot]),
     ylims=(0.0,1.0), xlims=(sol.Beam.λ[1], sol.Beam.λ[end]))
gui()

# plot the refractive index profile
plot(TMMOPlotNprofile(), sol)
gui()
