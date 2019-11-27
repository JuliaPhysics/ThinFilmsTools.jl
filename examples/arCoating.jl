# Load modules
using Plots, LaTeXStrings
pyplot()
using ThinFilmsTools

function main()
    # Define beam
    λ = LinRange(400,1000,600) # wavelength range [nm]
    θ = [0.] # angle of incidence [degrees]
    beam = PlaneWave(λ,θ)
    # Define layers with their parameters
    n = [RIdb.air(beam.λ) RIdb.silicon(beam.λ)]
    p1 = 0.89
    p2 = 0.70
    p3 = 0.41
    layers = [
        LayerTMMO(n[:,1]), # incident medium
        LayerTMMO(RI.looyenga([p1 1-p1],n); d=77.),
        LayerTMMO(RI.looyenga([p2 1-p2],n); d=56.),
        LayerTMMO(RI.looyenga([p3 1-p3],n), d=39.),
        LayerTMMO(n[:,2]),  # emergent (substrate) medium
    ]
    # solve
    return tmm_optics(beam, layers)
end

sol = main()

# plot the R, T and A spectra
plot(Spectrum1D(),
     sol.beam.λ,
     [sol.Spectra.Rp, sol.Spectra.Tp, 1.0 .- (sol.Spectra.Rp .+ sol.Spectra.Tp)],
     label=["Reflectance" "Transmittance" "Absorbance"],
     line=([:solid :dash :dashdot]),
)
gui()

# plot the refractive index profile
plot(RIprofile(), sol)
gui()
