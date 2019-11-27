# https://www.rp-photonics.com/brewsters_angle.html

# Load modules
using Plots, LaTeXStrings
pyplot()
using ThinFilmsTools

function main()
    # Define beam
    λ = [1064.]
    θ = LinRange(0.01,90,901)
    beam = PlaneWave(λ,θ)
    # Define layers: notice that lambda is outside the range
    layers = [
        LayerTMMO(RIdb.dummy(beam.λ,1.0,0.0)),
        LayerTMMO(RIdb.dummy(beam.λ,1.5,0.0), d=10.),
        LayerTMMO(RIdb.dummy(beam.λ,1.5,0.0)),
    ]
    # call main script
    sol = tmm_optics(beam, layers)
    return sol
end

sol = main()

# plot the R, T and A spectra
plot(SpectrumAngle1D(),
    sol.beam.θ,
    [sol.Spectra.Rp[1,:], sol.Spectra.Rs[1,:], sol.Spectra.Rs[1,:]./sol.Spectra.Rp[1,:]./1000.],
    label=["p-wave" "s-wave" "(Rs/Rp)/1000"],
    line=([:solid :dash :dashdot]),
    xlims=(sol.beam.θ[1], sol.Beam.θ[end]),
    yaxis=("Reflectance", (0.,0.2)),
)
gui()

# plot the refractive index profile
plot(RIprofile(), sol)
gui()
