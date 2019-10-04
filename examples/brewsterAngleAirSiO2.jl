# https://www.rp-photonics.com/brewsters_angle.html

# Load modules
using Plots, LaTeXStrings
pyplot(reuse=false, grid=false)
closeall()
using ThinFilmsTools

function main()
    # Define beam
    λ = [1064.]
    θ = LinRange(0.01, 90, 901)
    beam = PlaneWave(λ, θ)
    # Define layers: notice that lambda is outside the range
    layers = [ LayerTMMO1DIso(RIdb.dummy(beam.λ, 1., 0.)),
               LayerTMMO1DIso(RIdb.dummy(beam.λ, 1.5, 0.), d=10.),
               LayerTMMO1DIso(RIdb.dummy(beam.λ, 1.5, 0.)) ]
    # call main script
    return TMMO1DIsotropic(beam, layers)
end

sol = main()

# plot the R, T and A spectra
plot(TMMOPlotSpectraAngle1D(),
     sol.Beam.θ, [sol.Spectra.Rp[1,:], sol.Spectra.Rs[1,:], sol.Spectra.Rs[1,:]./sol.Spectra.Rp[1,:]./1000.],
     label=["p-wave" "s-wave" "(Rs/Rp)/1000"], line=([:solid :dash :dashdot]), xlims=(sol.Beam.θ[1], sol.Beam.θ[end]), yaxis=("Reflectance", (0.,0.2)))
gui()

# plot the refractive index profile
plot(TMMOPlotNprofile(), sol)
gui()
