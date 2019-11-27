# Sensors and Actuators B 149 (2010) 189-193

# Load modules
using Plots, LaTeXStrings
pyplot()
using ThinFilmsTools

function main()
    # Define beam
    λ = LinRange(400,1000,500) # wavelength range [nm]
    θ = [0.] # angle of incidence [degrees]
    beam = PlaneWave(λ,θ)
    # Define layers with d=1/4 by default
    p1 = 0.86
    p2 = 0.54
    n = [RIdb.air(beam.λ) RIdb.silicon(beam.λ)]
    l0 = LayerTMMO(RIdb.air(beam.λ)) # incident medium
    l1 = LayerTMMO(RI.looyenga([p1 1-p1],n); type=:OT)
    l2 = LayerTMMO(RI.looyenga([p2 1-p2],n); type=:OT)
    l3 = LayerTMMO(RIdb.glass(beam.λ./1e3)) # emergent medium
    layers = vec([l0 l1 l2 l1 l2 l1 l2 l1 l2 repeat([l1], 1, 10) l2 l1 l2 l1 l2 l1 l2 l1 l3])
    # Reference wavelength
    λ0 = 730.0
    # call main script
    sol = tmm_optics(beam, layers; λ0=λ0, emfflag=true, h=10)
    return sol
end

sol = main()

# plot the R, T and A spectra
plot(Spectrum1D(),
    sol.beam.λ, [sol.Spectra.Rp, sol.Spectra.Tp, 1.0.-(sol.Spectra.Rp.+sol.Spectra.Tp)],
    label=["Reflectance" "Transmittance" "Absorbance"],
    line=([:solid :dash :dashdot]),
    xlims=(sol.beam.λ[1], sol.Beam.λ[end]),
    yaxis=("Transmittance", (0.,1.0)),
)
gui()

plot(EMF2D(), sol.beam.λ, sol.Misc.ℓ, sol.Field.emfp[:,1,:])
gui()

# plot the refractive index profile
plot(RIprofile(), sol)
gui()
