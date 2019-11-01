# Load modules
using Plots, LaTeXStrings
pyplot()
using ThinFilmsTools

function main()
    # Define beam
    λ = LinRange(200,1000,600) # wavelength range [nm]
    θ = [0.] # angle of incidence [degrees]
    beam = PlaneWave(λ, θ)
    # Define layers
    l0 = LayerTMMO(RIdb.air(beam.λ))
    l1 = LayerTMMO(RIdb.dummy(beam.λ,1.45,0.0); type=:OT, d=1/4.)
    l2 = LayerTMMO(RIdb.dummy(beam.λ,3.45,0.0); type=:OT, d=1/4.)
    layers = [l0, l1, l2, l1, l2, l1, l2, l1, l2, l0]
    # Reference wavelenth [nm]
    λ0 = 700.
    # Call main script
    sol = TMMOptics(beam, layers; λ0=λ0, emfflag=true, h=10, pbgflag=true)
    return sol
end

sol = main()

# Plot spectra
plot(Spectrum1D(),
    sol.Beam.λ, [sol.Spectra.Rp, sol.Spectra.Tp, 1.0.-(sol.Spectra.Rp.+sol.Spectra.Tp)], 
    label=["Reflectance" "Transmittance" "Absorbance"],
    line=([:solid :dash :dashdot]),
    ylims=(0.0,1.0),
    xlims=(sol.Beam.λ[1], sol.Beam.λ[end]),
);
gui()

# Plot the EMF pattern for normal incidence
plot(EMF2D(),
    sol.Beam.λ, sol.Misc.ℓ, log10.(sol.Field.emfp[:,1,:]),
    title=("Log of EMF intensity"),
);
gui()

# plot the refractive index profile
plot(RIprofile(), sol; plotemf=true)
gui()

# plot the photonic dispersion with custom function
plot(PBGDispersion1D(), sol.Bloch)
gui()

plot(PBGDispersion1Dalt(), sol.Bloch)
gui()

plot(PBGDispersion1Dimre(), sol.Bloch)
gui()
