# Science...

# Load modules
using Plots, LaTeXStrings
pyplot()
using ThinFilmsTools

# Substrate index of refraction taken from www.refractiveindex.info
function nacl(x::Float64)
    x /= 1000.0
    n = sqrt(complex(
            1.0 + 0.00055+0.19800/(1.0 - (0.050/x)^2) +
            0.48398/(1.0 - (0.100/x)^2) +
            0.38696/(1.0 - (0.128/x)^2) +
            0.25998/(1.0 - (0.158/x)^2) +
            0.08796/(1.0 - (40.50/x)^2) +
            3.17064/(1.0 - (60.98/x)^2) +
            0.30038/(1.0 - (120.34/x)^2))
    )
    return n
end

function plotSpec(sol, idx)
    plt = plot(Spectrum1D(),
            sol.beam.λ, [sol.Spectra.Rp[:,1,:], sol.Spectra.Rs[:,idx,:]],
            label=["p-wave" "s-wave"],
            yaxis=("Reflectance", (0., 1.)),
            xlims=(sol.beam.λ[1], sol.beam.λ[end]),
            title=("θ = $(sol.beam.θ[idx]) [°]"),
            line=([:solid :dashdot]));
    return plt
end

function main()
    # Define beam
    λ = LinRange(4900,15000,5000) # wavelength range [nm]
    θ = LinRange(0,80,80) # angle of incidence [degrees]
    beam = PlaneWave(λ,θ)
    # Define layers
    l0 = LayerTMMO(RIdb.air(beam.λ))
    l1 = LayerTMMO(RIdb.dummy(beam.λ,4.6,0.0); d=800.)
    l2 = LayerTMMO(RIdb.dummy(beam.λ,1.6,0.0); d=1650.)
    l3 = LayerTMMO(nacl.(beam.λ)) # notice the . in nacl to use the array input
    layers = [l0, l1, l2, l1, l2, l1, l2, l1, l2, l1, l3]
    # Reference wavelength
    λ0 = 12000.0
    # call main script
    sol = tmm_optics(beam, layers; λ0=λ0, pbgflag=true)
    return sol
end

sol = main()

# Plot the R, T and A spectra, at 0, 45 and 80 degress
t45 = Utils.find_closest(sol.beam.θ,45.0)
t45 = Utils.find_closest(sol.beam.θ,80.0)
plot(plotSpec(sol, 1), plotSpec(sol, t45), plotSpec(sol, t80), layout=(3,1))
gui()

# Plot EMF
p1 = plot(Spectrum2D(),
     sol.beam.λ, sol.beam.θ, sol.Spectra.Rs, title=("s-wave Reflectance"),
);
p2 = plot(Spectrum2D(),
     sol.beam.λ, sol.beam.θ, sol.Spectra.Rp, title=("p-wave Reflectance"),
)
plot(p1, p2, layout=(2,1))
gui()

# Plot the refractive index profile
plot(RIprofile(), sol)
gui()

# plot the photonic dispersion with custom function
plot(PBGDispersion2Dalt(), sol.Bloch)
gui()
