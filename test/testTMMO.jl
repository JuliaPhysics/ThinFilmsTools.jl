function main()
    # Define beam
    λ = LinRange(200,5000,1000) # wavelength range [nm]
    θ = LinRange(0,90,90)
    beam = PlaneWave(λ,θ)
    # Define layers
    l0 = LayerTMMO(RIdb.air(beam.λ))
    l1 = LayerTMMO(RIdb.dummy(beam.λ,1.45,0.0); type=:OT, d=1/4.)
    l2 = LayerTMMO(RIdb.dummy(beam.λ,3.45,0.0); type=:OT, d=1/4.)
    layers = [l0, l1, l2, l1, l2, l1, l2, l1, l2, l0]
    # Reference wavelenth
    λ0 = 700.0
    # call main script
    sol = TMMOptics(beam, layers; λ0=λ0, emfflag=true, h=10, pbgflag=true)
    return sol
end


sol = main()

