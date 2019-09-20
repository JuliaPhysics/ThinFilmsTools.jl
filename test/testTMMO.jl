# Define beam
λi = 200 # intial wavelength [nm]
λf = 5000 # final wavelength [nm]
λ = LinRange(λi, λf, 1000) # wavelength range [nm]
λ0 = 700. # reference wavelength
θi = 0
θf = 90
θ = LinRange(0, 90, θf-θi)
beam = PlaneWave(λ, θ)

# Define layers
l0 = LayerTMMO1DIso(RIdb.air(beam.λ))
l1 = LayerTMMO1DIso(RIdb.dummy(beam.λ, 1.45, 0.); type=:OT, d=1/4.)
l2 = LayerTMMO1DIso(RIdb.dummy(beam.λ, 3.45, 0.); type=:OT, d=1/4.)
layers = [l0, l1, l2, l1, l2, l1, l2, l1, l2, l0]

# Reference wavelenth
λ0 = 700.0

# call main script
sol = TMMO1DIsotropic(beam, layers; λ0=λ0, emfflag=true, h=10, pbgflag=true)

