function dummy(λ::AbstractArray{T,N}, x::S, y::S) where {T<:Number, N, S<:Number}
    Nc = ones.(length(λ)) .* (x + im*y)
end # EOF dummy(...)

# Define beam
λi = 400 # intial wavelength [nm]
λf = 1000 # final wavelength [nm]
λ = LinRange(λi, λf, λf-λi+1) # wavelength range [nm]
λ0 = 700. # reference wavelength
θ = [0.] # angle of incidence [degrees]
beam = PlaneWave(λ, λ0, θ)

# Define layers
l0 = LayerTMMO1DIso(type=:GT, n=dummy(beam.λ, 1., 0.), d=0.)
l1 = LayerTMMO1DIso(type=:GT, n=dummy(beam.λ, 2.5, 0.5), d=150.)
l2 = LayerTMMO1DIso(type=:GT, n=dummy(beam.λ, 1.5, 0), d=0.)
layers = [l0, l1, l2, l1, l2, l1, l2, l1, l2, l0]

# call main script
sol = TMMO1DIsotropic(beam, layers; emfflag=true, h=10, pbgflag=true)

