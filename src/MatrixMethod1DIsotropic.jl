module MatrixMethod1DIsotropic

using Statistics
using LinearAlgebra
using ..CommonStructures: getBeamParameters, PlaneWave, LayerTMMO1DIso

export TMMO1DIsotropic

"""Wrap the output into different types."""
abstract type Output end
struct Spectra{T1, T2} <: Output where {T1<:Float64, T2<:ComplexF64}
    Rp::Array{T1}; Rs::Array{T1}; Tp::Array{T1}; Ts::Array{T1};
    ρp::Array{T2}; ρs::Array{T2}; τp::Array{T2}; τs::Array{T2};
end
struct Field{T1} <: Output where {T1<:Float64}
    emfp::Array{T1}; emfs::Array{T1};
end
struct Bloch{T1, T2} <: Output where {T1<:ComplexF64, T2<:Float64}
    κp::Array{T1}; κs::Array{T1}; ω::Array{T2}; Λ::T2; ωl::T2; ωh::T2; qz::Array
end
struct Misc{T1, T2} <: Output where {T1<:Float64, T2<:ComplexF64}
    d::Array{T1}; ℓ::Array{T1}; nλ0::Array{T1}; layers_n::Array{T2}; λ0::T1
end
struct AdmPhase{T1} <: Output where {T1<:ComplexF64}
    ηp::Array{T1}; ηs::Array{T1}; δ::Array{T1};
end
struct TMMO1DIsotropic{T1, T2, T3, T4, T5, T6, T7} <: Output where {T1<:Spectra, T2<:Field, T3<:Bloch, T4<:Misc, T5<:AdmPhase, T6<:PlaneWave, T7<:LayerTMMO1DIso}
    Spectra::T1; Field::T2; Bloch::T3; Misc::T4; AdmPhase::T5; Beam::T6; Layers::Array{T7};
end

"""

    Data conditioning and main calculations.

"""
function TMMO1DIsotropic(beam::T1, layers::Array{T2,N2}; emfflag::T3=false, h::T4=10, pbgflag::T3=false, λ0::T5=-eps()) where {T1<:PlaneWave, T2<:LayerTMMO1DIso, N2, T3<:Bool, T4<:Int64, T5<:Float64}
    # Check if λ0 ∈ λ
    λ0 = checkλ0(λ0, beam)
    # Find beam.λ closest to λ0
    idxλ0 = findmin(abs.(beam.λ .- λ0))[2][1]
    # Build the sequence of n and d depending on the input
    nLen = length(vec(layers))
    d = Array{Float64,1}(undef, nLen)
    nλ0 = Array{Float64,1}(undef, nLen)
    layers_n = Array{ComplexF64,2}(undef, length(beam.λ), nLen)
    buildArrays!(d, layers_n, nλ0, idxλ0, λ0, vec(layers), beam)
    # Provide the multilayer depth considering the h division
    _ℓ::Array{Float64,2} = (d[2:end-1]/h)*ones.(1, h) # outer product
    ℓ = cumsum([0; _ℓ[:]], dims=1)[1:end-1] # remove last from cumsum
    # Get the parameters for beam
    beam_, _ = getBeamParameters(beam)
    # Call transfer matrix method
    tmmout = emfflag ? transferMatrixEMF(layers_n, d, beam_, h) : TransferMatrix(layers_n, d, beam_)
    # Photonic band gap for crystals without defects
    if pbgflag & (nLen > 3)
        ω = Array{Float64,1}(undef, length(beam_.λ))
        κp = Array{ComplexF64,2}(undef, (length(beam_.λ), length(beam_.θ)))
        κs = similar(κp)
        qz::Array{Float64,1} = sin.(beam_.θ).*π./2.0 # parallel wavevector qz
        Λ::Float64 = sum(d[2:3])
        d1, d2 = d[2], d[3]
        n0, n1, n2 = nλ0[1], nλ0[2], nλ0[3]
        ωh::Float64 = Λ/π/(d1*n1 + d2*n2)*acos(-abs(n1 - n2)/(n1 + n2))
        ωl::Float64 = Λ/π/(d2*sqrt(n2^2 - n0^2) + d1*sqrt(n1^2 - n0^2))*acos(abs((n1^2*sqrt(n2^2 - n0^2) - n2^2*sqrt(n1^2 - n0^2))/(n1^2*sqrt(n2^2 - n0^2) + n2^2*sqrt(n1^2 - n0^2))))
        photonicDispersion!(κp, κs, ω, Λ, beam_, layers_n[idxλ0, 2:3], d[2:3])
    else
        κp = []; κs = []; ω = []; Λ = 0.0; ωl = 0.0; ωh = 0.0; qz = []
    end
    # Return results
    TMMO1DIsotropic(tmmout[1], tmmout[2], Bloch(κp, κs, ω, Λ, ωl, ωh, qz), Misc(d, ℓ, nλ0, layers_n, λ0), tmmout[3], beam, layers)
end

"""

    Check input λ0 againts beam.λ.

"""
function checkλ0(λ0, beam)
    return isequal(λ0, -eps()) ? mean(beam.λ) : λ0
end

"""

    Build the sequence of indices of refraction and thicknesses depending on the input. It follows some logic depending on whether nλ0 was input.

"""
function buildArrays!(d::Array{T0,1}, layers_n::Array{T1,2}, nλ0::Array{T0,1}, idxλ0::T4, λ0::T0, layers::Array{T2,1}, beam::T5) where {T0<:Float64, T1<:ComplexF64, T2<:LayerTMMO1DIso, T3<:Bool, T4<:Int64, T5<:PlaneWave}
    @inbounds for s in eachindex(layers)
        # Refractive index
        @views layers_n[:, s] = layers[s].n
        # nλ0 depending on the input
        if !isequal(layers[s].nλ0[1], -eps()) # nλ0 specified
            nλ0[s] = real(layers[s].nλ0[1])
        else # nλ0 not specified
            nλ0[s] = real(layers_n[idxλ0, s])
        end
        # Build thickness depending on the input
        d[s] = layers[s].d
        if layers[s].type == :OT
            d[s] *= λ0/nλ0[s]
        end
    end
    return nothing
end

"""

    Computes the reflection and transmission coefficients and spectra with the transfer matrix method.

"""
function TransferMatrix(layers_n::Array{T1,2}, d::Array{T2,1}, beam::T3; λLen::T4=length(beam.λ), θLen::T4=length(beam.θ)) where {T1<:ComplexF64, T2<:Float64, T3<:PlaneWave, T4<:Int64}
    # λ, λLen, θ, θLen = getBeamParameters(beam)
    nLen::Int64 = size(layers_n, 2)
    τs = Array{ComplexF64,2}(undef, (λLen, θLen))
    τp = similar(τs); ρs = similar(τs); ρp = similar(τs)
    δ = Array{ComplexF64,3}(undef, (λLen, θLen, nLen))
    ηs = similar(δ); ηp = similar(ηs)
    # Calculation of complex coefficients of reflection, transmission and emf
    fresnelCoefficients!(ηs, ηp, δ, ρs, ρp, τs, τp, layers_n, d, beam, nLen)
    return (Spectra(abs2.(ρp), abs2.(ρs), real(@view(ηp[:, :, 1]).*@view(ηp[:, :, end])).*abs2.(τp), real(@view(ηs[:, :, 1]).*@view(ηs[:, :, end])).*abs2.(τs), ρp, ρs, τp, τs), Field([], []), AdmPhase(ηp, ηs, δ))
end

"""

    Compute the Fresnel coefficients.

"""
function fresnelCoefficients!(ηs::Array{T1,3}, ηp::Array{T1,3}, δ::Array{T1,3}, ρs::Array{T1,2}, ρp::Array{T1,2}, τs::Array{T1,2}, τp::Array{T1,2}, layers_n::Array{T1,2}, d::Array{T2,1}, beam::T3, nLen::T4) where {T1<:ComplexF64, T2<:Float64, T3<:PlaneWave, T4<:Int64}
    cosϕ = Array{ComplexF64,1}(undef, nLen)
    I2d = Array{ComplexF64,2}(I, 2, 2) # Identity 2x2 matrix
    @inbounds for a in eachindex(beam.θ), l in eachindex(beam.λ)
        # Angle of incidence inside each layer according to the cosine Snell law, to avoid cutoff of total internal reflection with complex angles
        cosϕ[1] = cos(beam.θ[a])
        for c = 2:nLen
            cosϕ[c] = cosϑ(layers_n[l, c-1], layers_n[l, c], cosϕ[c-1])
        end
        # Admittance of the first medium for both polarizations
        @. ηp[l, a, :] = ζp(@view(layers_n[l, :]), cosϕ)
        @. ηs[l, a, :] = ζs(@view(layers_n[l, :]), cosϕ)
        # Phase shift of the thickness
        @. δ[l, a, :] = 2.0*π*@view(layers_n[l, :])*d*cosϕ/beam.λ[l]
        # Total transfer matrix
        Ψp::Array{ComplexF64,2} = reduce(*, [[I2d]; Φ.(@view(δ[l, a, 2:end-1]), @view(ηp[l, a, 2:end-1]))])
        Ψs::Array{ComplexF64,2} = reduce(*, [[I2d]; Φ.(@view(δ[l, a, 2:end-1]), @view(ηs[l, a, 2:end-1]))])
        # Compute the Fresnell coefficients
        ρs[l, a] = ρ(ηs[l, a, 1], ηs[l, a, end], Ψs)
        ρp[l, a] = ρ(ηp[l, a, 1], ηp[l, a, end], Ψp)
        τs[l, a] = τ(ηs[l, a, 1], ηs[l, a, end], Ψs)
        τp[l, a] = τ(ηp[l, a, 1], ηp[l, a, end], Ψp)
    end
    return nothing
end

"""

    Computes the reflection and transmission coefficients, spectra and the electromagnetic field with the transfer matrix method.

"""
function transferMatrixEMF(layers_n::Array{T1,2}, d::Array{T2,1}, beam::T3, h::T4; λLen::T4=length(beam.λ), θLen::T4=length(beam.θ)) where {T1<:ComplexF64, T2<:Float64, T3<:PlaneWave, T4<:Int64}
    # λ, λLen, θ, θLen = getBeamParameters(beam)
    nLen::Int64 = size(layers_n, 2)
    τs = Array{ComplexF64,2}(undef, (λLen, θLen))
    τp = similar(τs); ρs = similar(τs); ρp = similar(τs)
    δ = Array{ComplexF64,3}(undef, (λLen, θLen, nLen))
    ηs = similar(δ); ηp = similar(ηs)
    emfs = Array{Float64,3}(undef, (λLen, θLen, (length(d) - 2)*h)); emfp = similar(emfs)
    # Calculation of complex coefficients of reflection, transmission and emf
    fresnelCoefficientsEMF!(ηs, ηp, δ, ρs, ρp, τs, τp, emfs, emfp, layers_n, d, beam, nLen, h)
    return (Spectra(abs2.(ρp), abs2.(ρs), real(@view(ηp[:, :, 1]).*@view(ηp[:, :, end])).*abs2.(τp), real(@view(ηs[:, :, 1]).*@view(ηs[:, :, end])).*abs2.(τs), ρp, ρs, τp, τs), Field(emfp, emfs), AdmPhase(ηp, ηs, δ))
end

"""

    Compute the Fresnel coefficients and the electromagnetic field.

"""
function fresnelCoefficientsEMF!(ηs::Array{T1,3}, ηp::Array{T1,3}, δ::Array{T1,3}, ρs::Array{T1,2}, ρp::Array{T1,2}, τs::Array{T1,2}, τp::Array{T1,2}, emfs::Array{T2,3}, emfp::Array{T2,3}, layers_n::Array{T1,2}, d::Array{T2,1}, beam::T3, nLen::T4, h::T4) where {T1<:ComplexF64, T2<:Float64, T3<:PlaneWave, T4<:Int64}
    cosϕ = Array{ComplexF64,1}(undef, nLen)
    I2d = Array{ComplexF64,2}(I, 2, 2) # Identity 2x2 matrix
    @inbounds for a in eachindex(beam.θ), l in eachindex(beam.λ)
        # Angle of incidence inside each layer according to the cosine Snell law, to avoid cutoff of total internal reflection with complex angles
        cosϕ[1] = cos(beam.θ[a])
        for c = 2:nLen
            cosϕ[c] = cosϑ(layers_n[l, c-1], layers_n[l, c], cosϕ[c-1])
        end
        # Admittance of the first medium for both polarizations
        @. ηp[l, a, :] = ζp(@view(layers_n[l, :]), cosϕ)
        @. ηs[l, a, :] = ζs(@view(layers_n[l, :]), cosϕ)
        # Phase shift of the thickness
        @. δ[l, a, :] = 2.0*π*@view(layers_n[l, :])*d*cosϕ/beam.λ[l]
        # Total transfer matrix
        Ψp::Array{ComplexF64,2} = reduce(*, [[I2d]; Φ.(@view(δ[l, a, 2:end-1]), @view(ηp[l, a, 2:end-1]))])
        Ψs::Array{ComplexF64,2} = reduce(*, [[I2d]; Φ.(@view(δ[l, a, 2:end-1]), @view(ηs[l, a, 2:end-1]))])
        # Compute the Fresnell coefficients
        ρs[l, a] = ρ(ηs[l, a, 1], ηs[l, a, end], Ψs)
        ρp[l, a] = ρ(ηp[l, a, 1], ηp[l, a, end], Ψp)
        τs[l, a] = τ(ηs[l, a, 1], ηs[l, a, end], Ψs)
        τp[l, a] = τ(ηp[l, a, 1], ηp[l, a, end], Ψp)
        # Compute the emfield
        emfs[l, a, :] .= emfield(@view(layers_n[l, :]), d, @view(δ[l, a, :]), @view(ηs[l, a, :]), Ψs, nLen, h)
        emfp[l, a, :] .= emfield(@view(layers_n[l, :]), d, @view(δ[l, a, :]), @view(ηp[l, a, :]), Ψp, nLen, h)
    end
    return nothing
end

"""

    Calculates the optical transfer matrix of a layer. (2x2 optical tranfer matrix.)
        δ: phase shift of the layer
        η: admittance of the layer

"""
Φ(δ::T1, η::T1) where {T1<:ComplexF64} = [cos(δ) (-im*sin(δ)/η); (-im*sin(δ)*η) cos(δ)]

"""

    Computes the reflection and transmission coefficients given the admittance and transfer matrix of the whole structure per wavelenth and angle of incidence.

"""
# Reflection coefficient
ρ(η0::T1, ηm::T1, Ψ::Array{T1,2}) where {T1<:ComplexF64} = (η0*Ψ[1,1] - Ψ[2,1] + η0*ηm*Ψ[1,2] - ηm*Ψ[2,2]) / (η0*Ψ[1,1] + Ψ[2,1] + η0*ηm*Ψ[1,2] + ηm*Ψ[2,2])
# Transmission coefficient
τ(η0::T1, ηm::T1, Ψ::Array{T1,2}) where {T1<:ComplexF64} = 2 / (η0*Ψ[1,1] + Ψ[2,1] + η0*ηm*Ψ[1,2] + ηm*Ψ[2,2])

"""

    Computes the inverse total transfer matrix for the whole structure at each wavelenth and angle of incidence and return the field.

"""
function emfield(N::SubArray{T1,1}, d::Array{T2,1}, δ::SubArray{T1,1}, η::SubArray{T1,1}, Ψ::Array{T1,2}, nLen::T3, h::T3) where {T1<:ComplexF64, T2<:Float64, T3<:Int64}
    m0 = Array{ComplexF64,2}(undef, 2, 2)
    m1 = Array{ComplexF64,2}(I, 2, 2) # Identity 2x2 matrix
    g11 = Array{ComplexF64,1}(undef, (nLen - 2)*h); g12 = similar(g11)
    # Divide the phase shift by h but keep η as is for each layer
    mδ::Array{ComplexF64} = δ/h
    _m1::Array = Ξ.(mδ, η)
    fieldIntensity!(g11, g12, m0, _m1, m1, Ψ, nLen, h)
    return FI(g11, g12, η[1], η[end], Ψ)
end

"""

	Compute the electric field distribution.

"""
function fieldIntensity!(g11::Array{T1,1}, g12::Array{T1,1}, m0::Array{T1,2}, _m1::Array{Array{T1,2},1}, m1::Array{T1,2}, Ψ::Array{T1,2}, nLen::T2, h::T2) where {T1<:ComplexF64, T2<:Int64}
    _h = collect(1:h)
    @inbounds for c =2:nLen-1, j in eachindex(_h)
        k::Int64 = h*(c - 2) + j
        m1 = _m1[c]*m1
        m0 = m1*Ψ
        g11[k] = m0[1, 1]
        g12[k] = m0[1, 2]
    end
    return nothing
end
FI(g11::Array{T1,1}, g12::Array{T1,1}, η1::T1, ηf::T1, Ψ::Array{T1,2}) where {T1<:ComplexF64} = abs2.((g11 .+ ηf.*g12)./(0.25*(η1*Ψ[1,1] + Ψ[2,1] + η1*ηf*Ψ[1,2] + ηf*Ψ[2,2])))

"""

	Calculates the inverse of optical transfer matrix of a layer. φ:  phase shift of the layer, y: admittance of the layer, Τ: 2x2 optical tranfer matrix.

"""
Ξ(φ::T1, η::T1) where {T1<:ComplexF64} = [cos(φ) (im*sin(φ)/η); (im*sin(φ)*η) cos(φ)]

"""

    Computes the photonic dispersion of binary structures (crystals only) alternating two different dielectric layers.

"""
function photonicDispersion!(κp::Array{T1,2}, κs::Array{T1,2}, ω::Array{T2,1}, Λ::T2, beam::T3, n::Array{T1,1}, d::Array{T2,1}) where {T1<:ComplexF64, T2<:Float64, T3<:PlaneWave}
    @. ω = 2*π/beam.λ # Angular frequency
    # Angle of incidence of the second layer with Snell's law of cosine
    cosθ1::Vector{ComplexF64} = cos.(beam.θ)
    cosθ2::Vector{ComplexF64} = cosϑ.(n[1], n[2], cosθ1)
    # Prefactor for Bloch wavevector
    factor_s = admFactor.(ζs.(n[1], cosθ1), ζs.(n[2], cosθ2))
    factor_p = admFactor.(ζp.(n[1], cosθ1), ζp.(n[2], cosθ2))
    # Bloch wavevectors
    @inbounds for a in eachindex(cosθ1), b in eachindex(ω)
        κp[b, a] = acos(cosκ(d[1]*ω[b]*n[1]*cosθ1[a], d[2]*ω[b]*n[2]*cosθ2[a], factor_p[a]))/Λ
        κs[b, a] = acos(cosκ(d[1]*ω[b]*n[1]*cosθ1[a], d[2]*ω[b]*n[2]*cosθ2[a], factor_s[a]))/Λ
    end
    return nothing
end

"""

	Prefactor for bloch wavevector.

"""
admFactor(η1::T1, η2::T1) where {T1<:ComplexF64} = 0.5*(η1^2 + η2^2)/η1/η2

"""

	Bloch wavevector.

"""
cosκ(a1::T1, a2::T1, f::T1) where {T1<:ComplexF64} = cos(a1)*cos(a2) - f*sin(a1)*sin(a2)

"""

	Snell's law in cosine form. Returns the cosine already.

"""
cosϑ(n1::T1, n2::T1, cosθ::T2) where {T1<:ComplexF64, T2<:Number} = sqrt(1.0 - (n1/n2)^2*(1.0 - cosθ^2))

"""

	Admittance of the medium for p and s polarizations.

"""
ζp(n::T1, cosθ::T2) where {T1<:ComplexF64, T2<:Number} = n/cosθ
ζs(n::T1, cosθ::T2) where {T1<:ComplexF64, T2<:Number} = n*cosθ

end # module
