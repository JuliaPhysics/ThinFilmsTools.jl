module TransferMatrixMethod

using Statistics
using LinearAlgebra
using ..CommonStructures: _get_beam_parameters, PlaneWave, LayerTMMO

export TMMOptics, tmm_optics

# Wrap the output into different types
abstract type Output end
struct Spectra{T1, T2} <: Output where {T1<:Float64, T2<:ComplexF64}
    Rp::Array{T1}; Rs::Array{T1}
    Tp::Array{T1}; Ts::Array{T1}
    ρp::Array{T2}; ρs::Array{T2}
    τp::Array{T2}; τs::Array{T2}
end
struct Field{T1} <: Output where {T1<:Float64}
    emfp::Array{T1}
    emfs::Array{T1}
end
struct Bloch{T1, T2} <: Output where {T1<:ComplexF64, T2<:Float64}
    κp::Array{T1}
    κs::Array{T1}
    ω::Array{T2}
    Λ::T2
    ωl::T2
    ωh::T2
    qz::Array
end
struct Misc{T1, T2} <: Output where {T1<:Float64, T2<:ComplexF64}
    d::Array{T1}
    ℓ::Array{T1}
    nλ0::Array{T1}
    layers_n::Array{T2}
    λ0::T1
end
struct AdmPhase{T1} <: Output where {T1<:ComplexF64}
    ηp::Array{T1}
    ηs::Array{T1}
    δ::Array{T1}
end
struct TMMOptics{T1, T2, T3, T4, T5, T6, T7} <: Output where {T1<:Spectra, T2<:Field, T3<:Bloch, T4<:Misc, T5<:AdmPhase, T6<:PlaneWave, T7<:LayerTMMO}
    Spectra::T1
    Field::T2
    Bloch::T3
    Misc::T4
    AdmPhase::T5
    Beam::T6
    Layers::Array{T7}
end

## Data conditioning and main calculations.
function tmm_optics(
    beam::T1, layers::Array{T2,N2};
    emfflag::T3=false, h::T4=10, pbgflag::T3=false, λ0::T5=-eps(),
) where {T1<:PlaneWave, T2<:LayerTMMO, N2, T3<:Bool, T4<:Int64, T5<:Float64}
    # Check if λ0 ∈ λ
    λ0 = _checkλ0(λ0, beam)
    # Find beam.λ closest to λ0
    idxλ0 = findmin(abs.(beam.λ .- λ0))[2][1]
    # Build the sequence of n and d depending on the input
    nLen = length(vec(layers))
    d = Vector{Float64}(undef, nLen)
    nλ0 = Vector{Float64}(undef, nLen)
    layers_n = Matrix{ComplexF64}(undef, length(beam.λ), nLen)
    _build_arrays!(d, layers_n, nλ0, idxλ0, λ0, vec(layers), beam)
    # Get the parameters for beam
    beam_, _ = _get_beam_parameters(beam)
    # Call transfer matrix method
    if emfflag
        tmmout = _transfer_matrix_emf(layers_n, d, beam_, h)
    else
        tmmout = transfer_matrix(layers_n, d, beam_)
    end
    # Photonic band gap for crystals without defects
    if pbgflag & (nLen > 3)
        ω = Vector{Float64}(undef, length(beam_.λ))
        κp = Matrix{ComplexF64}(undef, (length(beam_.λ), length(beam_.θ)))
        κs = similar(κp)
        qz::Vector{Float64} = sin.(beam_.θ).*π./2.0 # parallel wavevector qz
        Λ::Float64 = sum(d[2:3])
        d1, d2 = d[2], d[3]
        n0, n1, n2 = nλ0[1], nλ0[2], nλ0[3]
        ωh::Float64 = Λ/π/(d1*n1 + d2*n2)*acos(-abs(n1 - n2)/(n1 + n2))
        ωl::Float64 = Λ/π/(d2*sqrt(n2^2 - n0^2) + d1*sqrt(n1^2 - n0^2))*acos(abs((n1^2*sqrt(n2^2 - n0^2) - n2^2*sqrt(n1^2 - n0^2))/(n1^2*sqrt(n2^2 - n0^2) + n2^2*sqrt(n1^2 - n0^2))))
        _photonic_dispersion!(κp, κs, ω, Λ, beam_, layers_n[idxλ0, 2:3], d[2:3])
    else
        κp = []; κs = []; ω = []; Λ = 0.0; ωl = 0.0; ωh = 0.0; qz = []
    end
    # Provide the multilayer depth considering the h division
    _ℓ::Matrix{Float64} = (d[2:end-1]/h)*ones.(1, h) # outer product
    ℓ = cumsum([0; _ℓ[:]], dims=1)[1:end-1] # remove last from cumsum
    # Return results
    sol = TMMOptics(
            tmmout[1],
            tmmout[2],
            Bloch(κp, κs, ω, Λ, ωl, ωh, qz),
            Misc(d, ℓ, nλ0, layers_n, λ0),
            tmmout[3],
            beam,
            layers
    )
    return sol
end

"""

    Check input λ0 againts beam.λ.

"""
function _checkλ0(λ0, beam::PlaneWave)
    x = isequal(λ0, -eps()) ? mean(beam.λ) : λ0
    return x
end

"""

    Build the sequence of indices of refraction and thicknesses depending on the input. It follows some logic depending on whether nλ0 was input.

"""
function _build_arrays!(
    d::Vector{T0},
    layers_n::Matrix{T1},
    nλ0::Vector{T0},
    idxλ0::T4,
    λ0::T0,
    layers::Vector{T2},
    beam::T5,
) where {T0<:Float64, T1<:ComplexF64, T2<:LayerTMMO, T4<:Int64, T5<:PlaneWave}
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
function transfer_matrix(
    layers_n::Matrix{T1}, d::Vector{T2}, beam::T3;
    λLen::T4=length(beam.λ), θLen::T4=length(beam.θ),
) where {T1<:ComplexF64, T2<:Float64, T3<:PlaneWave, T4<:Int64}
    # λ, λLen, θ, θLen = _get_beam_parameters(beam)
    nLen::Int64 = size(layers_n, 2)
    τs = Matrix{ComplexF64}(undef, (λLen, θLen))
    τp = similar(τs); ρs = similar(τs); ρp = similar(τs)
    δ = Array{ComplexF64,3}(undef, (λLen, θLen, nLen))
    ηs = similar(δ); ηp = similar(ηs)
    # Calculation of complex coefficients of reflection, transmission and emf
    _fresnel_coefficients!(ηs, ηp, δ, ρs, ρp, τs, τp, layers_n, d, beam, nLen)
    x = (
         Spectra(abs2.(ρp),
                 abs2.(ρs),
                 real(@view(ηp[:, :, 1]).*@view(ηp[:, :, end])).*abs2.(τp),
                 real(@view(ηs[:, :, 1]).*@view(ηs[:, :, end])).*abs2.(τs), ρp, ρs, τp, τs),
         Field([], []),
         AdmPhase(ηp, ηs, δ),
    )
    return x
end

"""

    Compute the Fresnel coefficients.

"""
function _fresnel_coefficients!(
    ηs::Array{T1,3}, ηp::Array{T1,3},
    δ::Array{T1,3},
    ρs::Matrix{T1}, ρp::Matrix{T1},
    τs::Matrix{T1}, τp::Matrix{T1},
    layers_n::Matrix{T1},
    d::Vector{T2},
    beam::T3,
    nLen::T4,
) where {T1<:ComplexF64, T2<:Float64, T3<:PlaneWave, T4<:Int64}
    cosϕ = Vector{ComplexF64}(undef, nLen)
    I2d = Matrix{ComplexF64}(I,2,2) # Identity 2x2 matrix
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
        Ψp::Matrix{ComplexF64} = reduce(*, [[I2d]; Φ.(@view(δ[l, a, 2:end-1]), @view(ηp[l, a, 2:end-1]))])
        Ψs::Matrix{ComplexF64} = reduce(*, [[I2d]; Φ.(@view(δ[l, a, 2:end-1]), @view(ηs[l, a, 2:end-1]))])
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
function _transfer_matrix_emf(
    layers_n::Matrix{T1}, d::Vector{T2}, beam::T3, h::T4;
    λLen::T4=length(beam.λ), θLen::T4=length(beam.θ),
) where {T1<:ComplexF64, T2<:Float64, T3<:PlaneWave, T4<:Int64}
    # λ, λLen, θ, θLen = _get_beam_parameters(beam)
    nLen::Int64 = size(layers_n, 2)
    τs = Matrix{ComplexF64}(undef, (λLen, θLen))
    τp = similar(τs); ρs = similar(τs); ρp = similar(τs)
    δ = Array{ComplexF64,3}(undef, (λLen, θLen, nLen))
    ηs = similar(δ); ηp = similar(ηs)
    emfs = Array{Float64,3}(undef, (λLen, θLen, (length(d) - 2)*h)); emfp = similar(emfs)
    # Calculation of complex coefficients of reflection, transmission and emf
    _fresnel_coefficients_emf!(ηs, ηp, δ, ρs, ρp, τs, τp, emfs, emfp, layers_n, d, beam, nLen, h)
    x = (
         Spectra(abs2.(ρp),
                 abs2.(ρs),
                 real(@view(ηp[:, :, 1]).*@view(ηp[:, :, end])).*abs2.(τp),
                 real(@view(ηs[:, :, 1]).*@view(ηs[:, :, end])).*abs2.(τs), ρp, ρs, τp, τs),
         Field(emfp, emfs),
         AdmPhase(ηp, ηs, δ),
    )
    return x
end

"""

    Compute the Fresnel coefficients and the electromagnetic field.

"""
function _fresnel_coefficients_emf!(
    ηs::Array{T1,3}, ηp::Array{T1,3},
    δ::Array{T1,3},
    ρs::Matrix{T1}, ρp::Matrix{T1},
    τs::Matrix{T1}, τp::Matrix{T1},
    emfs::Array{T2,3}, emfp::Array{T2,3},
    layers_n::Matrix{T1},
    d::Vector{T2},
    beam::T3,
    nLen::T4,
    h::T4,
) where {T1<:ComplexF64, T2<:Float64, T3<:PlaneWave, T4<:Int64}
    cosϕ = Vector{ComplexF64}(undef, nLen)
    I2d = Matrix{ComplexF64}(I, 2, 2) # Identity 2x2 matrix
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
        Ψp::Matrix{ComplexF64} = reduce(*, [[I2d]; Φ.(@view(δ[l, a, 2:end-1]), @view(ηp[l, a, 2:end-1]))])
        Ψs::Matrix{ComplexF64} = reduce(*, [[I2d]; Φ.(@view(δ[l, a, 2:end-1]), @view(ηs[l, a, 2:end-1]))])
        # Compute the Fresnell coefficients
        ρs[l, a] = ρ(ηs[l, a, 1], ηs[l, a, end], Ψs)
        ρp[l, a] = ρ(ηp[l, a, 1], ηp[l, a, end], Ψp)
        τs[l, a] = τ(ηs[l, a, 1], ηs[l, a, end], Ψs)
        τp[l, a] = τ(ηp[l, a, 1], ηp[l, a, end], Ψp)
        # Compute the _emfield
        emfs[l, a, :] .= _emfield(@view(layers_n[l, :]), d, @view(δ[l, a, :]), @view(ηs[l, a, :]), Ψs, nLen, h)
        emfp[l, a, :] .= _emfield(@view(layers_n[l, :]), d, @view(δ[l, a, :]), @view(ηp[l, a, :]), Ψp, nLen, h)
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
function ρ(η0::T1, ηm::T1, Ψ::Matrix{T1}) where {T1<:ComplexF64}
    r = (η0*Ψ[1,1] - Ψ[2,1] + η0*ηm*Ψ[1,2] - ηm*Ψ[2,2])/(η0*Ψ[1,1] + Ψ[2,1] + η0*ηm*Ψ[1,2] + ηm*Ψ[2,2])
    return r
end
# Transmission coefficient
function τ(η0::T1, ηm::T1, Ψ::Matrix{T1}) where {T1<:ComplexF64}
    t = 2.0/(η0*Ψ[1,1] + Ψ[2,1] + η0*ηm*Ψ[1,2] + ηm*Ψ[2,2])
    return t
end

"""

    Computes the inverse total transfer matrix for the whole structure at each wavelenth and angle of incidence and return the field.

"""
function _emfield(
    N::SubArray{T1,1},
    d::Vector{T2},
    δ::SubArray{T1,1},
    η::SubArray{T1,1},
    Ψ::Matrix{T1},
    nLen::T3,
    h::T3,
) where {T1<:ComplexF64, T2<:Float64, T3<:Int64}
    m0 = Matrix{ComplexF64}(undef, 2, 2)
    m1 = Matrix{ComplexF64}(I, 2, 2) # Identity 2x2 matrix
    g11 = Vector{ComplexF64}(undef, (nLen - 2)*h); g12 = similar(g11)
    # Divide the phase shift by h but keep η as is for each layer
    mδ::Array{ComplexF64} = δ/h
    _m1::Array = Ξ.(mδ, η)
    _g_elements!(g11, g12, m0, _m1, m1, Ψ, nLen, h)
    fi = _field_intensity(g11, g12, η[1], η[end], Ψ)
    return fi
end

"""

    Compute the elements of the G matrix.

"""
function _g_elements!(
    g11::Vector{T1}, g12::Vector{T1},
    m0::Matrix{T1},
    _m1::Vector{Matrix{T1}},
    m1::Matrix{T1},
    Ψ::Matrix{T1},
    nLen::T2,
    h::T2,
) where {T1<:ComplexF64, T2<:Int64}
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

"""

    Compute the field intensity.

"""
function _field_intensity(
    g11::Vector{T1}, g12::Vector{T1}, η1::T1, ηf::T1, Ψ::Matrix{T1},
) where {T1<:ComplexF64}
    fi = abs2.((g11 .+ ηf.*g12)./(0.25*(η1*Ψ[1,1] + Ψ[2,1] + η1*ηf*Ψ[1,2] + ηf*Ψ[2,2])))
    return fi
end

"""

    Calculates the inverse of optical transfer matrix of a layer.

        T = Ξ(φ, η)

            φ:  phase shift of the layer
            η: admittance of the layer

            T: 2x2 optical tranfer matrix.

"""
Ξ(φ::T1, η::T1) where {T1<:ComplexF64} = [cos(φ) (im*sin(φ)/η); (im*sin(φ)*η) cos(φ)]

"""

    Computes the photonic dispersion of binary structures (crystals only) alternating
    two different dielectric layers.

"""
function _photonic_dispersion!(
    κp::Matrix{T1}, κs::Matrix{T1},
    ω::Vector{T2},
    Λ::T2,
    beam::T3,
    n::Vector{T1},
    d::Vector{T2},
) where {T1<:ComplexF64, T2<:Float64, T3<:PlaneWave}
    @. ω = 2*π/beam.λ # Angular frequency
    # Angle of incidence of the second layer with Snell's law of cosine
    cosθ1::Vector{ComplexF64} = cos.(beam.θ)
    cosθ2::Vector{ComplexF64} = cosϑ.(n[1], n[2], cosθ1)
    # Prefactor for Bloch wavevector
    factor_s = _adm_factor.(ζs.(n[1], cosθ1), ζs.(n[2], cosθ2))
    factor_p = _adm_factor.(ζp.(n[1], cosθ1), ζp.(n[2], cosθ2))
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
_adm_factor(η1::T1, η2::T1) where {T1<:ComplexF64} = 0.5*(η1^2 + η2^2)/η1/η2

"""

    Bloch wavevector.

"""
cosκ(a1::T1, a2::T1, f::T1) where {T1<:ComplexF64} = cos(a1)*cos(a2) - f*sin(a1)*sin(a2)

"""

    Snell's law in cosine form. Returns the cosine already.

"""
function cosϑ(n1::T1, n2::T1, cosθ::T2) where {T1<:ComplexF64, T2<:Number}
    x = sqrt(1.0 - (n1/n2)^2*(1.0 - cosθ^2))
    return x
end

"""

    Admittance of the medium for p and s polarizations.

"""
ζp(n::T1, cosθ::T2) where {T1<:ComplexF64, T2<:Number} = n/cosθ
ζs(n::T1, cosθ::T2) where {T1<:ComplexF64, T2<:Number} = n*cosθ

end # module
