module CommonStructures

using Statistics

export ModelFit,
       BoundariesFit,
       PlaneWave,
       LayerTMMO,
       meanSquaredError,
       getBeamParameters,
       FitProcedure

## Defines the abstract type for all the fitting.
abstract type FitProcedure end

## General structure to build the indices of refraction models for otimization.
struct ModelFit{T1,T2,T3} <: FitProcedure where{T1<:Symbol, T2<:Real, T3<:Real}
    model::T1
    N # effective medium equations, EMA
    df::T2
    dfm::T1
    c::T3
    function ModelFit(
        model::T1;
        N=([],[]), df::T2=1.0/3.0, dfm::T1=:spheres, c::T3=0.0,
    ) where {T1<:Symbol, T2<:Real, T3<:Real}
        length(N) == 2 || throw("length(N) == 2, one vector for each index of refraction:
        N=(n1, n2), where n1 and n2 are vectors.")
        return new{T1,T2,T3}(model, N, df, dfm, c)
    end
end

## Definition of the boundaries for the solution space grids.
struct BoundariesFit{T1,T2,T3,T4,T5} <: FitProcedure where {T1<:Real, T2<:Real, T3<:Real, T4<:Real, T5<:Int64}
    odlo::T1
    odup::T2
    plo::T3
    pup::T4
    Nod::T5
    Np::T5
    function BoundariesFit(
        odlo::T1, odup::T2, plo::T3, pup::T4;
        Nod::T5=20, Np::T5=20,
    ) where {T1<:Real, T2<:Real, T3<:Real, T4<:Real, T5<:Int64}
        return new{T1,T2,T3,T4,T5}(odlo, odup, plo, pup, Nod, Np)
    end
end

## Definition of the beam type.
abstract type LightSource end
struct PlaneWave{T1,T2,T3} <: LightSource where {T1<:Real,T2<:Real,T3<:Real}
    λ::AbstractVector{T1}
    θ::AbstractVector{T2}
    p::T3
    function PlaneWave(
        λ::AbstractVector{T1}, θ::AbstractVector{T2};
        p::T3=0.5,
    ) where {T1<:Real,T2<:Real,T3<:Real}
        0.0<=p<=1.0 || throw("The polarisation, p argument, must lie between 0 and 1: 0.0<=p<=1.0")
        return new{T1,T2,T3}(λ, θ, p)
    end
end

## Material type definition with geometrical (physical) and optical thickness.
abstract type Material end
struct LayerTMMO{T1,T2,T3,T4} <: Material where {T1<:ComplexF64, T2<:Symbol, T3<:Real, T4<:Number}
    n::Array{T1,1}
    type::T2
    d::T3
    nλ0::Array{T4}
    function LayerTMMO(
        n::AbstractArray{T1,1};
        type::T2=:GT,
        d::T3=1.0/4.0,
        nλ0::Array{T4}=[-eps()],
    ) where {T1<:ComplexF64, T2<:Symbol, T3<:Real, T4<:Number}
        return new{T1,T2,T3,T4}(n, type, d, nλ0)
    end
end

"""

    Returns the mean squared error.

        mse = meanSquaredError(X, Xexp, σ)

            X: model spectrum
            Xexp: experimental spectrum to fit
            σ: std of Xexp

            mse: mean squared error

"""
function meanSquaredError(
    X::Array{T1}, Xexp::Array{T1};
    σ::Array{T1}=ones.(length(Xexp),1),
) where {T1<:Float64}
    mse = mean(abs2.((X .- Xexp)./σ))
    return mse
end

## Returns the parameters with the correct format for the beam suitable for the calculations.
function getBeamParameters(beam::PlaneWave)
    beam_ = PlaneWave(float.(beam.λ), deg2rad.(beam.θ))
    x = (λLen=length(beam.λ), θLen=length(beam.θ))
    return beam_, x
end

end # module
