module CommonStructures

using Statistics
using ..RI

export ModelFit,
       BoundariesFit,
       PlaneWave,
       LayerTMMO1DIso,
       refractiveIndexMR,
       checkInput,
       checkInputAlpha,
       multilayerParameters!,
       buildArraysBragg,
       buildArraysFP1d,
       monotonicLin!,
       meanSquaredError,
       getBeamParameters

# Defined the abstract type for all the fitting
abstract type FitProcedure end

"""General structure to build the indices of refraction models for otimization."""
struct ModelFit{T1} <: FitProcedure where{T1<:Symbol}
    model::T1;
    N # effective medium equations, EMA
    function ModelFit(model::T1; N=([],[])) where {T1<:Symbol}
        length(N) == 2 || throw("length(N) == 2, one vector for each index of refraction: N=(n1, n2), where n1 and n2 are vectors.")
        return new{T1}(model, N)
    end
end

"""Definition of the boundaries for the solution space grids."""
struct BoundariesFit{T1, T2, T3, T4, T5} <: FitProcedure where {T1<:Real, T2<:Real, T3<:Real, T4<:Real, T5<:Int64}
    odlo::T1; odup::T2; plo::T3; pup::T4; Nod::T5; Np::T5
    function BoundariesFit(odlo::T1, odup::T2, plo::T3, pup::T4; Nod::T5=20, Np::T5=20) where {T1<:Real, T2<:Real, T3<:Real, T4<:Real, T5<:Int64}
        return new{T1, T2, T3, T4, T5}(odlo, odup, plo, pup, Nod, Np)
    end
end

"""Definition of the beam type."""
abstract type LightSource end
struct PlaneWave{T1, T2, T3} <: LightSource where {T1<:Real, T2<:Real, T3<:Real}
    λ::AbstractVector{T1}; θ::AbstractVector{T2}; p::T3
    function PlaneWave(λ::AbstractVector{T1}, θ::AbstractVector{T2}; p::T3=0.5) where {T1<:Real, T2<:Real, T3<:Real}
        0.0<=p<=1.0 || throw("The polarisation, p argument, must lie between 0 and 1: 0.0<=p<=1.0")
        return new{T1, T2, T3}(λ, θ, p)
    end
end

"""Material type definition with geometrical (physical) and optical thickness."""
abstract type Material end
struct LayerTMMO1DIso{T1, N1, T2, T3, T4} <: Material where {T1<:ComplexF64, N1, T2<:Symbol, T3<:Float64, T4<:Number}
    n::Array{T1,1}; type::T2; d::T3; nλ0::Array{T4}
    function LayerTMMO1DIso(n::AbstractArray{T1,N1}; type::T2=:GT, d::T3=0.0, nλ0::Array{T4}=[-eps()]) where {T1<:ComplexF64, N1, T2<:Symbol, T3<:Float64, T4<:Number}
        return new{T1, N1, T2, T3, T4}(n, type, d, nλ0)
    end
end

"""

    Returns the index of refraction for the selected model. RI is imported from RefractiveIndicesModels.jl.

        N = refractiveIndexMR(layer, argin, λ)

            layer: input layers with ModelFit structure
            argin: input with parameters to optimize
            λ: input range of wavelengths [nm] to calculate the index of refraction of non-EMA models

"""
function refractiveIndexMR(layer::T1, argin::AbstractArray{T2}, λ::Array{T2,1}) where {T1<:ModelFit, T2<:Float64}
    if layer.model == :bruggemanspheresbin
        errorParametersEMA(layer)
        return RI.bruggemanspheresbin(argin[1], layer.N[1], layer.N[2])
    elseif layer.model == :looyengacylindersbin
        errorParametersEMA(layer)
        return RI.looyengacylinders(argin[1], layer.N[1], layer.N[2])
    elseif layer.model == :looyengaspheresbin
        errorParametersEMA(layer)
        return RI.looyengaspheresbin(argin[1], layer.N[1], layer.N[2])
    elseif layer.model == :moneckebin
        errorParametersEMA(layer)
        return RI.moneckebin(argin[1], layer.N[1], layer.N[2])
    elseif layer.model == :maxwellgarnettspheresbin
        errorParametersEMA(layer)
        return RI.maxwellgarnettspheresbin(argin[1], layer.N[1], layer.N[2])
    elseif layer.model == :gedfbin
        errorParametersEMA(layer)
        return RI.gedfbin(argin, layer.N[1], layer.N[2])
    elseif layer.model == :gembin
        errorParametersEMA(layer)
        return RI.gembin(argin, layer.N[1], layer.N[2])
    elseif layer.model == :sellmeier
        return RI.sellmeier(argin, λ./1e3) # Sellmeier eq takes λ in μm
    elseif layer.model == :cauchyurbach
    	return RI.cauchyurbach(argin, λ)
    elseif layer.model == :lorentzdispersion
    	return RI.lorentzdispersion(argin, λ)
    elseif layer.model == :tauclorentz
    	return RI.tauclorentz(argin, λ)
    elseif layer.model == :forouhibloomer
    	return RI.forouhibloomer(argin, λ)
    elseif layer.model == :forouhibloomermodified
    	return RI.forouhibloomermodified(argin, λ)
    end
end

"""

    Error if EMA models do not have the two columns indices of refraction.

"""
function errorParametersEMA(layer::ModelFit)
    !|(isempty(layer.N[1]), isempty(layer.N[2])) || throw("To use the $(layer.model) model you need to specify the two indices of refraction of the components in ModelFit.")
    return nothing
end

"""

    Check the input sizes and the seed with alpha.

"""
function checkInput(xinit, layers)
    length(xinit) == length(findall(isa.(layers, ModelFit))) || throw("There should be seeds for all active layers (ModelFit).")
    return nothing
end

function checkInputAlpha(xinit, layers)
    length(xinit) == length(findall(isa.(layers, ModelFit)))+1 || throw("There should be seeds for all active layers (ModelFit) plus the alpha in the last one.")
    return nothing
end

"""

    Construct the multilayers parameters to optimise depending on their type (LayerTMMO1DIso, ModelFit).

"""
function multilayerParameters!(N::Array{T1,2}, d::Array{T2,1}, x::Array{Array{T2,1},1}, beam::T3, layers::Array) where{T1<:ComplexF64, T2<:Float64, T3<:PlaneWave}
    k::Int64 = 1
    for i in eachindex(vec(layers))
        if isa(layers[i], LayerTMMO1DIso)
            d[i] = layers[i].d
            @views N[:, i] = layers[i].n
        elseif isa(layers[i], ModelFit)
            d[i] = x[k][1] # the thickness is defined first
            @views N[:, i] = refractiveIndexMR(layers[i], x[k][2:end], beam.λ)
            k += 1
        end
    end
    return nothing
end

"""

    Build arrays according to the number of periods.

"""
function buildArraysBragg(N::Array{T1,2}, d::Array{T2,1}, L::T3) where {T1<:ComplexF64, T2<:Float64, T3<:Int64}
    N = hcat(@view(N[:, 1]), repeat([@view(N[:, 2]) @view(N[:, 3])], 1, L), @view(N[:, end]))
    d = hcat(d[1], repeat([d[2] d[3]], 1, L), d[end])
    return N, vec(d)
end

"""

    Build arrays according to the number of periods and one defect.

"""
function buildArraysFP1d(N::Array{T1,2}, d::Array{T2,1}, L::Array{T3,1}, Ld::Array{T3,1}) where {T1<:ComplexF64, T2<:Float64, T3<:Int64}
    N = hcat(@view(N[:, 1]), repeat([@view(N[:, 2]) @view(N[:, 3])], 1, L[1]), repeat(@view(N[:, 2]), 1, Ld[1]), repeat([@view(N[:, 3]) @view(N[:, 2])], 1, L[2]), @view(N[:, end]))
    d = hcat(d[1], repeat([d[2] d[3]], 1, L[1]), repeat([d[2]], 1, Ld[1]), repeat([d[3] d[2]], 1, L[2]), d[end])
    return N, vec(d)
end

"""

    Linear decrease of thickness.

"""
function monotonicLin!(d::Array{T1,1}, α::T1) where {T1<:Float64}
    # d .*= α.^(1:length(d))
    # Only for the second active layer until the one before the substrate
    d .*= [1.0; 1.0; α.^(3:length(d)-1); 1.0]
    return nothing
end

"""

    Returns the mean squared error.

        mse = meanSquaredError(X, Xexp, σ)

            X: model spectrum
            Xexp: experimental spectrum to fit
            σ: std of Xexp
            mse: mean squared error

"""
function meanSquaredError(X::Array{T1}, Xexp::Array{T1}; σ::Array{T1}=ones.(length(Xexp), 1)) where {T1<:Float64}
    return mean(abs2.((X .- Xexp)./σ))
end

"""

    Returns the parameters with the correct format for the beam suitable for the calculations.

        beam_, x = getBeamParameters(beam::PlaneWave)

            beam: Planewave structure
            beam_: Planewave structure with λ as float type, and θ in radians as float type
            x: NamedTuple with
                x.λLen: length(λ)
                x.θLen: length(θ)

"""
function getBeamParameters(beam::PlaneWave)
    return PlaneWave(float.(beam.λ), deg2rad.(beam.θ)), (λLen=length(beam.λ), θLen=length(beam.θ))
end

end # module
