module ThreeOmegaMethod

using QuadGK

export ThreeOmega,
       LayerTOM,
       HeaterGeometry,
       Source

abstract type LayerTOMInformation end
struct LayerTOM{T1} <: LayerTOMInformation where {T1<:Float64}
    ky::T1
    kxy::T1
    d::T1
    ρC::T1
end

abstract type HeaterInformation end
struct HeaterGeometry{T1} <: HeaterInformation where {T1<:Float64}
    b::T1
    l::T1
    ρh::Array{T1}
end

abstract type HeaterSource end
struct Source{T1} <: HeaterSource where {T1<:Float64}
    p::T1
    f::Array{T1}
end

function ThreeOmega(
    Layers::Array{T1,N1},
    heater::T2,
    source::T3,
    thresistances::Array{T4,N4};
    int_limit::T5=1.0e6,
) where {T1<:LayerTOM, N1, T2<:HeaterGeometry, T3<:Source, T4<:Number, N4, T5<:Number}
    numLayers::Int64 = size(Layers, 1)
    # Check input
    numLayers > 1 || throw(DimensionMismatch("There should be at least one layer below the heater."))
    length(thresistances) == numLayers-1 || throw(DimensionMismatch("The number of interface-resistances should be lower than the number of layers: length(thresistances) == size(layers,2)-1."))
    # Power per unit length
    p_l::Float64 = source.p / heater.l
    # Reverse order of input for bottom-up calculation
    Layers = reverse(Layers, dims=1)
    thresistances = reverse(vec(thresistances), dims=1)
    # Lambda parameters
    λ = [lambda(Layers[i].ky, Layers[i].kxy, Layers[i].ky / Layers[i].ρC) for i = 1 : numLayers]
    # Construct B parameter for the system of LayerTOMs
    B = concatenateB(λ, Layers, thresistances, numLayers)
    # Construct integration term. The half term comes from the heater at the top of the stack
    F::Function = (η,ω) -> (1.0./λ[end](η,ω).*(B(η,ω)[1] + B(η,ω)[2])./(0.5*B(η,ω)[2] - 0.5*B(η,ω)[1]).*(sin.(heater.b*η)./(heater.b*η)).^2)
    # return integration
    sol = integrateTemperature(
            source.f,
            int_limit,
            F,
            heater.ρh,
            p_l/π/2,
            p_l/2/heater.b,
            im*4*π.*source.f,
    )
    return sol
end # ThreeOmegaMethod

"""Lambda parameter."""
lambda(ky::T1, kxy::T1, α::T1) where {T1<:Float64} = (η,ω) -> (ky*sqrt.(kxy*η^2 + im*2*ω/α))

"""Interface matrix."""
function Λ(λ1::T0, λ2::T0, r::T1) where {T0<:Function, T1<:Float64}
    x11 = (η,ω) -> λ2(η,ω) + λ1(η,ω) + λ2(η,ω).*λ1(η,ω).*r
    x12 = (η,ω) -> λ2(η,ω) - λ1(η,ω) - λ2(η,ω).*λ1(η,ω).*r
    x21 = (η,ω) -> λ2(η,ω) - λ1(η,ω) + λ2(η,ω).*λ1(η,ω).*r
    x22 = (η,ω) -> λ2(η,ω) + λ1(η,ω) - λ2(η,ω).*λ1(η,ω).*r
    X = (η,ω) -> (0.5./λ2(η,ω).*[x11(η,ω) x12(η,ω);
                                 x21(η,ω) x22(η,ω)])
    return X
end

"""Transfer matrix."""
function U(λ::T0, ky::T1, d::T1) where {T0<:Function, T1<:Float64}
    X = (η,ω) -> ([exp.(-λ(η,ω)/ky*d) 0.;
                   0. exp.(λ(η,ω)/ky*d)])
    return X
end

"""Integration of the temperature term."""
function integrateTemperature(
    f::Array{T0,N0},
    int_limit::T0,
    F::T1,
    ρh::Array{T0,N2},
    plπ::T0,
    plb::T0,
    h_param::Array{T2,N0},
) where {T0<:Float64, N0, T1<:Function, T2<:ComplexF64, N2}
    ΔTh = Array{ComplexF64,1}(undef, length(f))
    int_error = Array{Float64,1}(undef, length(f))
    @inbounds for i in eachindex(f)
        temp0::Tuple{ComplexF64,Float64} = quadgk((η) -> F(η, 2*π*f[i]), 0, int_limit, rtol=sqrt.(eps()))
        temp1::ComplexF64 = temp0[1]
        int_error[i] = temp0[2]
        ΔTh[i] = ( plπ*temp1 + ρh[1]*plb ) ./ ( 1 + ρh[2]*h_param[i] * (ρh[1] + temp1/plb)  )
    end
    sol = (f=f, ΔThreal=real(ΔTh), ΔThimag=imag(ΔTh), int_error=int_error)
    return sol
end # function integrateTemperature()

"""

    Concatenate B terms depending on the number of LayerTOMs. The boundary condition for
    the substrate is semi-infinite.

"""
function concatenateB(
    λ::Array{T0,N0},
    Layers::Array{T1,N1},
    thresistances::Array{T2,N2},
    numLayers::T3,
) where {T0<:Function, N0, T1<:LayerTOM, N1, T2<:Number, N2, T3<:Int64}
    B = Λ(λ[end-1], λ[end], thresistances[end])
    @inbounds for i = (numLayers-1):-1:2
        tempU::Function = U(λ[i], Layers[i].ky, Layers[i].d)
        tempΛ::Function = Λ(λ[i-1], λ[i], thresistances[i])
        B = let B = B
            (η,ω) -> B(η,ω) * tempU(η,ω) * tempΛ(η,ω)
        end # B = let B = B
    end # i = numLayers-1 : -1 : 2
    # Add semi-infinite boundary condition
    B = let B = B
        (η,ω) -> B(η,ω) * [0.; 1.]
    end # B = let B = B
    return B
end # function concatenateB()

end # module
