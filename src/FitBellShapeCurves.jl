module FitBellShapeCurves

using Optim
using ..Utils: Gaussian, Lorentzian, Voigtian, flattenArrays, arrayArrays
using ..CommonStructures: meanSquaredError

export FitCurveModel

abstract type FitShapeCurves end
struct FitCurveModel{T0, T1} <: FitShapeCurves where {T0<:Symbol, T1<:Float64}
    model::T0;
    xdata::Vector{T1}; ydata::Vector{T1};
    seed::Array{Array{T1,1},1}; lb::Array{Array{T1,1},1}; ub::Array{Array{T1,1},1};
    solution; ymodel::Vector{T1}; optParams::Array{Array{T1,1},1}; sigma::Array{T1,1}
end

function FitCurveModel(model::T0, xdata::AbstractVector{T1}, ydata::AbstractVector{T1}, seed::Array{Array{T1,1},1}; lb::Array{Array{T1,1},1}=-2.0.*arrayArrays((abs.(flattenArrays(seed)).+1.0), seed), ub::Array{Array{T1,1},1}=2.0.*arrayArrays((abs.(flattenArrays(seed)).+1.0), seed), σ::Array{T1,1}=ones(length(ydata)), alg=SAMIN(), options=Optim.Options(g_abstol=1e-8, g_reltol=1e-8, iterations=10^5, show_trace=false, store_trace=true)) where {T0<:Symbol, T1<:Real}
    isequal(length(xdata), length(ydata)) || throw("Input xdata and ydata must have the same length.")
    sol = optimize(p->fitMSE(p, model, seed, xdata, ydata, σ),
                   flattenArrays(lb), flattenArrays(ub), flattenArrays(seed), alg, options)
    xopt = arrayArrays(sol.minimizer, seed)
    return FitCurveModel(model, float.(xdata), float.(ydata), float.(seed), float.(lb), float.(ub), sol, eval(model)(xdata, xopt), xopt, σ)
end

function fitMSE(p, model, seed, xdata, ydata, σ)
    ŷ = eval(model)(xdata, arrayArrays(p, seed))
    return meanSquaredError(vec(ŷ), ydata; σ=σ)
end

end # FitBellShapeCurves
