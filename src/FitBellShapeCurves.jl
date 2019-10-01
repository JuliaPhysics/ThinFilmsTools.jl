module FitBellShapeCurves

using LsqFit
using ..Utils: Gaussian1, Gaussian2, Gaussian3, Lorentzian1, Lorentzian2, Lorentzian3, Voigt1, Voigt2, Voigt3

export CurveModelSolution,
       FitCurveModel

abstract type FitShapeCurves end
struct FitCurveModel{T0, T1} <: FitShapeCurves where {T0<:Symbol, T1<:Float64}
    model::T0;
    xdata::Vector{T1}; ydata::Vector{T1};
    p::Vector{T1}; lb::Vector{T1}; ub::Vector{T1};
    solution; sigma::Vector{T1}; merror::Vector{T1}; confinter;
    ymodel::Vector{T1};
    alpha::T1
end

function FitCurveModel(model::T0, xdata::AbstractVector{T1}, ydata::AbstractVector{T1}, p::Vector{T1}; lb::Vector{T1}=-2.0.*abs.(seed), ub::Vector{T1}=2.0.*abs.(seed), α::T1=0.05) where {T0<:Symbol, T1<:Real}
    sol = curve_fit(eval(model), float.(xdata), float.(ydata), p, lower=lb, upper=ub)
    return FitCurveModel(model, float.(xdata), float.(ydata), float.(p), float.(lb), float.(ub), sol, stderror(sol), margin_error(sol, α), confidence_interval(sol, α), eval(model)(xdata, sol.param), float(α))
end

end # FitBellShapeCurves
