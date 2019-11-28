module FitBellCurves

using Optim
using ..Utils: gaussian, lorentzian, voigtian, flatten_arrays, array_arrays
using ..FTFSpectrum: _objective_function, MeanAbs

abstract type FitShapeCurves end
struct FitCurveModel{T0, T1} <: FitShapeCurves where {T0<:Symbol, T1<:Float64}
    model::T0
    xdata::Vector{T1}
    ydata::Vector{T1}
    seed::Array{Array{T1,1},1}
    lb::Array{Array{T1,1},1}
    ub::Array{Array{T1,1},1}
    solution
    ymodel::Vector{T1}
    optParams::Array{Array{T1,1},1}
    sigma::Array{T1,1}
end

export FitCurveModel, fit_curve_model

function fit_curve_model(
    model::T0,
    xdata::AbstractVector{T1},
    ydata::AbstractVector{T1},
    seed::Array{Array{T1,1},1};
    lb::Array{Array{T1,1},1}=-2.0.*array_arrays((abs.(flatten_arrays(seed)).+1.0), seed),
    ub::Array{Array{T1,1},1}=2.0.*array_arrays((abs.(flatten_arrays(seed)).+1.0), seed),
    σ::Array{T1,1}=ones(length(ydata)),
    alg=SAMIN(),
    options=Optim.Options(
            g_abstol=1e-8,
            g_reltol=1e-8,
            iterations=10^5,
            show_trace=false,
            store_trace=true,
    ),
) where {T0<:Symbol, T1<:Real}
    isequal(length(xdata), length(ydata)) || throw("Input xdata and ydata must have the same length.")
    oftype = MeanAbs()
    sol = optimize(p -> _fit_mse(p, model, seed, xdata, ydata, σ, oftype),
                   flatten_arrays(lb), flatten_arrays(ub), flatten_arrays(seed), alg, options)
    xopt = array_arrays(sol.minimizer, seed)
    sol_ = FitCurveModel(
                model,
                float.(xdata),
                float.(ydata),
                float.(seed),
                float.(lb),
                float.(ub),
                sol,
                eval(model)(xdata, xopt),
                xopt,
                σ,
    )
    return sol_
end

function _fit_mse(p, model, seed, xdata, ydata, σ, oftype)
    ŷ = eval(model)(xdata, array_arrays(p, seed))
    mse = _objective_function(oftype, vec(ŷ), ydata; σ=σ)
    return mse
end

end # module
