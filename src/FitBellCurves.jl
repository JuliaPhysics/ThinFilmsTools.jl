module FitBellCurves

using Optim
using ..Utils: gaussian, lorentzian, voigtian, flattenArrays, arrayArrays
using ..CommonStructures: meanSquaredError

export FitCurveModel

function FitCurveModel(
    model::T0,
    xdata::AbstractVector{T1},
    ydata::AbstractVector{T1},
    seed::Array{Array{T1,1},1};
    lb::Array{Array{T1,1},1}=-2.0.*arrayArrays((abs.(flattenArrays(seed)).+1.0), seed),
    ub::Array{Array{T1,1},1}=2.0.*arrayArrays((abs.(flattenArrays(seed)).+1.0), seed),
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
    sol = optimize(p->fitMSE(p, model, seed, xdata, ydata, σ),
                   flattenArrays(lb), flattenArrays(ub), flattenArrays(seed), alg, options)
    xopt = arrayArrays(sol.minimizer, seed)
    sol_ = (
        model = model,
        xdata = float.(xdata),
        ydata = float.(ydata),
        seed = float.(seed),
        lb = float.(lb),
        ub = float.(ub),
        solution = sol,
        ymodel = eval(model)(xdata, xopt),
        optParams = xopt,
        sigma = σ,
    )
    return sol_
end

function fitMSE(p, model, seed, xdata, ydata, σ)
    ŷ = eval(model)(xdata, arrayArrays(p, seed))
    mse = meanSquaredError(vec(ŷ), ydata; σ=σ)
    return mse
end

end # module
