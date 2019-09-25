module FitThinFilmPattern

using Optim
using Statistics
using ..CommonStructures: ModelFit, BoundariesFit, PlaneWave, LayerTMMO1DIso, refractiveIndexMR, checkInput, checkInputAlpha, multilayerParameters!, buildArraysBragg, buildArraysFP1d, monotonicLin!, meanSquaredError, getBeamParameters
using ..Utils: build_interpolation, averagePolarisation, flattenArrays, arrayArrays
using ..MatrixMethod1DIsotropic: TransferMatrix

export SpaceSolutionEMA,
       TheoreticalSpectrum,
       NormalizeReflectance,
       FitTMMO1DIsotropic,
       Reflectance,
       Transmittance,
       Generic,
       DBR,
       DBRAlpha,
       MC1d,
       MC1dAlpha

abstract type FitProcedure end

"""Results of the solution space parameters."""
struct SpaceSolution2D{T1, T2} <: FitProcedure where {T1<:Float64, T2<:PlaneWave}
    solSpace::Array{T1}; objfunMin::T1; optOD::T1; optParams::T1; optThickness::T1; spectrumFit::Array; spectrumExp::Array; x1; x2; beam::T2
end

struct SpaceSolutionEMA{T1, T2} <: FitProcedure where {T1<:Float64, T2<:PlaneWave}
    solSpace::Array{T1}; objfunMin::T1; optOD::T1; optParams::T1; optThickness::T1; spectrumFit::Array; spectrumExp::Array; od; p; beam::T2
end

"""Solution of the optimization with the Optim solver for multilayers."""
struct FitTMMO1DIsotropic{T1, T2, T3} <: FitProcedure where {T1<:Float64, T2<:PlaneWave, T3<:NamedTuple}
    spectrumFit::Array{T1}; spectrumExp::Array{T1}; optParams::Array{Array{T1,1},1}; objfunMin::T1; beam::T2; layers::Array; fitOptions::T3
end

"""Definition of subtypes to avoid repitition of functions."""
struct FitReflectance <: FitProcedure end
Reflectance() = FitReflectance()
struct FitTransmittance <: FitProcedure end
Transmittance() = FitTransmittance()
struct FitGeneric <: FitProcedure end
Generic() = FitGeneric()
struct FitDBR <: FitProcedure end
DBR() = FitDBR()
struct FitDBRAlpha <: FitProcedure end
DBRAlpha() = FitDBRAlpha()
struct FitMC1d <: FitProcedure end
MC1d() = FitMC1d()
struct FitMC1dAlpha <: FitProcedure end
MC1dAlpha() = FitMC1dAlpha()

"""

    Normalize the experimental reflectance dividing by the reference reflectance and multiplying by the theoretical one. The experimental and the reference reflectances must have the same scale.

        R = NormalizeReflectance(λ, Rexp, Rthe, Rref)

            λ: wavelength of interest [nm]
            Rexp: 2-columns array with wavelength in first one and experimental reflectance in the second
            Rthe: 2-columns array with wavelength in first one and theoretical reflectance in the second
            Rref: 2-columns array with wavelength in first one and reference reflectance in the second

    The output has the same scale as the theoretical reflectance with the same length as λ.

"""
function NormalizeReflectance(λ::AbstractArray{T0,1}, Rexp::Array{T1,2}, Rthe::Array{T1,2}, Rref::Array{T1,2}) where {T0<:Real, T1<:Float64}
    # Interpolate for the λ range
    spl_exp = build_interpolation(Rexp)
    spl_ref = build_interpolation(Rref)
    spl_the = build_interpolation(Rthe)
    exp2ref = spl_exp.(λ) ./ spl_ref.(λ)
    # Find out if the sample is overshooting the reference
    (maximum(exp2ref) < 1.0) || return (exp2ref./ maximum(exp2ref) .* spl_the.(λ))
    return exp2ref .* spl_the.(λ)
end

"""

    Returns the polarisation averaged theoretical spectrum calculated with the transfer matrix method. It consiers light hitting the first medium (incidentRI).

        X = TheoreticalSpectrum(specType, beam, incidentRI, emergentRI)

            specType: type of spectrum to compute, Reflectance() or Transmittance()
            beam: a PlaneWave structure
            incidentRI: Array containing the index of refraction of the incident material for several wavelengths. For instance, RIdb.air(beam.λ).
            emergentRI: Array containing the index of refraction of the emergent material for several wavelengths. For instance, RIdb.silicon(beam.λ).
            X: averaged spectrum = beam.p*Xp + (1.0 - beam.p)*Xs

"""
function TheoreticalSpectrum(specType::T0, beam::T1, N1::Array{T2,1}, N2::Array{T2,1}) where {T0<:FitProcedure, T1<:PlaneWave, T2<:ComplexF64}
    beam_, _ = getBeamParameters(beam)
    return computeTransferMatrix(specType, vec([0. 100.0 0.]), [N1 N2 N2], beam_)
end

"""

    Computes the solution space for a window range of two parameters, the optical thickness and porosity. Only works with these two parameters and for a single layer between incident and emergent media.

    It uses the mean(abs(Xexp - X)) objective function, where Xexp and X are the experimental and modelled spectra.

        sol = SpaceSolutionEMA(specType, b, beam, Xexp, layers)

            specType: type of spectrum to compute, Reflectance() or Transmittance()
            b: lower and upper boundaries for the optical thickness and porosity (e.g. [odl, odu,
               pl pu; Nod, Np]), where Nod and Np indicates the number of points to search in the solution space as optional parameters
            beam: structure from PlaneWave
            Xexp: absolute experimental spectrum
            layers: columnwise array of three LayerTMMO1dIso and ModelFit with information about the layers.
                    For SpaceSolutionEMA: length(vec(layers)) = 3
                        The typeof first and third layers = LayerTMMO1DIso
                        The typeof second layers = ModelFit
            sol: structure with the solution
                solSpace: matrix with objective function values for the whole space
                objfunMin: minimum value of the objective function found
                optOD: optimum value of optical thickness
                optParams: optimum value of the other parameter
                optThickness: optimum value of thickness
                spectrumFit: optimum spectrum obtained with dopt and xopt
                spectrumExp: experimental spectrum
                od: range of optical thicknesses explored
                p: range of porosities explored
                beam: information of the light used

"""
function SpaceSolutionEMA(specType::T0, b::T1, beam::T2, Xexp::Array{T3}, layers::Array) where{T0<:FitProcedure, T1<:BoundariesFit, T2<:PlaneWave, T3<:Float64}
    # Generate grids
    _d = LinRange(b.odlo, b.odup, b.Nod)
    _p = LinRange(b.plo, b.pup, b.Np)
    d = similar(_d)
    solSpace = zeros(Float64, b.Nod, b.Np)
    beam_, _ = getBeamParameters(beam)
    Xexp = vec(Xexp)
    @inbounds for j in eachindex(_p), i in eachindex(_d)
        neff = refractiveIndexMR(layers[2], [_p[j]], beam_.λ)
        d[i] = _d[i] / mean(real.(neff))
        X = computeTransferMatrix(specType, vec([0. d[i] 0.]), [layers[1].n neff layers[3].n], beam_)
        solSpace[i, j] = meanSquaredError(vec(X), Xexp)
    end
    smin = findmin(solSpace)
    X = computeTransferMatrix(specType, vec([0. d[smin[2][1]] 0.]), [layers[1].n refractiveIndexMR(layers[2], [_p[smin[2][2]]], beam_.λ) layers[3].n], beam_)
    return SpaceSolutionEMA(solSpace, smin[1], _d[smin[2][1]], _p[smin[2][2]], d[smin[2][1]], X, Xexp, _d, _p, beam_)
end

"""

    Optimize the parameters of the thin films in a multilayer stack using the selected models and the experimental spectrum.

        sol = FitTMMO1DIsotropic(specType, xinit, beam, Rexp, layers; options=Optim.Options(), alg=NelderMead(), lb=0.5.*xinit, ub=1.5.*xinit, arrange=Generic(), L=[1], Ld=[1], σ=ones.(length(beam.λ), 2))

            specType: type of spectrum to fit, Reflectance() or Transmittance()
            xinit: array with the initial parameters for optimization. For instance, to optimise two layers in the system you need to wrap the seeds for the algorithm as follow: xinit = [seed1, seed2], even if it is only one single layer to fit, xinit = [seed1]. The first parameter of each seed must be always the thickness, the rest are the parameters for the selected model. The number of seeds must match that of the number of ModelFit layers.
            beam: structure from PlaneWave
            Xexp: absolute experimental spectrum
            layers: columnwise array of LayerTMMO1dIso and ModelFit with information about the layers.
                σ: array with the standard deviation of the spectrum for each wavelength, by default is ones.(length(beam.λ), 1).
                options: optional Optim.Options structure
                alg: optional algorithm method selected, by default takes SAMIN(). You can pass the options inside as well, for instance, SAMIN(rt=0.1). Right now, NelderMead() and SAMIN() are supported from Optim.jl. If you select SAMIN() you need to input lb and ub, otherwise will be set as 0.5 and 1.5 times the xinit argument, respectively.
                lb: lower bounds for the optimisation variables, by default lb=0.5.*xinit
                ub: upper bounds for the optimisation variables, by default ub=1.5.*xinit
                arange: type of structure considered for the fitting procedure. Posible values:
                    Generic(): the general way of fitting, default
                    DBR(): custom optional optimisation of periodic structures, a.k.a, Bragg reflectors = FitDBR()
                    DBRAlpha(): custom optional optimisation of periodic structures, a.k.a, Bragg reflectors, with monotonic decrease in the thicknesses of the layers = FitDBRAlpha()
                    MC(): custom optional optimisation of Fabry-Perot (MicroCavity) type strucutures, i.e., a DBR with a defect in the center = FitMC1d()
                    MCAlpha(): custom optional optimisation of Fabry-Perot (MicroCavity) type strucutures, i.e., a DBR with a defect in the center, with monotonic decrease in the thicknesses of the layers = FitMC1dAlpha()
                L: period of each of the mirrors. In the case stack=:bragg this array has only one element while for stack=:fp1d with one defect, this has two elements.
                Ld: number of times the defect(s) layers are repeated. The 4th element of layers argument is considered to be the first defect layer, the fifth element of layers is the second defect and so on.
            sol: structure with results, with fields
                OptimSolverInfo: status info from the Optim solver
                spectrumFit: Spectrum obtained with the model and the optimal parameters
                spectrumExp: Input experimental spectrum
                optParams: array of arrays with optimal parameters
                objfunMin: optimum value of objective function MSE
                beam: structure of the light used
                fitOptions: NamedTuple with information of the fitting procedure

"""
function FitTMMO1DIsotropic(specType::T0, xinit::Array{Array{T5,1},1}, beam::T1, Xexp::Array{T2}, layers::Array; σ::Array{T2}=ones.(length(beam.λ)), options=Optim.Options(), alg=SAMIN(), lb::Array{Array{T6,1},1}=0.5.*xinit, ub::Array{Array{T6,1},1}=1.5.*xinit, arrange::T3=Generic(), L::Array{T4}=[1], Ld::Array{T4}=[1]) where{T0<:FitProcedure, T1<:PlaneWave, T2<:Float64, T3<:FitProcedure, T4<:Int64, T5<:Real, T6<:Real}
    if isa(arrange, FitDBR) || isa(arrange, FitDBRAlpha)
        length(vec(layers)) == 4 || throw("Bragg stack option accepts only 4 layers: incident medium, first lactive layer, second active layer and emergent medium.")
    end
    if isa(arrange, FitMC1d) || isa(arrange, FitMC1dAlpha)
        length(vec(layers)) == 4 || throw("Fabry-Perot stack option accepts only 4 layers: incident medium, first lactive layer, second active layer and the emergent medium. The defects layers are composed by the first active layer.")
        length(L) == 2 || throw("Fabry-Perot stack option accepts only 2 mirrors (length(L)=2) and 1 defect (length(Ld)=1).")
    end
    # Get the parameters for beam
    beam_, _ = getBeamParameters(beam)
    N = Array{ComplexF64,2}(undef, length(beam_.λ), length(vec(layers)))
    d = Array{Float64,1}(undef, length(vec(layers)))
    # Run the solver depending on the algorithm
    solution = runFittingProcedure(alg, arrange, specType, float.(xinit), beam_, vec(Xexp), σ, layers, options, float.(lb), float.(ub), d, N; L=Integer.(vec(L)), Ld=Integer.(vec(Ld)))
    return FitTMMO1DIsotropic(solution.X, Xexp, solution.xfinal, solution.solmin, beam_, layers, (fitType=specType, arrange=arrange, optimAlgorithm=Symbol(summary(solution.sol)), optimSolverInfo=solution.sol))
end

"""

    Run the fitting procedure.

        results = runFittingProcedure(alg, arrange, specType, xinit, beam, Xexp, σ, layers, options, lb, ub, d, N; kwargs...)

            alg: optional algorithm method selected
            arrange: type of fitting procedure to compute
            specType: type of spectrum to compute
            xinit: array with the initial parameters for optimization. For instance, to optimise two layers in the system you need to wrap the seeds for the algorithm as follow: xinit = [seed1, seed2], even if it is only one single layer to fit, xinit = [seed1]. The first parameter of each seed must be always the thickness, the rest are the parameters for the selected model. The number of seeds must match that of the number of ModelFit layers.
            beam: structure from PlaneWave
            Xexp: absolute experimental spectrum
            σ: std of the experimental spectrum
            layers: columnwise array of LayerTMMO1dIso and ModelFit with information about the layers.
            options: Optim.Options() structure
            d: thicknesses of the whole system to be updated
            N: index of refraction of the whole system to be updated
            kwargs: depending on the type of arrange
                L: period of each of the mirrors. For arrange=DBR()=DBRAlpha(), this array has only one element while for arrange=MC1d()=MC1dAlpha(), this has two elements.
                Ld: number of times the defect(s) layers are repeated. The 4th element of layers argument is considered to be the first defect layer, the fifth element of layers is the second defect and so on.
            results: NamedTuple with fields
                sol: Optim solver information
                X: modelled spectrum obtained with the optimal parameters
                xfinal: optimal parameters
                solmin: minimum of the objective function

"""
function runFittingProcedure end

# SAMIN, Generic
function runFittingProcedure(alg::T0, arrange::T1, specType::T2, xinit::Array{Array{T4,1},1}, beam::T3, Xexp::Array{T4,1}, σ::Array{T4,1}, layers::Array, options, lb::Array{Array{T4,1},1}, ub::Array{Array{T4,1},1}, d::Array{T4,1}, N::Array{T5,2}; kwargs...) where {T0<:SAMIN, T1<:FitGeneric, T2<:FitProcedure, T3<:PlaneWave, T4<:Float64, T5<:ComplexF64}
    checkInput(xinit, layers)
    sol = optimize(x->fitMSE(x, arrange, specType, beam, Xexp, σ, layers, xinit, d, N), flattenArrays(lb), flattenArrays(ub), flattenArrays(xinit), alg, options)
    xfinal = arrayArrays(sol.minimizer, xinit)
    multilayerParameters!(N, d, xfinal, beam, layers)
    X = computeTransferMatrix(specType, d, N, beam)
    return (sol=sol, X=X, xfinal=xfinal, solmin=sol.minimum)
end

# SAMIN, DBR
function runFittingProcedure(alg::T0, arrange::T1, specType::T2, xinit::Array{Array{T4,1},1}, beam::T3, Xexp::Array{T4,1}, σ::Array{T4,1}, layers::Array, options, lb::Array{Array{T4,1},1}, ub::Array{Array{T4,1},1}, d::Array{T4,1}, N::Array{T5,2}; L::Array{T6,1}=[1], kwargs...) where {T0<:SAMIN, T1<:FitDBR, T2<:FitProcedure, T3<:PlaneWave, T4<:Float64, T5<:ComplexF64, T6<:Int64}
    checkInput(xinit, layers)
    sol = optimize(x->fitMSE(x, arrange, specType, beam, Xexp, σ, layers, xinit, d, N, L[1]), flattenArrays(lb), flattenArrays(ub), flattenArrays(xinit), alg, options)
    xfinal = arrayArrays(sol.minimizer, xinit)
    multilayerParameters!(N, d, arrayArrays(sol.minimizer, xfinal), beam, layers)
    N, d = buildArraysBragg(N, d, L[1])
    X = computeTransferMatrix(specType, d, N, beam)
    return (sol=sol, X=X, xfinal=xfinal, solmin=sol.minimum)
end

# SAMIN, DBRAlpha
function runFittingProcedure(alg::T0, arrange::T1, specType::T2, xinit::Array{Array{T4,1},1}, beam::T3, Xexp::Array{T4,1}, σ::Array{T4,1}, layers::Array, options, lb::Array{Array{T4,1},1}, ub::Array{Array{T4,1},1}, d::Array{T4,1}, N::Array{T5,2}; L::Array{T6,1}=[1], kwargs...) where {T0<:SAMIN, T1<:FitDBRAlpha, T2<:FitProcedure, T3<:PlaneWave, T4<:Float64, T5<:ComplexF64, T6<:Int64}
    checkInputAlpha(xinit, layers)
    sol = optimize(x->fitMSE(x, arrange, specType, beam, Xexp, σ, layers, xinit, d, N, L[1]), flattenArrays(lb), flattenArrays(ub), flattenArrays(xinit), alg, options)
    xfinal = arrayArrays(sol.minimizer[1:end-1], xinit[1:end-1])
    push!(xfinal, [sol.minimizer[end]])
    multilayerParameters!(N, d, arrayArrays(sol.minimizer[1:end-1], xfinal[1:end-1]), beam, layers)
    N, d = buildArraysBragg(N, d, L[1])
    monotonicLin!(d, sol.minimizer[end])
    X = computeTransferMatrix(specType, d, N, beam)
    return (sol=sol, X=X, xfinal=xfinal, solmin=sol.minimum)
end

# SAMIN, MC1d
function runFittingProcedure(alg::T0, arrange::T1, specType::T2, xinit::Array{Array{T4,1},1}, beam::T3, Xexp::Array{T4,1}, σ::Array{T4,1}, layers::Array, options, lb::Array{Array{T4,1},1}, ub::Array{Array{T4,1},1}, d::Array{T4,1}, N::Array{T5,2}; L::Array{T6,1}=[1], Ld::Array{T6,1}=[1]) where {T0<:SAMIN, T1<:FitMC1d, T2<:FitProcedure, T3<:PlaneWave, T4<:Float64, T5<:ComplexF64, T6<:Int64}
    checkInput(xinit, layers)
    sol = optimize(x->fitMSE(x, arrange, specType, beam, Xexp, σ, layers, xinit, d, N, L, Ld), flattenArrays(lb), flattenArrays(ub), flattenArrays(xinit), alg, options)
    xfinal = arrayArrays(sol.minimizer, xinit)
    multilayerParameters!(N, d, arrayArrays(sol.minimizer, xfinal), beam, layers)
    N, d = buildArraysFP1d(N, d, L, Ld)
    X = computeTransferMatrix(specType, d, N, beam)
    return (sol=sol, X=X, xfinal=xfinal, solmin=sol.minimum)
end

# SAMIN, MC1dAlpha
function runFittingProcedure(alg::T0, arrange::T1, specType::T2, xinit::Array{Array{T4,1},1}, beam::T3, Xexp::Array{T4,1}, σ::Array{T4,1}, layers::Array, options, lb::Array{Array{T4,1},1}, ub::Array{Array{T4,1},1}, d::Array{T4,1}, N::Array{T5,2}; L::Array{T6,1}=[1], Ld::Array{T6,1}=[1]) where {T0<:SAMIN, T1<:FitMC1dAlpha, T2<:FitProcedure, T3<:PlaneWave, T4<:Float64, T5<:ComplexF64, T6<:Int64}
    checkInputAlpha(xinit, layers)
    sol = optimize(x->fitMSE(x, arrange, specType, beam, Xexp, σ, layers, xinit, d, N, L, Ld), flattenArrays(lb), flattenArrays(ub), flattenArrays(xinit), alg, options)
    xfinal = arrayArrays(sol.minimizer[1:end-1], xinit[1:end-1])
    push!(xfinal, [sol.minimizer[end]])
    multilayerParameters!(N, d, arrayArrays(sol.minimizer[1:end-1], xfinal[1:end-1]), beam, layers)
    N, d = buildArraysFP1d(N, d, L, Ld)
    monotonicLin!(d, sol.minimizer[end])
    X = computeTransferMatrix(specType, d, N, beam)
    return (sol=sol, X=X, xfinal=xfinal, solmin=sol.minimum)
end

# NelderMead, Generic
function runFittingProcedure(alg::T0, arrange::T1, specType::T2, xinit::Array{Array{T4,1},1}, beam::T3, Xexp::Array{T4,1}, σ::Array{T4,1}, layers::Array, options, lb::Array{Array{T4,1},1}, ub::Array{Array{T4,1},1}, d::Array{T4,1}, N::Array{T5,2}; kwargs...) where {T0<:NelderMead, T1<:FitGeneric, T2<:FitProcedure, T3<:PlaneWave, T4<:Float64, T5<:ComplexF64}
    checkInput(xinit, layers)
    sol = optimize(x->fitMSE(x, arrange, specType, beam, Xexp, σ, layers, xinit, d, N), flattenArrays(xinit), alg, options)
    xfinal = arrayArrays(sol.minimizer, xinit)
    multilayerParameters!(N, d, xfinal, beam, layers)
    X = computeTransferMatrix(specType, d, N, beam)
    return (sol=sol, X=X, xfinal=xfinal, solmin=sol.minimum)
end

# NelderMead, DBR
function runFittingProcedure(alg::T0, arrange::T1, specType::T2, xinit::Array{Array{T4,1},1}, beam::T3, Xexp::Array{T4,1}, σ::Array{T4,1}, layers::Array, options, lb::Array{Array{T4,1},1}, ub::Array{Array{T4,1},1}, d::Array{T4,1}, N::Array{T5,2}; L::Array{T6,1}=[1], kwargs...) where {T0<:NelderMead, T1<:FitDBR, T2<:FitProcedure, T3<:PlaneWave, T4<:Float64, T5<:ComplexF64, T6<:Int64}
    checkInput(xinit, layers)
    sol = optimize(x->fitMSE(x, arrange, specType, beam, Xexp, σ, layers, xinit, d, N, L[1]), flattenArrays(xinit), alg, options)
    xfinal = arrayArrays(sol.minimizer, xinit)
    multilayerParameters!(N, d, arrayArrays(sol.minimizer, xfinal), beam, layers)
    N, d = buildArraysBragg(N, d, L[1])
    X = computeTransferMatrix(specType, d, N, beam)
    return (sol=sol, X=X, xfinal=xfinal, solmin=sol.minimum)
end

# NelderMead, DBRAlpha
function runFittingProcedure(alg::T0, arrange::T1, specType::T2, xinit::Array{Array{T4,1},1}, beam::T3, Xexp::Array{T4,1}, σ::Array{T4,1}, layers::Array, options, lb::Array{Array{T4,1},1}, ub::Array{Array{T4,1},1}, d::Array{T4,1}, N::Array{T5,2}; L::Array{T6,1}=[1, L::Array{T6,1}, Ld::Array{T5,1}, N::Array{T6,2}, d::Array{T4,1}], kwargs...) where {T0<:NelderMead, T1<:FitDBRAlpha, T2<:FitProcedure, T3<:PlaneWave, T4<:Float64, T5<:ComplexF64, T6<:Int64}
    checkInputAlpha(xinit, layers)
    sol = optimize(x->fitMSE(x, arrange, specType, beam, Xexp, σ, layers, xinit, d, N, L[1]), flattenArrays(xinit), alg, options)
    xfinal = arrayArrays(sol.minimizer[1:end-1], xinit[1:end-1])
    push!(xfinal, [sol.minimizer[end]])
    multilayerParameters!(N, d, arrayArrays(sol.minimizer[1:end-1], xfinal[1:end-1]), beam, layers)
    N, d = buildArraysBragg(N, d, L[1])
    monotonicLin!(d, sol.minimizer[end])
    X = computeTransferMatrix(specType, d, N, beam)
    return (sol=sol, X=X, xfinal=xfinal, solmin=sol.minimum)
end

# NelderMead, MC1d
function runFittingProcedure(alg::T0, arrange::T1, specType::T2, xinit::Array{Array{T4,1},1}, beam::T3, Xexp::Array{T4,1}, σ::Array{T4,1}, layers::Array, options, lb::Array{Array{T4,1},1}, ub::Array{Array{T4,1},1}, d::Array{T4,1}, N::Array{T5,2}; L::Array{T6,1}, Ld::Array{T6,1}) where {T0<:NelderMead, T1<:FitMC1d, T2<:FitProcedure, T3<:PlaneWave, T4<:Float64, T5<:ComplexF64, T6<:Int64}
    checkInput(xinit, layers)
    sol = optimize(x->fitMSE(x, arrange, specType, beam, Xexp, σ, layers, xinit, d, N, L, Ld), flattenArrays(xinit), alg, options)
    xfinal = arrayArrays(sol.minimizer, xinit)
    multilayerParameters!(N, d, arrayArrays(sol.minimizer, xfinal), beam, layers)
    N, d = buildArraysFP1d(N, d, L, Ld)
    X = computeTransferMatrix(specType, d, N, beam)
    return (sol=sol, X=X, xfinal=xfinal, solmin=sol.minimum)
end

# NelderMead, MC1dAlpha
function runFittingProcedure(alg::T0, arrange::T1, specType::T2, xinit::Array{Array{T4,1},1}, beam::T3, Xexp::Array{T4,1}, σ::Array{T4,1}, layers::Array, options, lb::Array{Array{T4,1},1}, ub::Array{Array{T4,1},1}, d::Array{T4,1}, N::Array{T5,2}; L::Array{T6,1}, Ld::Array{T6,1}) where {T0<:NelderMead, T1<:FitMC1dAlpha, T2<:FitProcedure, T3<:PlaneWave, T4<:Float64, T5<:ComplexF64, T6<:Int64}
    checkInputAlpha(xinit, layers)
    sol = optimize(x->fitMSE(x, arrange, specType, beam, Xexp, σ, layers, xinit, d, N, L, Ld), flattenArrays(xinit), alg, options)
    xfinal = arrayArrays(sol.minimizer[1:end-1], xinit[1:end-1])
    push!(xfinal, [sol.minimizer[end]])
    multilayerParameters!(N, d, arrayArrays(sol.minimizer[1:end-1], xfinal[1:end-1]), beam, layers)
    N, d = buildArraysFP1d(N, d, L, Ld)
    monotonicLin!(d, sol.minimizer[end])
    X = computeTransferMatrix(specType, d, N, beam)
    return (sol=sol, X=X, xfinal=xfinal, solmin=sol.minimum)
end

"""

    Mean squared error (MSE) objective function used for the optimization process. (IU only)

        objfun = fitMSE(x, arrange, specType, beam, Xexp, σ, layers, d, N, L, Ld)

            x: array with parameters to optimze. The first one is the thickness
            arrange: Type of structure to fit, Generic(), DBR(), DBRAlpha(), MC(), MCAlpha()
            specType: type of spectrum to compute, Reflectance() or Transmittance()
            beam: PlaneWave structure
            Xexp: absolute experimental spectrum
            σ: std of the experimental spectrum
            layers: columnwise array of LayerTMMO1dIso and ModelFit with information about the layers
            d: thicknesses of the whole system to be updated
            N: index of refraction of the whole system to be updated
            L: period of each of the mirrors surrounding the defect(s)
            Ld: number of times the defect(s) layers are repeated. The 4th element of layers argument is considered to be the first defect layer, the fifth element of layers is the second defect and so on
            objfun: objective function value

"""
function fitMSE end

# Generic
function fitMSE(x::Array{T0,1}, arrange::T1, specType::T2, beam::T3, Xexp::Array{T0,1}, σ::Array{T0,1}, layers::Array, xinit::Array{Array{T0,1},1}, d::Array{T0,1}, N::Array{T4,2}) where {T0<:Float64, T1<:FitGeneric, T2<:FitProcedure, T3<:PlaneWave, T4<:ComplexF64}
    multilayerParameters!(N, d, arrayArrays(x, xinit), beam, layers)
    X = computeTransferMatrix(specType, d, N, beam)
    return meanSquaredError(vec(X), Xexp; σ=σ)
end

# DBR
function fitMSE(x::Array{T0,1}, arrange::T1, specType::T2, beam::T3, Xexp::Array{T0,1}, σ::Array{T0,1}, layers::Array, xinit::Array{Array{T0,1},1}, d::Array{T0,1}, N::Array{T4,2}, L::T5=1) where {T0<:Float64, T1<:FitDBR, T2<:FitProcedure, T3<:PlaneWave, T4<:ComplexF64, T5<:Int64}
    multilayerParameters!(N, d, arrayArrays(x, xinit), beam, layers)
    N, d = buildArraysBragg(N, d, L)
    X = computeTransferMatrix(specType, d, N, beam)
    return meanSquaredError(vec(X), Xexp; σ=σ)
end

# DBRAlpha
function fitMSE(x::Array{T0,1}, arrange::T1, specType::T2, beam::T3, Xexp::Array{T0,1}, σ::Array{T0,1}, layers::Array, xinit::Array{Array{T0,1},1}, d::Array{T0,1}, N::Array{T4,2}, L::T5=1) where {T0<:Float64, T1<:FitDBRAlpha, T2<:FitProcedure, T3<:PlaneWave, T4<:ComplexF64, T5<:Int64}
    multilayerParameters!(N, d, arrayArrays(x[1:end-1], xinit[1:end-1]), beam, layers)
    N, d = buildArraysBragg(N, d, L)
    monotonicLin!(d, x[end])
    X = computeTransferMatrix(specType, d, N, beam)
    return meanSquaredError(vec(X), Xexp; σ=σ)
end

# MC1d
function fitMSE(x::Array{T0,1}, arrange::T1, specType::T2, beam::T3, Xexp::Array{T0,1}, σ::Array{T0,1}, layers::Array, xinit::Array{Array{T0,1},1}, d::Array{T0,1}, N::Array{T4,2}, L::Array{T5,1}=[1], Ld::Array{T5,1}=[1]) where {T0<:Float64, T1<:FitMC1d, T2<:FitProcedure, T3<:PlaneWave, T4<:ComplexF64, T5<:Int64}
    multilayerParameters!(N, d, arrayArrays(x, xinit), beam, layers)
    N, d = buildArraysFP1d(N, d, L, Ld)
    X = computeTransferMatrix(specType, d, N, beam)
    return meanSquaredError(vec(X), Xexp; σ=σ)
end

# MC1dAlpha
function fitMSE(x::Array{T0,1}, arrange::T1, specType::T2, beam::T3, Xexp::Array{T0,1}, σ::Array{T0,1}, layers::Array, xinit::Array{Array{T0,1},1}, d::Array{T0,1}, N::Array{T4,2}, L::Array{T5,1}=[1], Ld::Array{T5,1}=[1]) where {T0<:Float64, T1<:FitMC1dAlpha, T2<:FitProcedure, T3<:PlaneWave, T4<:ComplexF64, T5<:Int64}
    multilayerParameters!(N, d, arrayArrays(x[1:end-1], xinit[1:end-1]), beam, layers)
    N, d = buildArraysFP1d(N, d, L, Ld)
    monotonicLin!(d, x[end])
    X = computeTransferMatrix(specType, d, N, beam)
    return meanSquaredError(vec(X), Xexp; σ=σ)
end

"""

    Compute the transfer matrix and return the specified spectra.

        X = computeTransferMatrix(specType, d, N, beam)

            specType: type of spectrum to compute, Reflectance() or Transmittance()
            d: thicknesses of the whole system to be updated
            N: index of refraction of the whole system to be updated
            beam: PlaneWave structure
            X: polarisation averaged modelled spectrum

"""
function computeTransferMatrix end

# Transmittance
function computeTransferMatrix(specType::T0, d::Array{T1,1}, N::Array{T2,2}, beam::T3) where {T0<:FitTransmittance, T1<:Float64, T2<:ComplexF64, T3<:PlaneWave}
    spectra = TransferMatrix(N, d, beam)[1]
    return averagePolarisation(beam.p, spectra.Tp, spectra.Ts)
end

# Reflectance
function computeTransferMatrix(specType::T0, d::Array{T1,1}, N::Array{T2,2}, beam::T3) where {T0<:FitReflectance, T1<:Float64, T2<:ComplexF64, T3<:PlaneWave}
    spectra = TransferMatrix(N, d, beam)[1]
    return averagePolarisation(beam.p, spectra.Rp, spectra.Rs)
end

end # module
