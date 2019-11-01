module FTFSpectrum

using Optim
using Statistics
using ..CommonStructures: ModelFit, BoundariesFit, PlaneWave, LayerTMMO, FitProcedure
using ..CommonStructures: getBeamParameters
using ..Utils: build_interpolation, averagePolarisation, flattenArrays, arrayArrays, oscillatorsInput
using ..TransferMatrixMethod: TransferMatrix
using ..RI

export SpaceSolutionEMA,
       TheoreticalSpectrum,
       NormalizeReflectance,
       FitTMMOptics,
       Reflectance,
       Transmittance,
       Ellipsometry,
       NoAlpha,
       UseAlpha,
       MeanAbs,
       SumAbs,
       SumMeanAbs

## Definition of subtypes
mutable struct Ellipsometry{T1,T2} <: FitProcedure where {T1<:Real, T2<:Real}
    X::Array{T1,2}
    Xexp::Array{T2,2}
    function Ellipsometry(
        X::Array{T1,2};
        Xexp::Array{T2,2}=zeros(size(X)),
    ) where {T1<:Real, T2<:Real}
      return new{T1,T2}(X, Xexp)
    end
end
Ellipsometry() = Ellipsometry([]; Xexp=zeros(1,2))

struct Reflectance{T1} <: FitProcedure where {T1<:Real}
    Xexp::Array{T1,1}
    function Reflectance(Xexp::Array{T1,1}) where {T1<:Real}
      return new{T1}(Xexp)
    end
end
Reflectance() = Reflectance([0.0])

struct Transmittance{T1} <: FitProcedure where {T1<:Real}
    Xexp::Array{T1,1}
    function Transmittance(Xexp::Array{T1,1}) where {T1<:Real}
      return new{T1}(Xexp)
    end
end
Transmittance() = Transmittance([0.0])

struct MeanAbsObjFun <: FitProcedure end
MeanAbs() = MeanAbsObjFun()
struct SumAbsObjFun <: FitProcedure end
SumAbs() = SumAbsObjFun()
struct SumMeanAbsObjFun <: FitProcedure end
SumMeanAbs() = SumMeanAbsObjFun()

# Non-exported
struct DoNotAlpha <: FitProcedure end
NoAlpha() = DoNotAlpha()
struct DoAlpha <: FitProcedure end
UseAlpha() = DoAlpha()


"""

    Normalize the experimental reflectance dividing by the reference reflectance and
    multiplying by the theoretical one. The experimental and the reference reflectances
    must have the same scale.

        R = NormalizeReflectance(λ, Rexp, Rthe, Rref)

            λ: wavelength of interest [nm]
            Rexp: 2-columns array with wavelength in first one and experimental
                reflectance in the second
            Rthe: 2-columns array with wavelength in first one and theoretical
                reflectance in the second
            Rref: 2-columns array with wavelength in first one and reference reflectance
                in the second

    The output has the same scale as the theoretical reflectance with the same length as λ.

"""
function NormalizeReflectance(
    λ::AbstractArray{T0,1}, Rexp::Array{T1,2}, Rthe::Array{T1,2}, Rref::Array{T1,2},
) where {T0<:Real, T1<:Float64}
    # Interpolate for the λ range
    spl_exp = build_interpolation(Rexp)
    spl_ref = build_interpolation(Rref)
    spl_the = build_interpolation(Rthe)
    exp2ref = spl_exp.(λ) ./ spl_ref.(λ)
    # Find out if the sample is overshooting the reference
    (maximum(exp2ref) < 1.0) || return (exp2ref./maximum(exp2ref).*spl_the.(λ))
    return exp2ref.*spl_the.(λ)
end

"""

    Returns the polarisation averaged theoretical spectrum calculated with the transfer
    matrix method. It consiers light hitting the first medium (incidentRI).

        X = TheoreticalSpectrum(specType, beam, incidentRI, emergentRI)

            specType: type of spectrum to compute, Reflectance() or Transmittance().
            beam: a PlaneWave structure.
            incidentRI: Array containing the index of refraction of the incident material
                for several wavelengths. For instance, RIdb.air(beam.λ).
            emergentRI: Array containing the index of refraction of the emergent material
                for several wavelengths. For instance, RIdb.silicon(beam.λ).
            X: averaged spectrum = beam.p*Xp + (1.0 - beam.p)*Xs

"""
function TheoreticalSpectrum(
    specType::T0, beam::T1, N1::Array{T2,1}, N2::Array{T2,1},
) where {T0<:FitProcedure, T1<:PlaneWave, T2<:ComplexF64}
    beam_, _ = getBeamParameters(beam)
    X = computeTransferMatrix(specType, vec([0. 100.0 0.]), [N1 N2 N2], beam_)
    return X
end

"""

    Computes the solution space for a window range of two parameters, the optical thickness
    and porosity. Only works with these two parameters and for a single layer between
    incident and emergent media.

        sol = SpaceSolutionEMA(specType, b, beam, Xexp, layers; oftype=MeanAbs())

            specType: Type of spectrum to fit, with accepted values:
                Reflectance(Xexp), Transmittance(Xexp) or Ellipsometry(Xexp).
                For the Reflectance() and Transmittance() you must pass the spectrum array
                Xexp (between 0 and 1) to be fitted. For the Ellipsometry() you have to
                pass the Ψ (first column) and Δ (second column) spectra in degrees.
            b: lower and upper boundaries for the optical thickness and porosity.
               (e.g. [odl, odu, pl pu; Nod, Np]), where Nod and Np indicates the number of
               points to search in the solution space as optional parameters.
            beam: structure from PlaneWave
            layers: columnwise array of three LayerTMMO and ModelFit with information
                about the layers. For SpaceSolutionEMA: length(vec(layers)) = 3. Notice that
                the type of the first and third layers is LayerTMMO while the type of the
                second layers is ModelFit.
                oftype: Metric for the objective function to optimise.
                    Can take values: MeanAbs() (default), SumAbs() and SumMeanAbs().
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
function SpaceSolutionEMA(
    specType::T0, b::T1, beam::T2, layers::Array;
    oftype::T4=MeanAbs(),
) where{T0<:FitProcedure, T1<:BoundariesFit, T2<:PlaneWave, T3<:Float64, T4<:FitProcedure}
    # Generate grids
    _d = LinRange(b.odlo, b.odup, b.Nod)
    _p = LinRange(b.plo, b.pup, b.Np)
    d = similar(_d)
    sol = zeros(Float64, b.Nod, b.Np)
    beam_, _ = getBeamParameters(beam)
    # If it is ellipsometry, then fit the norm of ρ
    if isa(specType, Ellipsometry)
        ρexp = @. tan(deg2rad(specType.X[:,1]))*exp(im*deg2rad(specType.X[:,2]))
        specType.Xexp = hcat(real.(ρexp), imag.(ρexp))
        # specType.Xexp = @. real(ρexp*conj(ρexp))
    end
    @inbounds for j in eachindex(_p), i in eachindex(_d)
        neff = refractiveIndexMR(layers[2], [_p[j]], beam_.λ)
        d[i] = _d[i] / mean(real.(neff))
        X = computeTransferMatrix(specType, vec([0. d[i] 0.]), [layers[1].n neff layers[3].n], beam_)
        sol[i, j] = objectiveFunction(oftype, X, specType.Xexp)
    end
    smin = findmin(sol)
    X = computeTransferMatrix(
            specType,
            vec([0. d[smin[2][1]] 0.]),
            [layers[1].n refractiveIndexMR(layers[2], [_p[smin[2][2]]], beam_.λ) layers[3].n],
            beam_,
    )
    sol_ = (
        solSpace = sol,
        objfunMin = smin[1],
        optOD = _d[smin[2][1]],
        optPorosity = _p[smin[2][2]],
        optThickness = d[smin[2][1]],
        spectrumFit = X,
        spectrumExp = specType.Xexp,
        od = _d,
        p = _p,
        Beam = (λ=beam_.λ, θ=beam_.θ, pol=beam_.p),
    )
    return sol_
end

"""

    Optimize the parameters of the thin films in a multilayer stack using the selected models and the experimental spectrum.

        sol = FitTMMOptics(
            specType, xinit, beam, Rexp, layers;
            order=1:length(layers), options=Optim.Options(), alg=SAMIN(),
            lb=0.5.*xinit, ub=1.5.*xinit, σ=ones(size(Xexp)), alpha=false,
            oftype=MeanAbs(),
        )

            specType: Type of spectrum to fit, with accepted values:
                Reflectance(Xexp), Transmittance(Xexp) or Ellipsometry(Xexp).
                For the Reflectance() and Transmittance() you must pass the spectrum array
                Xexp (0<=Xexp<=1) to be fitted. For the Ellipsometry() you have to
                pass the Ψ (first column) and Δ (second column) spectra in degrees, i.e.,
                Xexp = [Ψ Δ].
            xinit: Array with the initial parameters for optimization. For instance, to
                optimise two layers in the system you need to wrap the seeds for the
                algorithm as follow: xinit = [seed1, seed2], even if it is only one
                single layer to fit, xinit = [seed1]. The first parameter of each seed
                must be always the thickness, the rest are the parameters for the selected
                model. The number of seeds must match that of the number of ModelFit layers.
            beam: Structure from PlaneWave.
            layers: Columnwise array of LayerTMMO and ModelFit with information about the
                layers.
                order: Set the order of the layers (build the system) with the layers.
                    This parameter is directly bound to the (previously defined) layers
                    parameter. Order basically contains the indices of each layer inside
                    the parameter layers. Notice that the first and last indices in
                    order cannot be associated to ModelFit layers, as they are the
                    outter media.
                options: Optional Optim.Options structure.
                alg: Optional algorithm method selected, by default takes SAMIN().
                    You can pass the options inside as well, for instance, SAMIN(rt=0.1).
                    Right now, LBFGS(), NelderMead() and SAMIN() are supported from
                    Optim.jl. If you select SAMIN() you need to input lb and ub, otherwise
                    will be set as 0.5 and 1.5 times the xinit argument, respectively.
                    The LBFGS algorithm requires the jacobian for the optimisation but it
                    can be omitted and the algorithm will calculate the differentials
                    numerically.
                lb: Lower bounds for the optimisation variables, by default lb=0.5.*xinit
                ub: Upper bounds for the optimisation variables, by default ub=1.5.*xinit
                σ: Array with the standard deviation of the spectrum for each wavelength,
                   by default is ones(size(Xexp)).
                alpha: Tells the procedure to either use (alpha=true) a linear decrease in
                    the thicknesses of the modelling layers or not (alpha=false). The
                    default is false, do not use alpha.
                oftype: Metric for the objective function to optimise.
                    Can take values: MeanAbs() (default), SumAbs() and SumMeanAbs().

            sol: FitTMMOptics structure with results, with fields
                OptimSolverInfo: Status info from the solver
                spectrumFit: Spectrum obtained with the model and the optimal parameters
                spectrumExp: Input experimental spectrum
                optParams: Array of arrays with optimal parameters
                objfunMin: Optimum value of objective function MSE
                beam: Structure of the light used
                fitOptions: NamedTuple with information of the fitting procedure

"""
function FitTMMOptics(
    specType::T0,
    xinit,
    beam::T1,
    layers::Array;
    order::AbstractArray{T3,N3}=1:length(layers),
    options=Optim.Options(),
    alg=SAMIN(),
    lb=0.5.*xinit,
    ub=1.5.*xinit,
    σ::Array{T4}=ones(size(specType.Xexp)),
    alpha::T5=false,
    oftype::T7=MeanAbs(),
),
) where{T0<:FitProcedure, T1<:PlaneWave, T3<:Int64, N3, T4<:Real, T5<:Bool, T6<:Symbol, T7<:FitProcedure}
    isa(layers[order[1]], LayerTMMO) || throw("The first layer of the system cannot be ModelFit type. Check the order parameter.")
    isa(layers[order[end]], LayerTMMO) || throw("The last layer of the system cannot be ModelFit type. Check the order parameter.")
    if sum(isa.(layers, ModelFit)) < 1
        error("There should be at least one layer to modelling inside the system (ModelFit type).")
    end
    # Get the parameters for beam
    beam_, _ = getBeamParameters(beam)
    # Warm-up
    N = Array{ComplexF64,2}(undef, length(beam_.λ), length(vec(order)))
    d = Array{Float64,1}(undef, length(vec(order)))
    # Check alpha input
    α = alpha ? UseAlpha() : NoAlpha()
    # Build the indices vector for the seeds
    idxSeed = _indicesSeed(layers, vec(order))
    # If it is ellipsometry, then fit norm of ρ
    if isa(specType, Ellipsometry)
        ρexp = @. tan(deg2rad(specType.X[:,1]))*exp(im*deg2rad(specType.X[:,2]))
        specType.Xexp = hcat(real.(ρexp), imag.(ρexp))
        # specType.Xexp = @. real(ρexp*conj(ρexp))
    end
    # Run the solver depending on the algorithm
    if isa(alg, Symbol)
        solution = runFittingProcedureBBO(
            alg,
            specType,
            float.(xinit),
            beam_,
            float.(σ),
            layers,
            vec(collect(order)),
            optionsbbo,
            float.(lb),
            float.(ub),
            d,
            N,
            idxSeed,
            α,
            oftype,
        )
    else
        solution = runFittingProcedure(
            alg,
            specType,
            float.(xinit),
            beam_,
            float.(σ),
            layers,
            vec(collect(order)),
            options,
            float.(lb),
            float.(ub),
            d,
            N,
            idxSeed,
            α,
            oftype,
        )
    end
    solution_ = (
        spectrumFit = solution.X,
        spectrumExp = specType.Xexp,
        optParams = solution.xfinal,
        objfunMin = solution.solmin,
        Beam = (λ=beam_.λ, θ=beam_.θ, pol=beam_.p),
        layers = layers,
        fitOptions = (
            fitType = specType,
            optimSolverInfo = solution.sol,
        ),
    )
    return solution_
end

## Run the fitting procedure for Optim's solvers depending on input.
function runFittingProcedure end

# SAMIN
function runFittingProcedure(
    alg::T1,
    specType::T2,
    xinit::Array{Array{T3,1},1},
    beam::T4,
    σ::Array{T3},
    layers::Array,
    order::Array{T5,1},
    options,
    lb::Array{Array{T3,1},1},
    ub::Array{Array{T3,1},1},
    d::Array{T3,1},
    N::Array{T6,2},
    idxSeed::Array{T5,1},
    alpha::T7,
    oftype::T8,
) where {T1<:SAMIN, T2<:FitProcedure, T3<:Float64, T4<:PlaneWave, T5<:Int64, T6<:ComplexF64, T7<:DoNotAlpha, T8<:FitProcedure}
    _checkInput(xinit, layers)
    sol = optimize(
            x->fitOF(x, specType, beam, σ, layers, order, xinit, d, N, idxSeed, alpha, oftype),
            flattenArrays(lb),
            flattenArrays(ub),
            flattenArrays(xinit),
            alg,
            options,
    )
    xfinal = arrayArrays(sol.minimizer, xinit)
    multilayerParameters!(N, d, xfinal, beam, layers, order, idxSeed)
    X = computeTransferMatrix(specType, d, N, beam)
    sol_ = (sol=sol, X=X, xfinal=xfinal, solmin=sol.minimum)
    return sol_
end

# SAMIN, alpha
function runFittingProcedure(
    alg::T1,
    specType::T2,
    xinit::Array{Array{T3,1},1},
    beam::T4,
    σ::Array{T3},
    layers::Array,
    order::Array{T5,1},
    options,
    lb::Array{Array{T3,1},1},
    ub::Array{Array{T3,1},1},
    d::Array{T3,1},
    N::Array{T6,2},
    idxSeed::Array{T5,1},
    alpha::T7,
    oftype::T8,
) where {T1<:SAMIN, T2<:FitProcedure, T3<:Float64, T4<:PlaneWave, T5<:Int64, T6<:ComplexF64, T7<:DoAlpha, T8<:FitProcedure}
    _checkInputAlpha(xinit, layers)
    sol = optimize(
            x->fitOF(x, specType, beam, σ, layers, order, xinit, d, N, idxSeed, alpha, oftype),
            flattenArrays(lb),
            flattenArrays(ub),
            flattenArrays(xinit),
            alg,
            options,
    )
    xfinal = arrayArrays(sol.minimizer[1:end-1], xinit[1:end-1])
    push!(xfinal, [sol.minimizer[end]])
    multilayerParameters!(N, d, arrayArrays(sol.minimizer[1:end-1], xfinal[1:end-1]), beam, layers, order, idxSeed)
    monotonicLin!(d, sol.minimizer[end])
    X = computeTransferMatrix(specType, d, N, beam)
    sol_ = (sol=sol, X=X, xfinal=xfinal, solmin=sol.minimum)
    return sol_
end

# NelderMead
function runFittingProcedure(
    alg::T1,
    specType::T2,
    xinit::Array{Array{T3,1},1},
    beam::T4,
    σ::Array{T3},
    layers::Array,
    order::Array{T5,1},
    options,
    lb::Array{Array{T3,1},1},
    ub::Array{Array{T3,1},1},
    d::Array{T3,1},
    N::Array{T6,2},
    idxSeed::Array{T5,1},
    alpha::T7,
    oftype::T8,
) where {T1<:NelderMead, T2<:FitProcedure, T3<:Float64, T4<:PlaneWave, T5<:Int64, T6<:ComplexF64, T7<:DoNotAlpha, T8<:FitProcedure}
    _checkInput(xinit, layers)
    sol = optimize(
            x->fitOF(x, specType, beam, σ, layers, order, xinit, d, N, idxSeed, alpha, oftype),
            flattenArrays(xinit),
            alg,
            options,
    )
    xfinal = arrayArrays(sol.minimizer, xinit)
    multilayerParameters!(N, d, xfinal, beam, layers, order, idxSeed)
    X = computeTransferMatrix(specType, d, N, beam)
    sol_ = (sol=sol, X=X, xfinal=xfinal, solmin=sol.minimum)
    return sol_
end

# NelderMead, alpha
function runFittingProcedure(
    alg::T1,
    specType::T2,
    xinit::Array{Array{T3,1},1},
    beam::T4,
    σ::Array{T3},
    layers::Array,
    order::Array{T5,1},
    options,
    lb::Array{Array{T3,1},1},
    ub::Array{Array{T3,1},1},
    d::Array{T3,1},
    N::Array{T6,2},
    idxSeed::Array{T5,1},
    alpha::T7,
    oftype::T8,
) where {T1<:NelderMead, T2<:FitProcedure, T3<:Float64, T4<:PlaneWave, T5<:Int64, T6<:ComplexF64, T7<:DoAlpha, T8<:FitProcedure}
    _checkInputAlpha(xinit, layers)
    sol = optimize(
            x->fitOF(x, specType, beam, σ, layers, order, xinit, d, N, idxSeed, alpha, oftype),
            flattenArrays(xinit),
            alg,
            options,
    )
    xfinal = arrayArrays(sol.minimizer[1:end-1], xinit[1:end-1])
    push!(xfinal, [sol.minimizer[end]])
    multilayerParameters!(N, d, arrayArrays(sol.minimizer[1:end-1], xfinal[1:end-1]), beam, layers, order, idxSeed)
    monotonicLin!(d, sol.minimizer[end])
    X = computeTransferMatrix(specType, d, N, beam)
    sol_ = (sol=sol, X=X, xfinal=xfinal, solmin=sol.minimum)
    return sol_
end

# LBFGS
function runFittingProcedure(
    alg::T1,
    specType::T2,
    xinit::Array{Array{T3,1},1},
    beam::T4,
    σ::Array{T3},
    layers::Array,
    order::Array{T5,1},
    options,
    lb::Array{Array{T3,1},1},
    ub::Array{Array{T3,1},1},
    d::Array{T3,1},
    N::Array{T6,2},
    idxSeed::Array{T5,1},
    alpha::T7,
    oftype::T8,
) where {T1<:LBFGS, T2<:FitProcedure, T3<:Float64, T4<:PlaneWave, T5<:Int64, T6<:ComplexF64, T7<:DoNotAlpha, T8<:FitProcedure}
    _checkInput(xinit, layers)
    sol = optimize(
            x->fitOF(x, specType, beam, σ, layers, order, xinit, d, N, idxSeed, alpha, oftype),
            flattenArrays(xinit),
            alg,
            options,
    )
    xfinal = arrayArrays(sol.minimizer, xinit)
    multilayerParameters!(N, d, xfinal, beam, layers, order, idxSeed)
    X = computeTransferMatrix(specType, d, N, beam)
    sol_ = (sol=sol, X=X, xfinal=xfinal, solmin=sol.minimum)
    return sol_
end

# LBFGS, alpha
function runFittingProcedure(
    alg::T1,
    specType::T2,
    xinit::Array{Array{T3,1},1},
    beam::T4,
    σ::Array{T3},
    layers::Array,
    order::Array{T5,1},
    options,
    lb::Array{Array{T3,1},1},
    ub::Array{Array{T3,1},1},
    d::Array{T3,1},
    N::Array{T6,2},
    idxSeed::Array{T5,1},
    alpha::T7,
    oftype::T8,
) where {T1<:LBFGS, T2<:FitProcedure, T3<:Float64, T4<:PlaneWave, T5<:Int64, T6<:ComplexF64, T7<:DoAlpha, T8<:FitProcedure}
    _checkInputAlpha(xinit, layers)
    sol = optimize(
            x->fitOF(x, specType, beam, σ, layers, order, xinit, d, N, idxSeed, alpha, oftype),
            flattenArrays(xinit),
            alg,
            options,
    )
    xfinal = arrayArrays(sol.minimizer[1:end-1], xinit[1:end-1])
    push!(xfinal, [sol.minimizer[end]])
    multilayerParameters!(N, d, arrayArrays(sol.minimizer[1:end-1], xfinal[1:end-1]), beam, layers, order, idxSeed)
    monotonicLin!(d, sol.minimizer[end])
    X = computeTransferMatrix(specType, d, N, beam)
    sol_ = (sol=sol, X=X, xfinal=xfinal, solmin=sol.minimum)
    return sol_
end

# NewtonTrustRegion
function runFittingProcedure(
    alg::T1,
    specType::T2,
    xinit::Array{Array{T3,1},1},
    beam::T4,
    σ::Array{T3},
    layers::Array,
    order::Array{T5,1},
    options,
    lb::Array{Array{T3,1},1},
    ub::Array{Array{T3,1},1},
    d::Array{T3,1},
    N::Array{T6,2},
    idxSeed::Array{T5,1},
    alpha::T7,
    oftype::T8,
) where {T1<:NewtonTrustRegion, T2<:FitProcedure, T3<:Float64, T4<:PlaneWave, T5<:Int64, T6<:ComplexF64, T7<:DoNotAlpha, T8<:FitProcedure}
    _checkInput(xinit, layers)
    sol = optimize(
            x->fitOF(x, specType, beam, σ, layers, order, xinit, d, N, idxSeed, alpha, oftype),
            flattenArrays(xinit),
            alg,
            options,
    )
    xfinal = arrayArrays(sol.minimizer, xinit)
    multilayerParameters!(N, d, xfinal, beam, layers, order, idxSeed)
    X = computeTransferMatrix(specType, d, N, beam)
    sol_ = (sol=sol, X=X, xfinal=xfinal, solmin=sol.minimum)
    return sol_
end

# NewtonTrustRegion, alpha
function runFittingProcedure(
    alg::T1,
    specType::T2,
    xinit::Array{Array{T3,1},1},
    beam::T4,
    σ::Array{T3},
    layers::Array,
    order::Array{T5,1},
    options,
    lb::Array{Array{T3,1},1},
    ub::Array{Array{T3,1},1},
    d::Array{T3,1},
    N::Array{T6,2},
    idxSeed::Array{T5,1},
    alpha::T7,
    oftype::T8,
) where {T1<:NewtonTrustRegion, T2<:FitProcedure, T3<:Float64, T4<:PlaneWave, T5<:Int64, T6<:ComplexF64, T7<:DoAlpha, T8<:FitProcedure}
    _checkInputAlpha(xinit, layers)
    sol = optimize(
            x->fitOF(x, specType, beam, σ, layers, order, xinit, d, N, idxSeed, alpha, oftype),
            flattenArrays(xinit),
            alg,
            options,
    )
    xfinal = arrayArrays(sol.minimizer[1:end-1], xinit[1:end-1])
    push!(xfinal, [sol.minimizer[end]])
    multilayerParameters!(N, d, arrayArrays(sol.minimizer[1:end-1], xfinal[1:end-1]), beam, layers, order, idxSeed)
    monotonicLin!(d, sol.minimizer[end])
    X = computeTransferMatrix(specType, d, N, beam)
    sol_ = (sol=sol, X=X, xfinal=xfinal, solmin=sol.minimum)
    return sol_
end

# Check the input sizes and the seed with alpha.
function _checkInput(xinit, layers)
    length(xinit) == length(findall(isa.(layers, ModelFit))) || throw("There should be seeds for all active layers (ModelFit).")
    return nothing
end

function _checkInputAlpha(xinit, layers)
    length(xinit) == length(findall(isa.(layers, ModelFit)))+1 || throw("There should be seeds for all active layers (ModelFit) plus the last element with the alpha.")
    return nothing
end

## Fit objective function used for the optimization process.
function fitOF end

# No alpha
function fitOF(
    x::Array{T1,1},
    specType::T2,
    beam::T3,
    σ::Array{T1},
    layers::Array,
    order::Array{T4,1},
    xinit::Array{Array{T1,1},1},
    d::Array{T1,1},
    N::Array{T5,2},
    idxSeed::Array{T4,1},
    alpha::T6,
    oftype::T7,
) where {T1<:Float64, T2<:FitProcedure, T3<:PlaneWave, T4<:Int64, T5<:ComplexF64, T6<:DoNotAlpha, T7<:FitProcedure}
    multilayerParameters!(N, d, arrayArrays(x, xinit), beam, layers, order, idxSeed)
    X = computeTransferMatrix(specType, d, N, beam)
    mse = objectiveFunction(oftype, X, specType.Xexp; σ=σ)
    return mse
end

# Alpha
function fitOF(
    x::Array{T1,1},
    specType::T2,
    beam::T3,
    σ::Array{T1},
    layers::Array,
    order::Array{T4,1},
    xinit::Array{Array{T1,1},1},
    d::Array{T1,1},
    N::Array{T5,2},
    idxSeed::Array{T4,1},
    alpha::T6,
    oftype::T7,
) where {T1<:Float64, T2<:FitProcedure, T3<:PlaneWave, T4<:Int64, T5<:ComplexF64, T6<:DoAlpha, T7<:FitProcedure}
    multilayerParameters!(N, d, arrayArrays(x[1:end-1], xinit[1:end-1]), beam, layers, order, idxSeed)
    monotonicLin!(d, x[end])
    X = computeTransferMatrix(specType, d, N, beam)
    mse = objectiveFunction(oftype, X, specType.Xexp; σ=σ)
    return mse
end

## Compute the transfer matrix and return the specified spectra.
function computeTransferMatrix end

# Transmittance
function computeTransferMatrix(
    specType::T0, d::Array{T1,1}, N::Array{T2,2}, beam::T3,
) where {T0<:Transmittance, T1<:Float64, T2<:ComplexF64, T3<:PlaneWave}
    spectra = TransferMatrix(N, d, beam)[1]
    avg = averagePolarisation(beam.p, spectra.Tp, spectra.Ts)
    return avg
end

# Reflectance
function computeTransferMatrix(
    specType::T0, d::Array{T1,1}, N::Array{T2,2}, beam::T3,
) where {T0<:Reflectance, T1<:Float64, T2<:ComplexF64, T3<:PlaneWave}
    spectra = TransferMatrix(N, d, beam)[1]
    avg = averagePolarisation(beam.p, spectra.Rp, spectra.Rs)
    return avg
end

# Ellipsometry
function computeTransferMatrix(
    specType::T0, d::Array{T1,1}, N::Array{T2,2}, beam::T3,
) where {T0<:Ellipsometry, T1<:Float64, T2<:ComplexF64, T3<:PlaneWave}
    spectra = TransferMatrix(N, d, beam)[1]
    ρ::Array{ComplexF64,1} = vec(spectra.ρp./spectra.ρs)
    return hcat(-real.(ρ), imag.(ρ))
    # ρnorm = @. real(ρ*conj(ρ))
    # return ρnorm
end

## Metrics (objective function) for the optimisation process.

# MeanAbs() = MeanAbsObjFun
function objectiveFunction(
    oftype::T0, X::Array{T1}, Xexp::Array{T1};
    σ::Array{T1}=ones(size(Xexp)),
) where {T0<:MeanAbsObjFun, T1<:Float64}
    _absq = abs2.(X .- Xexp)./σ
    _mse = mean(_absq, dims=1)
    objfun = sum(_mse)
    return objfun
end

# SumAbs() = SumAbsObjFun
function objectiveFunction(
    oftype::T0, X::Array{T1}, Xexp::Array{T1};
    σ::Array{T1}=ones(size(Xexp)),
) where {T0<:SumAbsObjFun, T1<:Float64}
    _absq = abs2.(X .- Xexp)./σ
    _mse = sum(_absq, dims=1)
    objfun = sum(_mse)
    return objfun
end

# SumMeanAbs() = SumMeanAbsObjFun
function objectiveFunction(
    oftype::T0, X::Array{T1}, Xexp::Array{T1};
    σ::Array{T1}=ones(size(Xexp)),
) where {T0<:SumMeanAbsObjFun, T1<:Float64}
    _absq = abs2.(X .- Xexp)./σ
    _mse = mean(_absq, dims=1)
    _sse = sum(_absq, dims=1)
    # objfun = sqrt(sum(_mse + _sse)/size(Xexp, 1))
    objfun = sum(_mse + _sse)/size(Xexp, 1)
    return objfun
end

## Returns the index of refraction for the selected model. RI is imported from RefractiveIndicesModels.jl.
function refractiveIndexMR(
    layer::T1, argin::AbstractArray{T2,1}, λ::Array{T2,1},
) where {T1<:ModelFit, T2<:Float64}
    if layer.model == :bruggeman
        _errorParametersEMA(layer)
        return RI.bruggeman(argin[1], [layer.N.ninc layer.N.nhost]; df=layer.df)
    elseif layer.model == :looyenga
        _errorParametersEMA(layer)
        return RI.looyenga([argin[1] 1.0-argin[1]], [layer.N.ninc layer.N.nhost]; df=layer.df)
    elseif layer.model == :monecke
        _errorParametersEMA(layer)
        return RI.monecke(argin[1], [layer.N.ninc layer.N.nhost]; dfm=layer.dfm, c=layer.c)
    elseif layer.model == :maxwell
        _errorParametersEMA(layer)
        return RI.maxwell([argin[1] 1.0-argin[1]], [layer.N.ninc layer.N.nhost]; df=layer.df)
    elseif layer.model == :gedf
        _errorParametersEMA(layer)
        return RI.gedf(argin, [layer.N.ninc layer.N.nhost])
    elseif layer.model == :gem
        _errorParametersEMA(layer)
        return RI.gem(argin, [layer.N.ninc layer.N.nhost])
    elseif layer.model == :lorentzlorenz
        _errorParametersEMA(layer)
        return RI.lorentzlorenz(argin[1], [layer.N.ninc layer.N.nhost])
    elseif layer.model == :sellmeier
        return RI.sellmeier(argin, λ./1e3) # Sellmeier eq takes λ in μm
    elseif layer.model == :cauchyurbach
    	return RI.cauchyurbach(argin, λ)
    elseif layer.model == :drudelorentz
        argin_ = [argin[1], argin[2:end]]
        x = oscillatorsInput(argin, 3)
        return RI.drudelorentz(x, 1240.0./λ)
    elseif layer.model == :tauclorentz
        argin_ = [argin[1:2], argin[3:end]]
        x = oscillatorsInput(argin_, 3)
    	return RI.tauclorentz(x, 1240.0./λ)
    elseif layer.model == :forouhibloomer
        argin_ = [argin[1:2], argin[3:end]]
        x = oscillatorsInput(argin_, 3)
        return RI.forouhibloomer(x, 1240.0./λ)
    elseif layer.model == :lorentzplasmon
        argin_ = [argin[1:3], argin[4:end]]
        x = oscillatorsInput(argin_, 3)
        return RI.lorentzplasmon(x, 1240.0./λ)
    elseif layer.model == :codylorentz
        argin_ = [argin[1:5], argin[6:end]]
        x = oscillatorsInput(argin_, 3)
        return RI.codylorentz(x, 1240.0./λ; cme=layer.cme)
    end
end

## Error if EMA models do not have the two columns indices of refraction.
function _errorParametersEMA(layer::ModelFit)
    !|(isempty(layer.N[1]), isempty(layer.N[2])) || throw("To use the $(layer.model) model you need to specify the two indices of refraction of the components in ModelFit.")
    return nothing
end

## Construct the multilayers parameters to optimise depending on their type (LayerTMMO, ModelFit).
function multilayerParameters!(
    N::Array{T1,2},
    d::Array{T2,1},
    x::Array{Array{T2,1},1},
    beam::T3,
    layers::Array,
    order::Array{T4,1},
    idxSeed::Array{T4,1},
) where{T1<:ComplexF64, T2<:Float64, T3<:PlaneWave, T4<:Int64}
    k::Int64 = 1
    for i in eachindex(order)
        idx = order[i]
        if isa(layers[idx], LayerTMMO)
            d[i] = layers[idx].d
            @views N[:,i] = layers[idx].n
        elseif isa(layers[idx], ModelFit)
            d[i] = x[idxSeed[k]][1] # the thickness is defined first
            @views N[:,i] = refractiveIndexMR(layers[idx], x[idxSeed[k]][2:end], beam.λ)
            k += 1
        end
    end
    return nothing
end

## Build the indices vector for the seeds
function _indicesSeed(layers, order)
    mflayers = findall(isa.(layers, ModelFit))
    idx = 1:length(mflayers)
    orderSeed = []
    for i in eachindex(order), j in eachindex(mflayers)
        if order[i] == mflayers[j]
            push!(orderSeed, idx[j])
        end
    end
    return Int.(orderSeed)
end

## Linear decrease of thickness.
function monotonicLin!(d::Array{T1,1}, α::T1) where {T1<:Float64}
    # d .*= α.^(1:length(d))
    # Only for the second active layer until the one before the substrate
    d .*= [1.0; 1.0; α.^(3:length(d)-1); 1.0]
    return nothing
end

end # module
