module Utils

using Interpolations
using HDF5
using SpecialFunctions
using LinearAlgebra # pinv

export build_interpolation,
	   _readh5_file,
	   multiple_reflections,
	   find_closest,
	   blackbody_radiation,
	   wavelength2rgb,
	   simple_moving_average,
	   savitzky_golay,
	   flatten_arrays,
	   array_arrays,
	   average_polarisation,
	   gaussian,
	   lorentzian,
	   voigtian,
	   _oscillators_input,
	   unfold_bnd,
	   Info

function Info()
    tmp1 = "\n " *
	"\n " *
    "\n Available functions from Utils module:" *
	"\n " *
	"\n     build_interpolation(X)" *
	"\n     multiple_reflections(n)" *
	"\n     find_closest(x, x0)" *
	"\n     blackbody_radiation(λ, T)" *
	"\n     wavelength2rgb(λ)" *
	"\n     simple_moving_average(x)" *
	"\n     savitzky_golay(m, n)" *
	"\n     savitzky_golay(m, n, x)" *
	"\n     flatten_arrays(x)" *
	"\n     array_arrays(x, y)" *
	"\n     average_polarisation(pol, Xp, Xs)" *
	"\n     gaussian(x, p)" *
	"\n     lorentzian(x, p)" *
	"\n     voigtian(x, p)" *
	"\n "
	"\n     To use any of these functions type: ?Utils.function" *
	"\n "
	return println(tmp1)
end

"""

	Find the index of the closest element in an array of reals x to a given real x0.

		idx = Utils.find_closest(x, x0)

			x: Data array
			x0: Element to search in x

			idx: Index

"""
function find_closest(x::AbstractArray{T1,N1}, x0::T1) where {T1<:Real, N1}
	idmin = findmin(abs.(x .- x0))[2][1]
	return idmin
end

"""Read the HDF5 file with the database."""
function _readh5_file(x::String, y::Symbol)
    datapath = joinpath(@__DIR__, "..", "data/")
	if isequal(y,:RI)
		fid = "RefractiveIndicesDB.h5"
	elseif isequal(y,:SP)
		fid = "ExampleSpectraDB.h5"
	end
    readf = h5open(joinpath(datapath,fid),"r") do file
        read(file,x)
    end
    return readf
end

"""

	Build interpolation objects from Matrix. Uses a Gridded(Linear()) grid.

		itp = Utils.build_interpolation(X)

			X: Matrix array with first column as independent variable and
               second column as dependent one

			itp: interpolation object

"""
function build_interpolation(X::Array{T1,2}) where {T1<:Float64}
	# Sort by ascending wavelength
	X = X[sortperm(X[:,1]),:]
    knots = (sort(vec(X[:,1])),)
	itp = interpolate(knots, vec(X[:,2]), Gridded(Linear()))
    return itp
end

"""

	Computes the reflection of a sequence of indices of refraction without
    interference and absorption effects.

		y = Utils.multiple_reflections(n)

			n: Vector with the indices of refraction of the involved layers

			y: Reflection

	author: Cuchu

"""
function multiple_reflections(n::AbstractArray{T1}) where {T1<:Real}
	r1::Float64 = 0.0
	r::Float64 = 0.0
	for i = 1:length(n)-1
		r2 = ((n[i] - n[i+1])/(n[i] + n[i+1]))^2
		r = (r1 + r2 - 2.0*r1*r2)/(1.0 - r1*r2)
		r1 = r
	end
	return r
end

"""

	Returns the spectral radiance (W·sr^-1·m^-3), given the wavelength (m)
    and the absolute temperature (K).

		B = Utils.blackbody_radiation(λ, T)

			λ: range of wavelengths [m]
			T: range of temperatures [K]

			B: spectral radiance [W·sr^-1·m^-3]

	Source: https://en.wikipedia.org/wiki/Planck%27s_law

"""
function blackbody_radiation(λ::AbstractArray{T1}, T::AbstractArray{T1}) where {T1<:Real}
	h::Float64 = 6.626070040e-34 # Planck's constant, J·s
	c::Float64 =  299792458.0 # speed of light, m·s
	kB::Float64 = 1.38064852e-23 # Boltzmann constant, J·K^-1
	R = Array{Float64,2}(undef,length(λ),length(T))
	num::Array = 2.0*h*c*c./λ.^5
	den::Array = h*c./(kB.*λ)
	for i in eachindex(T)
		R[:,i] .= num./expm1.(den./T[i])
	end
	return R
end

"""

	Returns an RGB matrix-color based on the wavelength input in nm.

		RBG = Utils.wavelength2rgb(λ)

			λ: range of wavelengths [nm]

	Wavelength range: 380nm and 780nm, black otherwise.

	Source: https://stackoverflow.com/questions/3407942/rgb-values-of-visible-spectrum

"""
function wavelength2rgb(λ::AbstractArray{<:Real})
    Γ::Float64 = 0.8
    # warm up
    R = Array{Float64, 1}(undef, length(λ))
    G = similar(R)
    B = similar(R)
    intensityCorrection = similar(R)
    R[380.0 .<= λ .<= 440.0] .= -(λ[380.0 .<= λ .<= 440.0] - 440.)./(440. - 380.)
    G[380.0 .<= λ .<= 440.0] .= 0.0
    B[380.0 .<= λ .<= 440.0] .= 1.0 .+ (0. * λ[380.0 .<= λ .<= 440.0])
    R[440.0 .< λ .<= 490.0] .= 0.0
    G[440.0 .< λ .<= 490.0] .= ((λ[440.0 .< λ .<= 490.0] - 440.)./(490. - 440.))
    B[440.0 .< λ .<= 490.0] .= 1.0
    R[490.0 .< λ .<= 510.0] .= 0.0
    G[490.0 .< λ .<= 510.0] .= 1.0
    B[490.0 .< λ .<= 510.0] .= (-(λ[490.0 .< λ .<= 510.0] - 510.)./(510. - 490.))
    R[510.0 .< λ .<= 580.0] .= ((λ[510.0 .< λ .<= 580.0] - 510.)./(580. - 510.))
    G[510.0 .< λ .<= 580.0] .= 1.0
    B[510.0 .< λ .<= 580.0] .= 0.0
    R[580.0 .< λ .<= 645.0] .= 1.0
    G[580.0 .< λ .<= 645.0] .= (-(λ[580.0 .< λ .<= 645.0] - 645.)./(645. - 580.))
    B[580.0 .< λ .<= 645.0] .= 0.0
    R[645.0 .< λ .<= 780.0] .= 1.0
    G[645.0 .< λ .<= 780.0] .= 0.0
    B[645.0 .< λ .<= 780.0] .= 0.0
    R[λ .> 780] .= 0.0
    R[λ .< 380] .= 0.0
    G[λ .> 780] .= 0.0
    G[λ .< 380] .= 0.0
    B[λ .> 780] .= 0.0
    B[λ .< 380] .= 0.0
    # LET THE INTENSITY SSS FALL OFF NEAR THE VISION LIMITS
    intensityCorrection[λ .> 700.0] .= 0.3 .+ 0.7.*(780.0 .- λ[λ .> 700.0])./(780. - 700.)
    intensityCorrection[λ .< 420.0] .= 0.3 .+ 0.7.*(λ[λ .< 420.0] - 380.0)./(420. - 380.)
    intensityCorrection[420.0 .< λ .< 700.0] .= 1.0.*λ[420.0 .< λ .< 700.0]
    # GAMMA ADJUST
    R .= (intensityCorrection.*R).^Γ
    G .= (intensityCorrection.*G).^Γ
    B .= (intensityCorrection.*B).^Γ
    # # Multiply by 255
    # R = R * 255
    # G = G * 255
    # B = B * 255
    return [R G B]
end

"""

	Smooth vector with simple moving average model.

		y = Utils.simple_moving_average(x; w=1)

			x: Vector to be smoothed out
			w: Integer with the window size

			y: Smoothed output

"""
function simple_moving_average(y::AbstractVector{T1}; w::T2=1) where {T1<:Real, T2<:Int64}
    y2 = [y[1]*ones(w); y; y[end]*ones(w)]
	m = 2*w + 1
	w = Int(w + 1)
	y_ = similar(y)
    for i in w:length(y2)-(w-1)
        y_[i-w+1] = sum(y2[i-w+1:i+w-1])
    end
    return y_ ./ m
end

"""

	Implementation of the Savitzky-Golay filter of window half-width M and degree N.

		y = Utils.savitzky_golay(m, n, x)
			m: is the number of points before and after to interpolate, the full width
			   of the window is 2m+1.
			n: is the polynomial order (must be lower than window size).
			x: unfiltered data vector.
			y: filtered data vector.

		sgf = savitzky_golay(m, n)
				m: is the number of points before and after to interpolate, the full width
				   of the window is 2m+1.
				n: is the polynomial order (must be lower than window size).
				sgf: Savitzky-Golay filter object.
					 Then you can use it as: sgf(x)

	Adapted from https://gist.github.com/jiahao/b8b5ac328c18b7ae8a17

	Notice: this version does not handle well the boundaries of the input vector.

"""
function savitzky_golay(m::T0, n::T0, x::Array{T1,1}) where {T0<:Int64, T1<:Float64}
	n < m || throw("The window size m must be larger than polynomial order n.")
    return savitzky_golayFilter{m,n}()(x)
end

function savitzky_golay(m::T0, n::T0) where {T0<:Int64}
	n < m || throw("The window size m must be larger than polynomial order n.")
    return savitzky_golayFilter{m,n}()
end

struct savitzky_golayFilter{M,N} end
@generated function (::savitzky_golayFilter{M,N})(data::AbstractVector{T}) where {M,N,T}
    # Create Jacobian matrix
    J = zeros(2M+1, N+1)
    for i=1:2M+1, j=1:N+1
        J[i, j] = (i-M-1)^(j-1)
    end
    e₁ = zeros(N+1)
    e₁[1] = 1.0
    # Compute filter coefficients
    C = J' \ e₁
    # Evaluate filter on data matrix
    To = typeof(C[1] * one(T)) # Calculate type of output
    expr = quote
        n = size(data, 1)
        smoothed = zeros($To, n)
        @inbounds for i in eachindex(smoothed)
            smoothed[i] += $(C[M+1])*data[i]
        end
        smoothed
    end
    for j=1:M
        insert!(expr.args[6].args[3].args[2].args, 1,
            :(if i - $j ≥ 1
                smoothed[i] += $(C[M+1-j])*data[i-$j]
            end)
        )
        push!(expr.args[6].args[3].args[2].args,
            :(if i + $j ≤ n
                smoothed[i] += $(C[M+1+j])*data[i+$j]
            end)
        )
    end
    return expr
end

"""

    Convert Array{Array{Real,1},1} to an Array{Float64,1}.

"""
function flatten_arrays(xin::Array{Array{T1,1},1}) where {T1<:Real}
    xout = []
    # Determine the total number of elements of xin
    for i in eachindex(1:length(xin)), j in eachindex(1:length(xin[i]))
        push!(xout, xin[i][j])
    end
    return float.(xout)
end

"""

    Convert Array{Real,1} to Array{Array{Float64,1},1}.

"""
function array_arrays(x::Array{T1,1}, _x::Array{Array{T1,1},1}) where{T1<:Real, T2<:Real}
    xfinal = deepcopy(_x)
    ki::Int64 = 1
    kf::Int64 = 0
    for i in eachindex(collect(1:length(xfinal)))
        kf += length(xfinal[i])
        xfinal[i][1] = x[ki]
        xfinal[i][2:end] = x[ki+1:kf]
        ki += length(xfinal[i])
    end
    return float.(xfinal)
end

"""

    Return the polarisation averaged spectrum.

        X = Utils.average_polarisation(pol, Xp, Xs)

            pol: indicates the polarisation
                 (pol=1.0 => p, pol=0.0 => s, 0.0<=pol<=1.0 => average)
            Xp: p/TM polarisation quantity (e.g. Rp)
            Xs: s/TE polarisation quantity (e.g. Rs)

            Xavg: unpolarized quantity = pol*Xp + (1.0 - pol)*Xs

"""
function average_polarisation(pol::T1, Xp::Array{T1}, Xs::Array{T1}) where {T1<:Float64}
	X = @. pol*Xp + (1.0 - pol)*Xs
    return X
end

"""

	Multipeak Gaussian PDF curve.

		modelg = Utils.gaussian(x, p)

			x = vector with data of the x-axis
			p = [
				  [p0], # offset
				  [A1, μ1, σ1], # first peak
				  [A2, μ2, σ2], # second peak
				  ..., # so on
				]
				p0: offset of the total curve
				Ai: amplitude of peak i
				μi: position of peak i (mean or expectation of the distribution)
				σi: standard deviation of peak i
				i = 1,...,n

	Source: https://en.wikipedia.org/wiki/Normal_distribution

"""
function gaussian(x::AbstractArray{T1}, p::Array{Array{T1,1},1}) where {T1<:Real}
    isequal(length(flatten_arrays(p)[2:end]), 3*(length(p)-1)) || throw("The inputs are not correct. For the Gaussian curve, each peak has three parameters plus the offset.")
    y = ones(length(x)).*p[1][1] # offset
    for i = 2:length(p)
        @. y += p[i][1]*exp(-((x - p[i][2])/p[i][3])^2)
    end
    return y
end

"""

	Multipeak Cauchy-Lorentz PDF curve.

		modell = Utils.lorentzian(x, p)

			x = vector with data of the x-axis
			p = [
				  [p0], # offset
				  [A1, μ1, Γ1], # first peak
				  [A2, μ2, Γ2], # second peak
				  ..., # so on
				]
				p0: offset of the total curve
				Ai: amplitude of peak i
				μi: position of peak i (location parameter)
				Γi: half-width at half-maximum (HWHM) of peak i
				i = 1,...,n

	Source: https://en.wikipedia.org/wiki/Cauchy_distribution

"""
function lorentzian(x::AbstractArray{T1}, p::Array{Array{T1,1},1}) where {T1<:Real}
    isequal(length(flatten_arrays(p)[2:end]), 3*(length(p)-1)) || throw("The inputs are not correct. For the Lorentzian curve, each peak has three parameters plus the offset.")
    y = ones(length(x)).*p[1][1] # offset
    for i = 2:length(p)
        @. y += p[i][1]/(1.0 + ((x - p[i][2])/p[i][3])^2)
    end
    return y
end

"""

	Single peak Voigt PDF curve. The function is computed using the Faddeeva's
    function through the scaled complementary error function of x.

		modelv = Utils.voigtian(x, p)

			x = vector with data of the x-axis
			p = [
				  [p0], # offset
				  [A1, μ1, σ1, Γ1], # first peak
				  [A2, μ2, σ2, Γ2], # second peak
				  ..., # so on
				]
				p0: offset of the total curve
				Ai: amplitude of the curve
				μi: position of peak i (location parameter)
				σi: Gaussian standard deviation of peak i
				Γi: Lorentzian half-width at half-maximum (HWHM) of peak i
				i = 1,...,n

		The μ parameter has been included in a custom fashion to allow the
        distribution to be uncentered.

	Source: https://en.wikipedia.org/wiki/Voigt_profile

"""
function voigtian(x::AbstractArray{T1}, p::Array{Array{T1,1},1}) where {T1<:Real}
    isequal(length(flatten_arrays(p)[2:end]), 4*(length(p)-1)) || throw("The inputs are not correct. For the Voigtian curve, each peak has four parameters plus the offset.")
    y = ones(length(x)).*p[1][1] # offset
    for i = 2:length(p)
        @. y += p[i][1]*real(erfcx(-im*((x - p[i][2] + im*p[i][3])/sqrt(2)/p[i][4])))
    end
    return y
end

"""

	Build the input for oscillators models. Each index of refraction associated with
	Lorentz can be modelled with arbitrary numbers of oscillators. This function provides a
	way to construct a valid input for them.

	The input for the following models can be constructed are:
		RI.drudelorentz, RI.tauclorentz, RI.forouhibloomer,
		RI.lorentzplasmon, RI.codylorentz

		x = _oscillators_input(argin, n)

			argin: Array with two elements.
				   The first element is an array with the parameters not associated to
				   oscillators, while the second element contains an array with all the
				   parameters of the oscillators.
			n: number of parameters in each oscillator
			x: Array to input in RI modelled with Lorentz.

	For instance, the Tauc-Lorentz model accepts inputs with two parameters
	(length(argin[1]) == 2) not associated with the oscillator and three (n=3) forming
	each oscillator. Say you want to pass a three oscillator model input to the
	Tauc-Lorentz:
		argin = [
			[4.0, 3.5],  # ϵinf, Eg
			[15.0, 6.0, 1.05, # osc 1, A1, E01, Γ1
			 8.41, 4.5, 0.2, # osc 2, A2, E02, Γ2
			 1.5, 4.6, 0.1, # osc 3, A3, E03, Γ3
			],
		]
		x = _oscillators_input(argin, 3)
		x = [
			  [4.0, 3.5],  # ϵinf, Eg
			  [15.0, 6.0, 1.05], # osc 1, A1, E01, Γ1
			  [8.41, 4.5, 0.2], # osc 2, A2, E02, Γ2
			  [1.5, 4.6, 0.1], # osc 3, A3, E03, Γ3
		]
		RI.tauclorentz(x, λ)

	Notice: Careful with the number of parameteres in each model.

"""
## Build the input for oscillators models.
# Each element of x has the parameters for the corresponding oscillator.
# function _oscillators_input(argin, n::Int64)
#     x = argin[1]
#     numel::Int64 = Int(length(argin[2:end])/n + 1)
#     k::Int64 = 0
#     for i = 2:numel
#         push!(x, argin[i+k:n+i-1+k])
#         k += n-1
#     end
#     return x
# end
function _oscillators_input(argin, n::Int64)
    x = [argin[1]]
    numel::Int64 = Int(length(argin[2])/n)
    k::Int64 = 0
    for i = 1:numel
        push!(x, argin[2][i+k:n+i-1+k])
        k += n-1
    end
    return x
end

"""

	unfold_bnd(lb, ub)

Returns two arrays with the lower and upper bounds in accepted format for BBO package.

"""
function unfold_bnd(lb, ub)
    return (collect(zip(flatten_arrays(lb), flatten_arrays(ub))))
end

end # Utils
