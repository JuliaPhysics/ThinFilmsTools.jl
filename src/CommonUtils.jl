module Utils

export multipleReflections,
	   bbRadiation,
	   wavelength2RGB,
	   movingAverage,
	   flattenArrays,
	   arrayArrays,
	   averagePolarisation,
	   Info

function Info()
    tmp1 = "\n " *
	"\n " *
    "\n Available functions from Utils module:" *
	"\n " *
	"\n     multipleReflections(n)" *
	"\n     bbRadiation(λ, T)" *
	"\n     wavelength2RGB(λ)" *
	"\n     movingAverage(x, w)" *
	"\n     flattenArrays(x)" *
	"\n     arrayArrays(x, y)" *
	"\n     averagePolarisation(pol, Xp, Xs)" *
	"\n "
	"\n     To use any of these functions type: Utils.function(args)." *
	"\n "
	return println(tmp1)
end

"""

	Computes the reflection of a sequence of indices of refraction without interference and absorption effects.

		y = multipleReflections(n)

			n: Vector with the indices of refraction of the involved layers

	authors: cuchu & lnacquaroli

"""
function multipleReflections(n::AbstractArray{T1}) where {T1<:Real}
  R1::Float64 = 0.0
  R::Float64 = 0.0
  for i = 1:length(n)-1
    R2 = ( (n[i] - n[i+1])/(n[i] + n[i+1]))^2
    R = (R1 + R2 - 2.0*R1*R2)/(1.0 - R1*R2)
    R1 = R
  end
  return R
end

"""

	Returns the spectral radiance (W·sr^-1·m^-3), given the wavelength (m) and the absolute temperature (K).

		B = bbRadiation(λ, T)

			λ: range of wavelengths [m]
			T: range of temperatures [K]
			B: spectral radiance [W·sr^-1·m^-3]

	Source: https://en.wikipedia.org/wiki/Planck%27s_law

"""
function bbRadiation(λ::AbstractArray{T1}, T::AbstractArray{T1}) where {T1<:Real}
    h::Float64 = 6.626070040e-34 # Planck's constant, J·s
    c::Float64 =  299792458.0 # speed of light, m·s
    kB::Float64 = 1.38064852e-23 # Boltzmann constant, J·K^-1
    R = Array{Float64, 2}(undef, length(λ), length(T))
	num::Array = 2.0*h*c*c./λ.^5
    den::Array = h*c./(kB.*λ)
    for i in eachindex(T)
        R[:,i] .= num./expm1.(den./T[i])
    end
    return R
end

"""

	Returns an RGB matrix-color based on the wavelength input in nm.

		RBG = wavelength2RGB(λ)

			λ: range of wavelengths [nm]

	Wavelength range: 380nm and 780nm, black otherwise.

	Source: https://stackoverflow.com/questions/3407942/rgb-values-of-visible-spectrum

"""
function wavelength2RGB(λ::AbstractArray{<:Real})
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

	Smooth vector with moving average model.

		y = movingAverage(x, w)

			x: Vector to be smoothed out
			w: Integer with the window size
			y: Smoothed output

"""
function movingAverage(x::AbstractVector{T1}, w::T2) where {T1<:Real, T2<:Int64}
	xLen::Int64 = length(x)
    y = Vector{Float64}(undef, xLen)
    for i in eachindex(x)
        lo = max(1, i - w)
        hi = min(xLen, i + w)
        y[i] = sum(x[lo:hi])./xLen
    end
    return y
end

"""

    Convert Array{Array{Real,1},1} to an Array{Real,1}.

"""
function flattenArrays(xin::Array{Array{T1,1},1}) where {T1<:Real}
    xout = []
    # Determine the total number of elements of xinit
    for i in eachindex(1:length(xin)), j in eachindex(1:length(xin[i]))
        push!(xout, xin[i][j])
    end
    return float.(xout)
end

"""

    Convert Array{Real,1} to Array{Array{Real,1},1}.

"""
function arrayArrays(x::Array{T1,1}, _x::Array{Array{T2,1},1}) where{T1<:Float64, T2<:Real}
    xfinal = deepcopy(_x)
    ki::Int64 = 1
    kf::Int64 = 0
    for i in eachindex(collect(1:length(xfinal)))
        kf += length(xfinal[i])
        xfinal[i][1] = x[ki]
        xfinal[i][2:end] = x[ki+1:kf]
        ki += length(xfinal[i])
    end
    return xfinal
end

"""

    Return the polarisation averaged spectrum.

        X = averagePolarisation(pol, Xp, Xs)

            pol: indicates the polarisation (1.0 = p, 0.0 = s, between 0.0 and 1.0 = average)
            Xp: p/TM polarisation quantity (e.g. Rp)
            Xs: s/TE polarisation quantity (e.g. Rs)
            Xavg: unpolarized quantity = pol*Rp + (1.0 - pol)*Rs

"""
function averagePolarisation(pol::T1, Xp::Array{T1}, Xs::Array{T1}) where {T1<:Float64}
    return pol.*Xp .+ (1.0 - pol).*Xs
end

end # CommonUtils
