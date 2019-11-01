module SpectraDB

using ..Utils: build_interpolation, readh5_file

export SL1ExpSpectrum,
	   SL1RefSpectrum,
	   SL2ExpSpectrum,
	   SL2RefSpectrum,
	   BraggSpectrum,
	   FPSpectrum,
	   HafniaSpectrum,
	   ScandiaSpectrum,
	   TantalaSpectrum,
	   Info

function Info()
	tmp1 = "\n " *
	"\n " *
	"\n Available functions from SpectraDB module:" *
	"\n " *
	"\n     SL1ExpSpectrum(λ)" *
	"\n     SL1RefSpectrum(λ)" *
	"\n     SL2ExpSpectrum(λ)" *
	"\n     SL2RefSpectrum(λ)" *
	"\n     BraggSpectrum(λ)" *
	"\n     FPSpectrum(λ)" *
	"\n     HafniaSpectrum(λ)" *
	"\n     ScandiaSpectrum(λ)" *
	"\n     TantalaSpectrum(λ)" *
	"\n "
	"\n     To use any of these functions type: SpectraDB.function(args)." *
	"\n "
	return println(tmp1)
end

"""

	Raw experimental reflectance spectrum of a porous silicon thin film on a crystalline
	silicon substrate. Interpolates for a valid input wavelength range.

		R = SL1ExpSpectrum(λ)

			λ: wavelength range (nm), ∈ [200, 1125]
			R: reflectance

"""
function SL1ExpSpectrum(λ::AbstractArray{T1,1}) where {T1<:Real}
	X = readh5_file("SL1sample", :SP)
    itp = build_interpolation(hcat(vec(X["lambda"]), vec(X["reflectance"])))
    return itp.(λ)
end

"""

	Reference experimental reflectance spectrum of a crystalline silicon substrate.
	Interpolates for a valid input wavelength range.

		R = SL1RefSpectrum(λ)

			λ: wavelength range (nm), ∈ [200, 1125]
			R: reflectance

"""
function SL1RefSpectrum(λ::AbstractArray{T1,1}) where {T1<:Real}
	X = readh5_file("SL1reference", :SP)
	itp = build_interpolation(hcat(vec(X["lambda"]), vec(X["reflectance"])))
    return itp.(λ)
end

"""

	Raw experimental reflectance spectrum of a porous silicon thin film on a crystalline
	silicon substrate. Interpolates for a valid input wavelength range.

		R = SL2ExpSpectrum(λ)

			λ: wavelength range (nm), ∈ [200, 1125]
			R: reflectance

"""
function SL2ExpSpectrum(λ::AbstractArray{T1,1}) where {T1<:Real}
	X = readh5_file("SL2sample", :SP)
	itp = build_interpolation(hcat(vec(X["lambda"]), vec(X["reflectance"])))
    return itp.(λ)
end

"""

	Reference experimental reflectance spectrum of a crystalline silicon substrate.
	Interpolates for a valid input wavelength range.

		R = SL2RefSpectrum(λ)

			λ: wavelength range (nm), ∈ [200, 1125]
			R: reflectance

"""
function SL2RefSpectrum(λ)
	X = readh5_file("SL2reference", :SP)
	itp = build_interpolation(hcat(vec(X["lambda"]), vec(X["reflectance"])))
    return itp.(λ)
end

"""

	Absolute experimental reflectance spectrum of a porous silicon thin film photonic
	crystal (DBR or Bragg stack) on a crystalline silicon substrate. Interpolates for
	a valid input wavelength range.

		R = BraggSpectrum(λ)

			λ: wavelength range (nm), ∈ [200, 1125]
			R: reflectance

"""
function BraggSpectrum(λ::AbstractArray{T1,1}) where {T1<:Real}
	X = readh5_file("Bragg", :SP)
	itp = build_interpolation(hcat(vec(X["lambda"]), vec(X["reflectance"])))
    return itp.(λ)
end

"""

	Absolute experimental reflectance spectrum of a porous silicon thin film Fabry-Perot
	(MicroCavity stack) on a glass substrate. Interpolates for a valid input wavelength range.

		R = FPSpectrum(λ)

			λ: wavelength range (nm), ∈ [200, 1125]
			R: reflectance

"""
function FPSpectrum(λ::AbstractArray{T1,1}) where {T1<:Real}
	X = readh5_file("MicroCavity", :SP)
	itp = build_interpolation(hcat(vec(X["lambda"]), vec(X["reflectance"])))
    return itp.(λ)
end

"""

	Absolute experimental transmitance spectrum of a hafnium oxide thin film on a fused
	silica UV substrate. Interpolates for a valid input wavelength range.

		T = HafniaSpectrum(λ)

			λ: wavelength range (nm), ∈ [190, 1100]
			T: transmitance

"""
function HafniaSpectrum(λ::AbstractArray{T1,1}) where {T1<:Real}
	X = readh5_file("Hafnia", :SP)
	itp = build_interpolation(hcat(vec(X["lambda"]), vec(X["transmitance"])))
    return itp.(λ)./100.0
end

"""

	Absolute experimental transmitance spectrum of a scandium oxide thin film on a fused
	silica UV substrate. Interpolates for a valid input wavelength range.

		T = ScandiaSpectrum(λ)

			λ: wavelength range (nm), ∈ [190, 1100]
			T: transmitance

"""
function ScandiaSpectrum(λ::AbstractArray{T1,1}) where {T1<:Real}
	X = readh5_file("Scandia", :SP)
	itp = build_interpolation(hcat(vec(X["lambda"]), vec(X["transmitance"])))
    return itp.(λ)./100.0
end

"""

	Absolute experimental ellipsomtry spectra of a tantalum oxide thin film on a fused
	silica UV substrate. Interpolates for a valid input wavelength range. Measured at 60
	degrees.

		Ψ, Δ = TantalaSpectrum(λ)

			λ: wavelength range (nm), ∈ [200, 2100]
			Ψ: psi spectrum
			Δ: delta spectrum

"""
function TantalaSpectrum(λ::AbstractArray{T1,1}) where {T1<:Real}
	X = readh5_file("Tantala", :SP)
	itp = build_interpolation(hcat(vec(X["lambda"]), vec(X["psi"])))
	psi = itp.(λ)
	itp = build_interpolation(hcat(vec(X["lambda"]), vec(X["delta"])))
	delta = itp.(λ)
    return psi, delta
end

end # module
