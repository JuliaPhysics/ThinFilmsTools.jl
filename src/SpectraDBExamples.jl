module SpectraDB

using Interpolations
using HDF5

export SL1ExpSpectrum,
	   SL1RefSpectrum,
	   SL2ExpSpectrum,
	   SL2RefSpectrum,
	   BraggSpectrum,
	   FPSpectrum,
	   HafniaSpectrum,
	   ScandiaSpectrum,
	   ScandiaSpectrum_2,
	   Info

function Info()
	tmp1 = "\n " *
	"\n " *
	"\n Available functions from Utils module:" *
	"\n " *
	"\n     SL1ExpSpectrum(λ)" *
	"\n     SL1RefSpectrum(λ)" *
	"\n     SL2ExpSpectrum(λ)" *
	"\n     SL2RefSpectrum(λ)" *
	"\n     BraggSpectrum(λ)" *
	"\n     FPSpectrum(λ)" *
	"\n     HafniaSpectrum(λ)" *
	"\n     ScandiaSpectrum(λ)" *
	"\n     ScandiaSpectrum_2(λ)" *
	"\n "
	"\n     To use any of these functions type: SpectraDB.function(args)." *
	"\n "
	return println(tmp1)
end

"""Read the HDF5 file with the database."""
function read_spectra_DB(x::String)
    datapath = joinpath(@__DIR__, "..", "data/")
    readf = h5open(joinpath(datapath, "ExampleSpectraDB.h5"), "r") do file
        read(file, x)
    end
    readf
end

"""

	Build interpolation objects from Matrix. Uses a Gridded(Linear()) grid.

		itp = build_interpolation(X)

			X: Matrix array with first column as independent variable and second column as dependent one
			itp: interpolation object

"""
function build_interpolation(X::Array{T1,2}) where {T1<:Float64}
	# Sort by ascending wavelength
	X = X[sortperm(X[:,1]), :]
    knots = (sort(vec(X[:,1])),)
    return interpolate(knots, vec(X[:,2]), Gridded(Linear()))
end

"""

	Raw experimental reflectance spectrum of a porous silicon thin film on a crystalline silicon substrate. Interpolates for a valid input wavelength range.

		R = SL1ExpSpectrum(λ)

			λ: wavelength range (nm), ∈ [200, 1125]
			R: reflectance

"""
function SL1ExpSpectrum(λ::AbstractArray{T1,1}) where {T1<:Real}
	X = read_spectra_DB("SL1sample")
    itp = build_interpolation(hcat(vec(X["lambda"]), vec(X["reflectance"])))
    return itp.(λ)
end

"""

	Reference experimental reflectance spectrum of a crystalline silicon substrate. Interpolates for a valid input wavelength range.

		R = SL1RefSpectrum(λ)

			λ: wavelength range (nm), ∈ [200, 1125]
			R: reflectance

"""
function SL1RefSpectrum(λ::AbstractArray{T1,1}) where {T1<:Real}
	X = read_spectra_DB("SL1reference")
	itp = build_interpolation(hcat(vec(X["lambda"]), vec(X["reflectance"])))
    return itp.(λ)
end

"""

	Raw experimental reflectance spectrum of a porous silicon thin film on a crystalline silicon substrate. Interpolates for a valid input wavelength range.

		R = SL2ExpSpectrum(λ)

			λ: wavelength range (nm), ∈ [200, 1125]
			R: reflectance

"""
function SL2ExpSpectrum(λ::AbstractArray{T1,1}) where {T1<:Real}
	X = read_spectra_DB("SL2sample")
	itp = build_interpolation(hcat(vec(X["lambda"]), vec(X["reflectance"])))
    return itp.(λ)
end

"""

	Reference experimental reflectance spectrum of a crystalline silicon substrate. Interpolates for a valid input wavelength range.

		R = SL2RefSpectrum(λ)

			λ: wavelength range (nm), ∈ [200, 1125]
			R: reflectance

"""
function SL2RefSpectrum(λ)
	X = read_spectra_DB("SL2reference")
	itp = build_interpolation(hcat(vec(X["lambda"]), vec(X["reflectance"])))
    return itp.(λ)
end

"""

	Absolute experimental reflectance spectrum of a porous silicon thin film photonic crystal (DBR or Bragg stack) on a crystalline silicon substrate. Interpolates for a valid input wavelength range.

		R = BraggSpectrum(λ)

			λ: wavelength range (nm), ∈ [200, 1125]
			R: reflectance

"""
function BraggSpectrum(λ::AbstractArray{T1,1}) where {T1<:Real}
	X = read_spectra_DB("Bragg")
	itp = build_interpolation(hcat(vec(X["lambda"]), vec(X["reflectance"])))
    return itp.(λ)
end

"""

	Absolute experimental reflectance spectrum of a porous silicon thin film Fabry-Perot (MicroCavity stack) on a glass substrate. Interpolates for a valid input wavelength range.

		R = FPSpectrum(λ)

			λ: wavelength range (nm), ∈ [200, 1125]
			R: reflectance

"""
function FPSpectrum(λ::AbstractArray{T1,1}) where {T1<:Real}
	X = read_spectra_DB("MicroCavity")
	itp = build_interpolation(hcat(vec(X["lambda"]), vec(X["reflectance"])))
    return itp.(λ)
end

"""

	Absolute experimental transmitance spectrum of a hafnium oxide thin film on a fused silica UV substrate. Interpolates for a valid input wavelength range.

		T = HafniaSpectrum(λ)

			λ: wavelength range (nm), ∈ [190, 1100]
			T: transmitance

"""
function HafniaSpectrum(λ::AbstractArray{T1,1}) where {T1<:Real}
	X = read_spectra_DB("Hafnia")
	itp = build_interpolation(hcat(vec(X["lambda"]), vec(X["transmitance"])))
    return itp.(λ)./100.0
end

"""

	Absolute experimental transmitance spectrum of a scandium oxide thin film on a fused silica UV substrate. Interpolates for a valid input wavelength range.

		T = ScandiaSpectrum(λ)

			λ: wavelength range (nm), ∈ [190, 1100]
			T: transmitance

"""
function ScandiaSpectrum(λ::AbstractArray{T1,1}) where {T1<:Real}
	X = read_spectra_DB("Scandia")
	itp = build_interpolation(hcat(vec(X["lambda"]), vec(X["transmitance"])))
    return itp.(λ)./100.0
end

"""

	Absolute experimental transmitance spectrum of a scandium oxide thin film on a fused silica UV substrate. Interpolates for a valid input wavelength range.

		T = ScandiaSpectrum_2(λ)

			λ: wavelength range (nm), ∈ [190, 1100]
			T: transmitance

"""
function ScandiaSpectrum_2(λ::AbstractArray{T1,1}) where {T1<:Real}
	X = read_spectra_DB("Scandia2")
	itp = build_interpolation(hcat(vec(X["lambda"]), vec(X["transmitance"])))
    return itp.(λ)./100.0
end

end
