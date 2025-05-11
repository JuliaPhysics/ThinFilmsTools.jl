module SpectraDB

using ..Utils: build_interpolation, _readh5_file

export sl1_exp_spectrum,
       sl1_ref_spectrum,
       sl2_exp_spectrum,
       sl2_ref_spectrum,
       bragg_spectrum,
       fp_spectrum,
       hafnia_spectrum,
       scandia_spectrum,
       tantala_spectrum,
       hafnia_ellips,
       Info

function Info()
    tmp1 = "\n " *
    "\n " *
    "\n Available functions from SpectraDB module:" *
    "\n " *
    "\n     sl1_exp_spectrum(λ)" *
    "\n     sl1_ref_spectrum(λ)" *
    "\n     sl2_exp_spectrum(λ)" *
    "\n     sl2_ref_spectrum(λ)" *
    "\n     bragg_spectrum(λ)" *
    "\n     fp_spectrum(λ)" *
    "\n     hafnia_spectrum(λ)" *
    "\n     scandia_spectrum(λ)" *
    "\n     tantala_spectrum(λ)" *
    "\n "
    "\n     To use any of these functions type: SpectraDB.function(args)." *
    "\n "
    return println(tmp1)
end

"""

    Raw experimental reflectance spectrum of a porous silicon thin film on a crystalline
    silicon substrate. Interpolates for a valid input wavelength range.

        R = sl1_exp_spectrum(λ)

            λ: wavelength range (nm), ∈ [200, 1125]
            R: reflectance

"""
function sl1_exp_spectrum(λ::AbstractVector{T1}) where {T1<:Real}
    X = _readh5_file("SL1sample", :SP)
    itp = build_interpolation(hcat(vec(X["lambda"]), vec(X["reflectance"])))
    return itp.(λ)
end

"""

    Reference experimental reflectance spectrum of a crystalline silicon substrate.
    Interpolates for a valid input wavelength range.

        R = sl1_ref_spectrum(λ)

            λ: wavelength range (nm), ∈ [200, 1125]
            R: reflectance

"""
function sl1_ref_spectrum(λ::AbstractVector{T1}) where {T1<:Real}
    X = _readh5_file("SL1reference", :SP)
    itp = build_interpolation(hcat(vec(X["lambda"]), vec(X["reflectance"])))
    return itp.(λ)
end

"""

    Raw experimental reflectance spectrum of a porous silicon thin film on a crystalline
    silicon substrate. Interpolates for a valid input wavelength range.

        R = sl2_exp_spectrum(λ)

            λ: wavelength range (nm), ∈ [200, 1125]
            R: reflectance

"""
function sl2_exp_spectrum(λ::AbstractVector{T1}) where {T1<:Real}
    X = _readh5_file("SL2sample", :SP)
    itp = build_interpolation(hcat(vec(X["lambda"]), vec(X["reflectance"])))
    return itp.(λ)
end

"""

    Reference experimental reflectance spectrum of a crystalline silicon substrate.
    Interpolates for a valid input wavelength range.

        R = sl2_ref_spectrum(λ)

            λ: wavelength range (nm), ∈ [200, 1125]
            R: reflectance

"""
function sl2_ref_spectrum(λ)
    X = _readh5_file("SL2reference", :SP)
    itp = build_interpolation(hcat(vec(X["lambda"]), vec(X["reflectance"])))
    return itp.(λ)
end

"""

    Absolute experimental reflectance spectrum of a porous silicon thin film photonic
    crystal (DBR or Bragg stack) on a crystalline silicon substrate. Interpolates for
    a valid input wavelength range.

        R = bragg_spectrum(λ)

            λ: wavelength range (nm), ∈ [200, 1125]
            R: reflectance

"""
function bragg_spectrum(λ::AbstractVector{T1}) where {T1<:Real}
    X = _readh5_file("Bragg", :SP)
    itp = build_interpolation(hcat(vec(X["lambda"]), vec(X["reflectance"])))
    return itp.(λ)
end

"""

    Absolute experimental reflectance spectrum of a porous silicon thin film Fabry-Perot
    (MicroCavity stack) on a glass substrate. Interpolates for a valid input wavelength range.

        R = fp_spectrum(λ)

            λ: wavelength range (nm), ∈ [200, 1125]
            R: reflectance

"""
function fp_spectrum(λ::AbstractVector{T1}) where {T1<:Real}
    X = _readh5_file("MicroCavity", :SP)
    itp = build_interpolation(hcat(vec(X["lambda"]), vec(X["reflectance"])))
    return itp.(λ)
end

"""

    Absolute experimental transmitance spectrum of a hafnium oxide thin film on a fused
    silica UV substrate. Interpolates for a valid input wavelength range.

        T = hafnia_spectrum(λ)

            λ: wavelength range (nm), ∈ [190, 1100]
            T: transmitance

"""
function hafnia_spectrum(λ::AbstractVector{T1}) where {T1<:Real}
    X = _readh5_file("Hafnia", :SP)
    itp = build_interpolation(hcat(vec(X["lambda"]), vec(X["transmitance"])))
    return itp.(λ)./100.0
end

"""

    Absolute experimental transmitance spectrum of a scandium oxide thin film on a fused
    silica UV substrate. Interpolates for a valid input wavelength range.

        T = scandia_spectrum(λ)

            λ: wavelength range (nm), ∈ [190, 1100]
            T: transmitance

"""
function scandia_spectrum(λ::AbstractVector{T1}) where {T1<:Real}
    X = _readh5_file("Scandia", :SP)
    itp = build_interpolation(hcat(vec(X["lambda"]), vec(X["transmitance"])))
    return itp.(λ)./100.0
end

"""

    Absolute experimental ellipsomtry spectra of a tantalum oxide thin film on a fused
    silica UV substrate. Interpolates for a valid input wavelength range. Measured at 60
    degrees.

        Ψ, Δ = tantala_spectrum(λ)

            λ: wavelength range (nm), ∈ [200, 2100]
            Ψ: psi spectrum
            Δ: delta spectrum

"""
function tantala_spectrum(λ::AbstractVector{T1}) where {T1<:Real}
    X = _readh5_file("Tantala", :SP)
    itp = build_interpolation(hcat(vec(X["lambda"]), vec(X["psi"])))
    psi = itp.(λ)
    itp = build_interpolation(hcat(vec(X["lambda"]), vec(X["delta"])))
    delta = itp.(λ)
    return psi, delta
end

"""

    Absolute experimental ellipsomtry spectra of a tantalum oxide thin film on a fused
    silica UV substrate. Interpolates for a valid input wavelength range. Measured at 60
    degrees.

        Ψ, Δ = hafnia_ellips(λ)

            λ: wavelength range (nm), ∈ [200, 2100]
            Ψ: psi spectrum
            Δ: delta spectrum

"""
function hafnia_ellips(λ::AbstractVector{T1}) where {T1<:Real}
    X = _readh5_file("HafniaEllips", :SP)
    itp = build_interpolation(hcat(1240.0./vec(X["lambda"]), vec(X["psi"])))
    psi = itp.(λ)
    itp = build_interpolation(hcat(vec(1240.0./X["lambda"]), vec(X["delta"])))
    delta = itp.(λ)
    return psi, delta
end

end # module
