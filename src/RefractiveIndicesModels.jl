module RI

export bruggemanspheresbin, looyengacylindersbin, looyengaspheresbin, lorentzlorenzbin, maxwellgarnettspheresbin, moneckebin, gedfbin, gembin, Info

function Info()
    tmp1 = "\n " *
	"\n " *
    "\n Available indices of refraction models:" *
	"\n " *
    "\n     bruggemanspheresbin(p, n1, n2)" *
    "\n     looyengacylindersbin(p, n1, n2)" *
    "\n     looyengaspheresbin(p, n1, n2)" *
    "\n     lorentzlorenzbin(p, n1, n2)" *
    "\n     maxwellgarnettspheresbin(p, n1, n2)" *
    "\n     moneckebin(p, n1, n2)" *
    "\n     gedfbin([p, β], n1, n2)" *
    "\n     gembin([ϕ, ϕc, tp], n1, n2)" *
    "\n     sellmeier([B1, B2, B3, C1, C2, C3], λ)" *
    "\n     cauchyurbach([A, B, C, D, E, E0], λ)" *
	"\n     lorentzdispersion([ϵinf, ħωp, ħωt, Γ], λ)" *
	"\n     tauclorentz([A, E0, C, Eg, ϵinf], λ)" *
	"\n     forouhibloomer([A, B, C, Egap, Ninf], λ)" *
	"\n     forouhibloomermodified([ninf, ħω0, Γ, f, Eg], λ)" *
	"\n "
	return println(tmp1)
end

"""

    Calculates the complex refractive index of a porous material using Bruggeman model, for nonmagnetic and isotropic materials. It considers spheres (n1) hosted in another material (n2).

        N = RI.bruggemanspheresbin(p, n1, n2)

            n1: refractive index of material 1
            n2: refractive index of material 2
            p: porosity parameter (real number, 0<p<1), proportion of n1 over total
            N: effective refractive index

    Source: PHYSICAL REVIEW B VOLUME 61, NUMBER 15 15 APRIL 2000-I.

"""
function bruggemanspheresbin(p::S, n1::Array{T}, n2::Array{T}) where {S<:Float64, T<:ComplexF64}
    df1 = n1.^2
    df2 = n2.^2
    return @. sqrt(0.25*(-df1 + 2.0*df2 + 3.0*df1*p - 3.0*df2*p + sqrt(8.0*df1*df2 + (-df1 + 2.0*df2 + 3.0*df1*p - 3.0*df2*p)^2)))
end

"""

    Calculates the complex refractive index of the porous material using Looyenga-Landau-Lifshitz model, for nonmagnetic and isotropic materials. It considers cylinders (n1) hosted in another material (n2).

       N = RI.looyengacylindersbin(p, n1, n2)

            n1: refractive index of material 1
            n2: refractive index of material 2
            p: porosity parameter (real number, 0<p<1), proportion of n1 over total
            N: effective refractive index

    Source: Langmuir 2013, 29, 2784−2789.

"""
function looyengacylindersbin(p::S, n1::Array{T}, n2::Array{T}) where {S<:Float64, T<:ComplexF64}
    return @. (1.0 - p)*n2^2 + n1^2*p
end

"""

    Calculates the complex refractive index of the porous material using Looyenga-Landau-Lifshitz model, for nonmagnetic and isotropic materials. It considers spheres (n1) hosted in another material (n2).

        N = RI.looyengaspheresbin(p, n1, n2)

            n1: refractive index of material 1
            n2: refractive index of material 2
            p: porosity parameter (real number, 0<p<1), proportion of n1 over total
            N: effective refractive index

    Source: IEEE TRANSACTIONS ON GEOSCIENCE AND REMOTE SENSING, VOL. 38, NO. 3, MAY 2000.

"""
function looyengaspheresbin(p::S, n1::Array{T}, n2::Array{T}) where {S<:Float64, T<:ComplexF64}
    df1 = n1.^2
    df2 = n2.^2
    return @. sqrt((((1.0 - p)*(df2^(1/3.))) + ((df1^(1/3.))*p))^3)
end

"""

    Returns the effective dielectric function of a binary liquid mixture. p is the proportion of component n1.

        N = RI.lorentzlorenzbin(p, n1, n2)

            n1: refractive index of material 1
            n2: refractive index of material 2
            p: proportion of n1 over total, parameter (real number, 0<p<1)
            N: effective refractive index

"""
function lorentzlorenzbin(p::S, n1::Array{T}, n2::Array{T}) where {S<:Float64, T<:ComplexF64}
    return @. sqrt((n2^2*(-n1^2 + 2.0*p^2 - 2.0 - 2.0*p*n1^2)) / (n1^2*(p - 1.0) - 2.0 - p*n2^2))
end

"""

    Calculates the complex refractive index of the porous material using Maxwell-Garnett model, for nonmagnetic and isotropic materials. It considers spheres (n1) hosted in another material (n2).

        N = RI.maxwellgarnettspheresbin(p, n1, n2)

            n1: refractive index of material 1
            n2: refractive index of material 2
            p: porosity parameter (real number, 0<p<1), proportion of n1 over total
            N: effective refractive index

    Source: IEEE TRANSACTIONS ON GEOSCIENCE AND REMOTE SENSING, VOL. 38, NO. 3, MAY 2000.

"""
function maxwellgarnettspheresbin(p::S, n1::Array{T}, n2::Array{T}) where {S<:Float64, T<:ComplexF64}
    df1 = n1.^2
    df2 = n2.^2
    return @. sqrt(df1 + ((2.0*p*df1)*(df2 - df1) / (df2 + df1 - p*(df2 - df1))))
end

"""

    Calculates the complex refractive index of the porous material using Monecke model, for nonmagnetic and isotropic materials.

        N = RI.moneckebin(p, n1, n2)

            n1: refractive index of material 1
            n2: refractive index of material 2
            p: porosity parameter (real number, 0<p<1), proportion of n1 over total
            N: effective refractive index

    Source: Phys. Rev. B, Vol. 61, Num. 15, 15 April, 2000-I.

"""
function moneckebin(p::S, n1::Array{T}, n2::Array{T}) where {S<:Float64, T<:ComplexF64}
    df1 = n1.^2
    df2 = n2.^2
    return @. sqrt((2.0*(p*df1 + (1.0 - p)*df2)^2) + (df1*df2) / ((1.0 + p)*df1) + ((2.0 - p)*df2))
end

"""

    Calculates the complex refractive index of a porous material using a general dielectric model, for nonmagnetic and isotropic materials. It considers n1) hosted in another material (n2).

        N = RI.gedfbin(x, n1, n2)

        	x: Array with
        		p: porosity parameter (real number, 0<p<1), proportion of n1 over total
	            β: either 1 or 3
            n1: refractive index of filling material
            n2: refractive index of host material
            N:  effective refractive index

    Source: Appl. Phys. Lett. 108, 102902 (2016); doi: 10.1063/1.4943639.

"""
function gedfbin(x::AbstractArray{S}, n1::Array{T}, n2::Array{T}) where {S<:Float64, T<:ComplexF64}
	p, β = x
    df1 = n1.^2
    df2 = n2.^2
    return @. sqrt(df2*(df1 + (β*df2) + (β*p*df1) - (β*p*df2)) / (df1 + (β*df2) - (p*df1) + (p*df2)))
end

"""

    Here we used equation (3) and the parameters ϕ, ϕc and tp as fitted. For a more detailed analysis and relation to the depolarization factors see the reference. It considers (n1) hosted in another material (n2).

        N = RI.gembin(x, n1, n2)

        	x: Array containing
        		ϕ: porosity parameter (real number, 0<ϕ<1), proportion of n1 over total
        		ϕc: critical volume fraction at which the high-index phase first percolates
        		tp: exponent related both to ϕc and to the shape of the grains
        	n1: refractive index of filling material (lower index)
            n2: refractive index of host material (higher index)
            N: effective refractive index

    Source: Carbon 40 (2002) 2801–2815.

"""
function gembin(x::AbstractArray{S}, n1::Array{T}, n2::Array{T}) where {S<:Float64, T<:ComplexF64}
	ϕ, ϕc, tp = x
    α = (1-ϕc)/ϕc
    β = 1/tp
    df1 = (n1.^2).^β
    df2 = (n2.^2).^β
    return @. sqrt(((α*df1 - ϕ*df1 - α*ϕ*df1 - df2 + ϕ*df2 + α*ϕ*df2 + sqrt(4.0*α*df1*df2 + ((α*(ϕ-1) + ϕ)*df1 + df2 - (1.0 + α)*ϕ*df2)^2))/(2.0*α))^tp)
end

"""

    Returns the index of refraction using the Sellmeier equation.

        N = RI.sellmeier(x, λ)

            x: Array containing the six values for the Sellmeier equation
               as described in https://en.wikipedia.org/wiki/Sellmeier_equation.
               (x = [B1, B2, B3, C1, C2, C3])
            λ: range of wavelengths in μm
            N: complex index of refraction

    Source: https://en.wikipedia.org/wiki/Sellmeier_equation

"""
function sellmeier(x::AbstractArray{T}, λ::Array{T}) where {T<:Float64}
    λ = λ.^2
    return @. sqrt(complex(1.0 + (x[1]*λ/(λ - x[4])) + (x[2]*λ/(λ - x[5])) + (x[3]*λ/(λ - x[6]))))
end

"""

    Returns the index of refraction using the Cauchy-Urbach model.

        N = RI.cauchyurbach(x, λ)

            x: Array containing the six values used in the model.
               (x = [An, Bn, Cn, Ak, Bk, Ck])
            λ: range of wavelengths in nm
            N: complex index of refraction

    Source: E. Krous. Characterization of scandium oxide thin films for use in interference coatings for high-power lasers operating in the near-infrared. Master’s thesis, Colorado State University, 2010. Page 60.

"""
function cauchyurbach(x::AbstractArray{T}, λ::Array{T}) where {T<:Float64}
	An, Bn, Cn, Ak, Bk, Ck = x
    return @. An + Bn*1e4/λ^2 + Cn*1e9/λ^4 + im*Ak*exp(Bk + (1240.0/λ/Ck))
end

"""

    Returns the index of refraction using the single Lorentz oscillator model.

        N = RI.lorentzdispersion(x, λ)

            x: Array containing the three values used in the model.
               (x = [ϵinf, ħωp, ħωt, Γ])
            λ: range of wavelengths in nm
            N: complex index of refraction

    Sources: F. Wooten, Optical properties of solids, Academic Press (1972)
			 http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Lorentz_Dispersion_Model.pdf

"""
function lorentzdispersion(x::AbstractArray{T}, λ::Array{T}) where {T<:Float64}
    ϵinf, ħωp, ħωt, Γ = x
	ϵ = 1240.0./λ
	return @. sqrt(ϵinf + (ħωp^2/(ħωt^2 - ϵ^2 + im*Γ*ϵ)))
end


"""

    Returns the index of refraction using the Tauc-Lorentz dispersion model.

        N = RI.tauclorentz(x, λ)

            x: Array containing the five values used in the model.
               (x = [A, E0, C, Eg, ϵinf], where E0, Eg and C are in eV units.)
			   Make sure that E0^2 - C^2/2 > 0!
            λ: range of wavelengths in nm
            N: complex index of refraction

    Sources: Appl. Phys. Lett. 69 (3), 1996
    		 Appl. Phys. Lett. 69 (14), 1996
             J. Phys.: Condens. Matter 20 015216, 2008 (doi.org/10.1088/0953-8984/20/01/015216)

"""
function tauclorentz(x::AbstractArray{T}, λ::Array{T}) where {T<:Float64}
	A, E0, C, Eg, ϵinf = x
	ħω = 1240.0 ./ λ
	ħω2 = ħω.^2; E02 = E0^2; Eg2 = Eg^2; C2 = C^2
	γ2 = real(sqrt(complex(E02 - C2/2)))^2
	α = real(sqrt(complex(4.0*E02 - C2)))
	ζ4 = @. (ħω2 - γ2)^2 + α^2*C2/4
	aL = @. (Eg2 - E02)*ħω2 + Eg2*C2 - E02*(E02 + 3.0*Eg2)
	aA = @. (ħω2 - E02)*(E02 + Eg2) + Eg2*C2
	ϵf1 = @. ϵinf +
			 (A*C*aL/2.0/π/ζ4/α/E0*log((E02 + Eg2 + α*Eg)/(E02 + Eg2 - α*Eg)) -
			 A*aA/π/ζ4/E0*(π - atan((2.0*Eg + α)/C) + atan((α - 2.0*Eg)/C)) +
			 2.0*A*E0*Eg/π/ζ4/α*(ħω2 - γ2)*(π + 2.0*atan(2.0/α/C*(γ2 - Eg2))) -
			 A*E0*C*(ħω2 + Eg2)/π/ζ4/ħω*log(abs(ħω - Eg)/(ħω + Eg)) +
			 2.0*A*E0*C*Eg/π/ζ4*log((abs(ħω - Eg)*(ħω + Eg))/sqrt((E02 - Eg2)^2 + Eg2*C2)))
	ϵf2 = @. (A*E0*C*(ħω - Eg)^2)/((ħω2 - E02)^2 + C2*ħω2)/ħω*(ħω > Eg)
	return @. sqrt(ϵf1 + im*ϵf2)
end

"""

    Returns the index of refraction using the Forouhi-Bloomer model with a single term.
    Example marerials: Ta2O5, HfO2, Si3N4, Al2O3, a-Si, AlN.

        N = RI.forouhibloomer(x, λ)

            x: Array containing the five values used in the model.
               (x = [A, B, C, Egap, Ninf], where B and Eg are in eV units, and C in eV^2.)
			   Make sure that 4*C - B^2 > 0!
            λ: range of wavelengths in nm
            N: complex index of refraction

    Source: Physical Review B, 38, 1865 (1988)

"""
function forouhibloomer(x::AbstractArray{T}, λ::Array{T}) where {T<:Float64}
	A, B, C, Egap, Ninf = x
	ħω = 1240.0./λ
	Q = 0.5*sqrt(4*C - B^2)
	B0 = A/Q*(-(0.5*B^2) + (Egap*B) - Egap^2 + C)
	C0 = A/Q*(((0.5*B)*(Egap^2 + C)) - (2.0*Egap*C))
	n = @. Ninf + (B0.*ħω + C0)/(ħω^2 - B*ħω + C)
	k = @. (A*(ħω - Egap)^2)/(ħω^2 - B*ħω + C)
	return @. n + im*k
end

"""

    Returns the index of refraction using a modified Forouhi-Bloomer model with a single term.
    Example marerials: Ta2O5, HfO2, Si3N4, Al2O3, a-Si, AlN.

        N = RI.forouhibloomermodified(x, λ)

            x: Array containing the five values used in the model.
               (x = [ninf, ħω0, Γ, f, Eg], where ħω0, Γ, Eg in eV units.)
            λ: range of wavelengths in nm
            N: complex index of refraction

    Sources: Physical Review B, 34, 7018 (1986)
             Physical Review B, 38, 1865 (1988)
             Spectroscopic ellipsometry user guide (Horiba Jobin Yvon)

"""
function forouhibloomermodified(x::AbstractArray{T}, λ::Array{T}) where {T<:Float64}
	ninf, ħω0, Γ, f, Eg = x
	ħω = 1240.0 ./ λ
	n = @. ninf + (f/Γ*(Γ^2 - (ħω0 - Eg)^2)*(ħω - ħω0) + 2.0*f*Γ*(ħω0 - Eg))/((ħω - ħω0)^2 + Γ^2)
	k = @. ((f*(ħω - Eg)^2)/((ħω - ħω0)^2 + Γ^2))*(ħω > Eg)
	return @. n + im*k
end

end # module RI
