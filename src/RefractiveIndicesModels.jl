module RI

import SpecialFunctions: expint

export lorentz_lorenz,
       bruggeman,
       looyenga,
       maxwell,
       monecke,
       gedf,
       gem,
       sellmeier,
       cauchy_urbach,
       drude_lorentz,
       tauc_lorentz,
       forouhi_bloomer,
       lorentz_plasmon,
       cody_lorentz,
       Info

function Info()
    tmp1 = "\n " *
    "\n " *
    "\n Available indices of refraction models:" *
    "\n " *
    "\n     bruggeman(p, n; df=1/3) -> EMA" *
    "\n     looyenga(p, n; df=1/3) -> EMA" *
    "\n     lorentz_lorenz(p, n) -> EMA" *
    "\n     maxwell(p, n; df=1/3) -> EMA" *
    "\n     monecke(p, n; dfm=:spheres, c=0.0) -> EMA" *
    "\n     gedf(x, n) -> EMA" *
    "\n     gem(x, n) -> EMA" *
    "\n     sellmeier(x, λ)" *
    "\n     cauchy_urbach(x, λ)" *
    "\n     drude_lorentz(x, ħω)" *
    "\n     tauc_lorentz(x, ħω)" *
    "\n     forouhi_bloomer(x, ħω)" *
    "\n     lorentz_plasmon(x, ħω)" *
    "\n     cody_lorentz(x, ħω; cme=:dipole)" *
    "\n " *
    "\n     To use any of these functions type: ?RI.function" *
    "\n " *
    "\n     EMA: Effective medium approximation." *
    "\n "
    return println(tmp1)
end

## Function that checks the input for some EMA.
function _checkInputs_np(p, n)
    length(p) == size(n,2) || throw("The number of proportion values should be
    equal to the number of columns of n: (length(p) == size(n,2).")
    # sum(p) == 1.0 || throw("The sum over the p array should be equal to 1.")
end

function _checkInputs_ndf(n, df)
    length(df) == size(n,2) || throw("The length of depolarisation factors should be
    equal to the number of columns of n: (length(df) == size(n,2).")
end

"""

    Returns the effective index of refraction of a mixture of two liquids.

        neff = RI.lorentz_lorenz(p, n)

            p: proportion (real number, 0<p<1) of material n[:,1].
            n: column-array with the complex indices of refraction of the materials in the
               mixture. Each column represents the index of refraction of a material for a
               range of wavelengths.

            n: effective refractive index

        Source: https://en.wikipedia.org/wiki/Clausius–Mossotti_relation
                Thin Solid Films 519 (2011) 2994-2997

"""
function lorentz_lorenz(
    p::T0, n::Matrix{T1},
) where {T0<:Real, T1<:ComplexF64}
    ϵ = n.^2
    ϵeff = @. @views (2.0*p*ϵ[:,1] + ϵ[:,1]*ϵ[:,2] - 2.0*p*ϵ[:,2] + 2.0*ϵ[:,2])/(ϵ[:,1] + p*(ϵ[:,2] - ϵ[:,1]) + 2.0)
    return sqrt.(ϵeff)
end

"""

    Calculates the complex refractive index of a porous material using Bruggeman
    model, for nonmagnetic materials.

    Supports two components due to its complexity for a general solution,
    with different types of inclusions' geometries. To input the geometrical
    inclusions choosing the optional parameter, the depolarisation factor, df.
    By default takes df = 1/3.

        neff = RI.bruggeman(p, n; df=1/3)

            p: proportion parameter (real number, 0<p<1), with the proportion of material 1.
            n: column-array with the complex indices of refraction of the materials in the
               mixture. The first column. n[:,1] has proportion p (inclusions), and n[:,2]
               proportion 1-p (host material). Each column represents the index of refraction of
               a material for a range of wavelengths.

               df: depolarisation factor (0.00001<=df<=0.99999, cannot be exactly zero or 1)
                  = 1/3 for spheres, default
                  = 1/2 for cylinders
                  = 0.99999 (close to one) for laminar structures with the fields directed parallel
                    to the direction along which layers alternate (stratification direction)
                  = 1/3<=df<=1 for prolate spheroidal shapes
                  = -1 for laminar structures for the fields directed perpendicularly
                   to the direction along which layers alternate (stratification direction)
            neff: effective index of refraction

    Source: Crystallography Reports, 2007, Vol. 52, No. 4, pp. 672-685
            Choy 2016, Effective medium theory: principles and applications
            Physical Review B, Vol 61, No 15, 2000-I, p 10437

"""
function bruggeman(
    p::T0, n::Matrix{T1};
    df=1.0/3.0,
) where {T0<:Real, T1<:ComplexF64}
    L = isequal(df, -1) ? :perpendicular : float(df)
    ϵeff = _bruggeman(float.([p, 1.0-p]), n.^2, L)
    return sqrt.(ϵeff)
end

function _bruggeman(
    p::Vector{T0}, ϵ::Matrix{T1}, df::T0,
) where {T0<:Float64, T1<:ComplexF64}
    β = @. @views ϵ[:,1]*(p[1]/df - 1.0) + ϵ[:,2]*(p[2]/df - 1.0)
    ϵeff = @. @views 0.5*(β + sqrt(β^2 + 4*ϵ[:,1]*ϵ[:,2]*(1/df - 1.0)))/(1.0/df - 1.0)
    return ϵeff
end

function _bruggeman(
    p::Vector{T0}, ϵ::Matrix{T1}, df::T2,
) where {T0<:Float64, T1<:ComplexF64, T2<:Symbol}
    ϵeff = @. @views ϵ[:,1]*p[1] + ϵ[:,2]*p[2]
    return ϵeff
end

"""

    Calculates the complex refractive index of the porous material using
    Looyenga-Landau-Lifshitz model, for nonmagnetic and isotropic materials.
    (See also Lichtenecker mixing formula: Stenzel, Optical characterization of
    thin solid films - 2018).

    Supports arbitrary number of components and the inclusions can be set
    through the polarisation factor df.

       neff = RI.looyenga(p, n; df=1.0/3.0)

           p: array with proportion parameters (real number, 0<p<1). p[i] is the proportion of
              material i respect to the total.
           n: column-array with the complex indices of refraction of the materials in the
              mixture. The first column. n[:,1] has proportion p[1], n[:,2] proportion p[2]
              and so on. Each column represents the index of refraction of a material for a
              range of wavelengths.
              df: depolarization factors
                 = 1/3 for spheres (default)
                 = 1/2 for cylinders
                 = 1 for laminar structures with the fields directed parallel to the direction
                   along which layers alternate (stratification direction)
                 = -1 for laminar structures with the fields directed perpendicular to the
                   direction along which layers alternate (stratification direction)
              The df entered as single values represent the same factors for all components, i.e.,
              there is only one geometrical type of inclusions for all of them.

           neff: effective index of refraction

    Sources: Physical Review B, Vol 61, No 15, 2000-I, p 10437
             IEEE Transactions on Geoscience and Remote Sensing, Vol. 38, 2000.

"""
function looyenga(
    p::Array{T0,N0}, n::Matrix{T1};
    df=1.0/3.0,
) where {T0<:Real, N0, T1<:ComplexF64}
    p = float.(vec(p))
    _checkInputs_np(p, n)
    ϵeff = _looyenga(p, n.^2, float.(df))
    return sqrt.(ϵeff)
end

function _looyenga(
    p::Vector{T0}, ϵ::Matrix{T1}, df::T0,
) where {T0<:Float64, T1<:ComplexF64}
    ϵeff = zeros(ComplexF64, size(ϵ,1))
    for i in collect(1:size(ϵ,2))
        @. @views ϵeff += ϵ[:,i]^(1.0 - 2.0*df)*p[i]
    end
    @. ϵeff = ϵeff^(1.0/(1.0 - 2.0*df))
    return ϵeff
end

function _looyenga(
    p::Vector{T0}, ϵ::Matrix{T1}, df::Vector{Vector{T0}},
) where {T0<:Float64, T1<:ComplexF64}
    _checkInputs_ndf(ϵ, df)
    ϵeff = zeros(ComplexF64, size(ϵ,1))
    for i in collect(1:size(ϵ,2)), j in collect(1:length(df))
        @. @views ϵeff += ϵ[:,i]^(1.0 - 2.0*df[i][j])*p[i]/length(df[i])
    end
    @. ϵeff = ϵeff^(1.0/(1.0 - 2.0*df))
    return ϵeff
end

"""

    Calculates the complex refractive index of a porous material using
    Maxwell-Garnett model, for nonmagnetic and isotropic materials.

        neff = RI.maxwell(p, n; df=1.0/3.0)

            p: array with proportion parameters (real number, 0<p<1). p[i] is the proportion of
               material i respect to the total.
            n: column-array with the complex indices of refraction of the materials in the
               mixture. The first column. n[:,1] has proportion p[1], and n[:,2] proportion p[2]
               and so on. NOTICE that the host material is taken as the last column in the array
               n. Each column represents the index of refraction of a material for a range of
               wavelengths.
               df: depolarization factors
                  = 1/3 for spheres (default)
                  = 1/2 for cylinders
                  = 1 for laminar structures with the fields directed parallel to the direction
                    along which layers alternate (stratification direction)
                  = -1 for laminar structures with the fields directed perpendicular to the
                    direction along which layers alternate (stratification direction)
               The df factors are the same for all components, i.e., there is only one
               geometrical type of inclusions for all of them.

            neff: effective index of refraction

    Source: IEEE Transactions on Geoscience and Remote Sensing, Vol. 38, 2000.
            Choy 2016, Effective medium theory: principles and applications
            Physical Review B, Vol 61, No 15, 2000-I, p 10437
            Crystallography Reports, 2007, Vol. 52, No. 4, pp. 672-685

"""
function maxwell(
    p::Array{T0,N0}, n::Matrix{T1};
    df::T0=1.0/3.0,
) where {T0<:Real, N0, T1<:ComplexF64}
    _checkInputs_np(p, n)
    ϵ = n.^2
    ϵeff = zeros(ComplexF64, size(n,1))
    for i in collect(1:size(n,2))
        @. @views ϵeff += p[i]*(ϵ[:,i] - ϵ[:,end])/(ϵ[:,end] + (ϵ[:,i] - ϵ[:,end])*df)
    end
    @. @views ϵeff = ϵ[:,end]*(ϵeff*(1.0 - df) + 1.0)/(1.0 - ϵeff*df)
    return sqrt.(ϵeff)
end

"""

    Calculates the complex refractive index of the porous material using
    Monecke model, for nonmagnetic and isotropic materials.

    Supports two components with different types of inclusions' geometries. To input the
    geometrical inclusions choosing the optional parameter dfm = :spheres or :cylinders.
    The last one are considered to be aligned infinite cylinders whose axes are perpendicular
    to the electric field. By default takes :spheres.

        neff = RI.monecke(p, n; dfm=:spheres, c=0.0)

            p: proportion parameter (real number, 0<p<1), with the proportion of n[:,1].
            n: column-array with the complex indices of refraction of the materials in the
               mixture. The first column. n[:,1] has proportion p (inclusions), and n[:,2]
               proportion 1-p (host material). Each column represents the index of refraction
               of a material for a range of wavelengths
               dfm: inclusions geometry
                   :spheres, default
                   :cylinders, aligned infinite cylinders whose axes are perpendicular to
                   the electric field
               c: parameter to define the concentration of the host see source 2, to be used
                  with the :cylinders geometry, is = 0.0 by default

            neff: effective index of refraction

    Source: Physical Review B, Vol 61, No 15, 2000-I, p 10437
            Physical Review B, Vol 55, No 11, 1997-I, p 6739

"""
function monecke(
    p::T0, n::Matrix{T1};
    dfm::T2=:spheres, c::T0=0.0,
) where {T0<:Real, T1<:ComplexF64, T2<:Symbol}
    x = float.([p c])
    ϵ = n.^2
    if isequal(dfm, :spheres)
        ϵeff = _monecke_spheres(x, ϵ)
    elseif isequal(dfm, :cylinders)
        ϵeff = _monecke_cylinders(x, ϵ)
    else
        error("The accepted df values are: :spheres or :cylinders.")
    end
    return sqrt.(ϵeff)
end

function _monecke_spheres(
    x::Matrix{T0}, ϵ::Matrix{T1},
) where {T0<:Float64, T1<:ComplexF64}
    p = x[1]
    β = @. @views (ϵ[:,2] - ϵ[:,1])/(2.0*ϵ[:,1] + ϵ[:,2])
    ϵeff = @. @views (ϵ[:,1]*β + 4.0*(1.0 - p)*ϵ[:,1] + 2.0*(1.0 - p)^2*(ϵ[:,2] - ϵ[:,1]))/(β + 1.0 - p)
    return ϵeff
end

function _monecke_cylinders(
    x::Matrix{T0}, ϵ::Matrix{T1},
) where {T0<:Float64, T1<:ComplexF64}
    p, c = x # Physical Review B, Vol 55, No 11, 1997-I, p 6739
    s = @. @views ϵ[:,1]/(ϵ[:,1] - ϵ[:,2])
    s0 = @. 0.5 - (1.0 - p)*0.5 + 1.5*(1.0 - p)^2 - (1.0 - p)^3
    η = @. 0.5*s0*p*(1.0 - p)
    ϵeff = @. @views ϵ[:,1]*(1.0 - (c - η)/s - η/(s - s0))
    return ϵeff
end

"""

    Calculates the complex refractive index of a porous material using a general
    dielectric model, for nonmagnetic and isotropic materials. It considers
    n[:,1] hosted in another material n[:,2].

        neff = RI.gedf(x, n)

            x: [p, β]
                p: porosity parameter (real number, 0<p<1), proportion of n[:,1] over total
                β: either 1 or 3
            n: column-array with the complex indices of refraction of the materials in the
               mixture. The first column. n[:,1] has proportion p (inclusions), and n[:,2]
               proportion 1-p (host material). Each column represents the index of refraction
               of a material for a range of wavelengths.

            neff: effective refractive index

    Source: Appl. Phys. Lett. 108, 102902 (2016); doi: 10.1063/1.4943639.

"""
function gedf(x::AbstractArray{S}, n::Matrix{T}) where {S<:Float64, T<:ComplexF64}
    p, β = x
    ϵ = n.^2
    ϵeff = @. @views ϵ[:,2]*(ϵ[:,1] + (β*ϵ[:,2]) + (β*p*ϵ[:,1]) - (β*p*ϵ[:,2])) / (ϵ[:,1] + (β*ϵ[:,2]) - (p*ϵ[:,1]) + (p*ϵ[:,2]))
    return sqrt.(ϵeff)
end

"""

    Here we used equation (3) and the parameters ϕ, ϕc and tp as fitted. For a
    more detailed analysis and relation to the depolarization factors see the
    reference. It considers n[:,1] hosted in another material n[:,2].

        neff = RI.gem(x, n)

            x: [ϕ, ϕc, tp]
                ϕ: porosity parameter (real number, 0<ϕ<1), proportion of n1 over total
                ϕc: critical volume fraction at which the high-index phase first percolates
                tp: exponent related both to ϕc and to the shape of the grains
            n: column-array with the complex indices of refraction of the materials in the
               mixture. The first column. n[:,1] has proportion p (inclusions), and n[:,2]
               proportion 1-p (host material). Each column represents the index of refraction
               of a material for a range of wavelengths.

            neff: effective refractive index

    Source: Carbon 40 (2002) 2801–2815.

"""
function gem(x::AbstractArray{S}, n::Matrix{T}) where {S<:Float64, T<:ComplexF64}
    ϕ, ϕc, tp = x
    α = (1-ϕc)/ϕc
    β = 1/tp
    ϵ = (n.^2).^β
    ϵeff = @. @views ((α*ϵ[:,1] - ϕ*ϵ[:,1] - α*ϕ*ϵ[:,1] - ϵ[:,2] + ϕ*ϵ[:,2] + α*ϕ*ϵ[:,2] + sqrt(4.0*α*ϵ[:,1]*ϵ[:,2] + ((α*(ϕ-1) + ϕ)*ϵ[:,1] + ϵ[:,2] - (1.0 + α)*ϕ*ϵ[:,2])^2))/(2.0*α))^tp
    return sqrt.(ϵeff)
end

"""

    Returns the index of refraction using the Sellmeier equation.

        n = RI.sellmeier(x, λ)

            x: Array containing the six values for the Sellmeier equation
               as described in https://en.wikipedia.org/wiki/Sellmeier_equation.
               (x = [B1, B2, B3, C1, C2, C3])
            λ: range of wavelengths in μm

            n: complex index of refraction

    Source: https://en.wikipedia.org/wiki/Sellmeier_equation

    Notice λ in μm.

"""
function sellmeier(x::AbstractArray{T}, λ::AbstractArray{T}) where {T<:Float64}
    λ = λ.^2
    n = @. sqrt(complex(1.0 + (x[1]*λ/(λ - x[4])) + (x[2]*λ/(λ - x[5])) + (x[3]*λ/(λ - x[6]))))
    return n
end

"""

    Returns the index of refraction using the Cauchy-Urbach model:
    n = Aₙ + Bₙ*1e4/λ^2 + Cₙ*1e9/λ^4 + im*α₀*exp((1240.0/λ - E₀)/Eᵤ).

        n = RI.cauchy_urbach(x, λ)

            x: Array containing the six values used in the model.
               (x = [
                     Aₙ, Bₙ, Cₙ, # parameters of the real part of the index of refraction
                     α₀, E₀, Eᵤ, # parameters of the imag part of the index of refraction
                ]
               )
            λ: range of wavelengths in nm

            n: complex index of refraction

    Source: E. Krous. Characterization of scandium oxide thin films for use in interference
            coatings for high-power lasers operating in the near-infrared. Master’s thesis,
            Colorado State University, 2010. Page 60.

"""
function cauchy_urbach(x::AbstractArray{T}, λ::AbstractArray{T}) where {T<:Float64}
    Aₙ, Bₙ, Cₙ, α₀, E₀, Eᵤ = x
    # n = @. Aₙ + Bₙ*1e4/λ^2 + Cₙ*1e9/λ^4 + im*Aₖ*exp(Bₖ + (1240.0/λ/Cₖ))
    # n = @. Aₙ + Bₙ*1e4/λ^2 + Cₙ*1e9/λ^4 + im*Aₖ*exp(Bₖ*(1.0/λ - 1.0/Cₖ))
    n = @. Aₙ + Bₙ*1e4/λ^2 + Cₙ*1e9/λ^4 + im*α₀*exp((1240.0/λ - E₀)/Eᵤ)
    return n
end

"""

    Returns the index of refraction using the Drude-Lorentz oscillator model.

        n = RI.drude_lorentz(x, ħω)

            x: Array of array containing the values used in the model.
               (x = [
                      [ϵinf], # dielectric function at infinite frequency (offset)
                      [f1, ħωp1, ħω01, Γ1], # first oscillator
                      [f2, ħωp2, ħω02, Γ2], # second oscillator
                      [...],  # kth oscillator
                ])
            ħω: range of energies (1240/λ(nm))

            n: complex index of refraction

    Sources: F. Wooten, Optical properties of solids, Academic Press (1972)
             http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Lorentz_Dispersion_Model.pdf
             Stenzel, Optical characterization of thin solid films - 2018

"""
function drude_lorentz(
    x::Vector{Vector{T0}}, ħω::AbstractArray{T1},
) where {T0<:Real, T1<:Real}
    ϵ = zeros(ComplexF64, length(vec(ħω)))
    ϵinf = x[1][1]
    for i = 2:length(x)
        f, ħωₚ, ħω₀, Γ = x[i]
        @. ϵ += (f*ħωₚ^2/(ħω₀^2 - ħω^2 - im*Γ*ħω))
    end
    @. ϵ += ϵinf
    return sqrt.(ϵ)
end


"""

    Returns the index of refraction using the Tauc-Lorentz dispersion model.

        n = RI.tauc_lorentz(x, ħω)

            x: Array of array containing the values used in the model.
               (x = [
                      [ϵinf, Eg],
                      [A1, E01, Γ1], # first oscillator
                      [A2, E02, Γ2], # second oscillator
                      [...],  # kth oscillator
                ])
               where A, E0 and Γ are in eV units.
            ħω: range of energies (1240/λ(nm))

            n: complex index of refraction

    Sources: Appl. Phys. Lett. 69 (3), 1996
             Appl. Phys. Lett. 69 (14), 1996
             J. Phys.: Condens. Matter 20 015216, 2008 (doi.org/10.1088/0953-8984/20/01/015216)

"""
function tauc_lorentz(
    x::Vector{Vector{T0}}, ħω::AbstractArray{T1},
) where {T0<:Real, T1<:Real}
    ϵ = zeros(ComplexF64, length(vec(ħω)))
    ϵinf, Eg = x[1]
    for i = 2:length(x)
        A, E₀, C = x[i]
        ħω2 = ħω.^2
        E02 = E₀^2
        Eg2 = Eg^2
        C2 = C^2
        γ2 = sqrt(complex(E02 - C2/2))^2
        α = sqrt(complex(4.0*E02 - C2))
        ζ4 = @. (ħω2 - γ2)^2 + α^2*C2/4
        aL = @. (Eg2 - E02)*ħω2 + Eg2*C2 - E02*(E02 + 3.0*Eg2)
        aA = @. (ħω2 - E02)*(E02 + Eg2) + Eg2*C2
        # ϵf1 = @. A*C*aL/2.0/π/ζ4/α/E₀*log((E02 + Eg2 + α*Eg)/(E02 + Eg2 - α*Eg)) -
        #        A*aA/π/ζ4/E₀*(π - atan((2.0*Eg + α)/C) + atan((α - 2.0*Eg)/C)) +
        #        2.0*A*E₀*Eg/π/ζ4/α*(ħω2 - γ2)*(π + 2.0*atan(2.0/α/C*(γ2 - Eg2))) -
        #        A*E₀*C*(ħω2 + Eg2)/π/ζ4/ħω*log(abs(ħω - Eg)/(ħω + Eg)) +
        #        2.0*A*E₀*C*Eg/π/ζ4*log((abs(ħω - Eg)*(ħω + Eg))/sqrt((E02 - Eg2)^2 + Eg2*C2))
        ϵf1 = @. A*C*aL/2.0/π/ζ4/α/E₀*log(complex((E02 + Eg2 + α*Eg)/(E02 + Eg2 - α*Eg))) -
                 A*aA/π/ζ4/E₀*(π - atan((2.0*Eg + α)/C) + atan((α - 2.0*Eg)/C)) +
                 2.0*A*E₀*Eg/π/ζ4/α*(ħω2 - γ2)*(π + 2.0*atan(2.0/α/C*(γ2 - Eg2))) -
                 A*E₀*C*(ħω2 + Eg2)/π/ζ4/ħω*log(abs(ħω - Eg)/(ħω + Eg)) +
                 2.0*A*E₀*C*Eg/π/ζ4*log(abs(ħω - Eg)*(ħω + Eg)/sqrt(complex((E02 - Eg2)^2 + Eg2*C2)))
        ϵf2 = @. (A*E₀*C*(ħω - Eg)^2)/((ħω2 - E02)^2 + C2*ħω2)/ħω*(ħω > Eg)
        @. ϵ += ϵf1 + im*ϵf2
    end
    @. ϵ += ϵinf
    return sqrt.(ϵ)
end

"""

    Returns the index of refraction using the Forouhi-Bloomer model.
    Example marerials: Ta2O5, HfO2, Si3N4, Al2O3, a-Si, AlN.

        n = RI.forouhi_bloomer(x, ħω)

        x: Array of array containing the values used in the model.
           (x = [
                  [ninf, Egap],
                  [A1, B1, C1], # first oscillator
                  [A2, B2, C2], # second oscillator
                  [...],  # kth oscillator
            ])
            where B and Egap are in eV units, and C in eV^2.
            ħω: range of energies (1240/λ(nm))

            n: complex index of refraction

    Source: Physical Review B, 34, 7018 (1986)
            Physical Review B, 38, 1865 (1988)

"""
function forouhi_bloomer(
    x::Vector{Vector{T0}}, ħω::AbstractArray{T1},
) where {T0<:Real, T1<:Real}
    N = zeros(ComplexF64, length(vec(ħω)))
    ninf, Egap = x[1]
    for i = 2:length(x)
        A, B, C = x[i]
        Q = real(0.5*sqrt(complex(4*C - B^2)))
        B₀ = A/Q*(-(0.5*B^2) + (Egap*B) - Egap^2 + C)
        C₀ = A/Q*(0.5*B*(Egap^2 + C) - (2.0*Egap*C))
        tmp1 = @. ħω^2 - B*ħω + C
        n = @. (B₀*ħω + C₀)/tmp1
        k = @. (A*(ħω - Egap)^2)/tmp1
        @. N += n + im*k
    end
    @. N += ninf
    return N
end

"""

    Returns the index of refraction using the Lorentz oscillator with plasmon damping.

        n = RI.lorentz_plasmon(x, ħω)

            x: Array of array containing the values used in the model.
                (x = [
                        [ϵinf, ħωp, Γ],
                        [ħωL1, ħωT1, γ1], # first oscillator
                        [ħωL2, ħωT2, γ2], # second oscillator
                        [...],  # kth oscillator
                ])
               where all the parameters are in eV units, and ϵinf non-dimensional.
            ħω: range of energies (1240./λ(nm))

            n: complex index of refraction

    Sources: Physical Review B, 61, 10437 (2000-I)
             J. Appl. Phys. 116, 233105 (2014)

"""
function lorentz_plasmon(
    x::Vector{Vector{T0}}, ħω::AbstractArray{T1},
) where {T0<:Real, T1<:Real}
    ϵ = zeros(ComplexF64, length(vec(ħω)))
    ϵinf, ħωp, Γ = x[1]
    for i = 2:length(x)
        ħωL, ħωT, γ = x[i]
        @. ϵ += (ħω^2 - ħωL^2 + im*γ*ħω)/(ħω^2 - ħωT^2 + im*γ*ħω)
    end
    @. ϵ -= ħωp^2/(ħω^2 + im*Γ*ħω)
    @. ϵ *= ϵinf
    return sqrt.(ϵ)
end

"""

    Returns the index of refraction using the Cody-Lorentz dispersion model.

        n = RI.cody_lorentz(x, ħω; cme=:dipole)

            x: Array of array containing the values used in the model.
               (x = [
                      [ϵinf, Eg, Et, Ep, Eu],
                      [A1, E01, Γ1], # first oscillator
                      [A2, E02, Γ2], # second oscillator
                      [...],  # kth oscillator
                ])
               where all the parameters are in eV units.
            ħω: range of energies [eV] (1240./λ(nm))
            cme: There are two versions of Cody-Lorentz as defined by Ferlauto et al.
                depending on the constant matrix element (cme): :dipole (default) or
                :momentum. If you choose cme=:dipole, then the G function is calculated
                based on the Cody constant dipole matrix element. If you set cme=:momentum
                then G is calculated based on the Tauc theory with constant momentum matrix
                element.

            n: complex index of refraction

    Sources: Ferlauto et al. J. Appl. Phys., Vol. 92, No. 5, 2002
             N. Malkova et al. Thin Solid Films 595 (2015) 32-35

    Notes:

    1. This function was not tested much since there is a scarce of biliography with data.

    2. This function uses the expint function (exponential integral function) from the
       Bridge.jl package developed by @stevenjg (and modified by @mschauer) since it is not
       yet included in the SpecialFunctions.jl module.
       For more, check:
        https://github.com/stevengj/18S096-iap17/blob/master/pset3/pset3-solutions.ipynb
        https://github.com/JuliaMath/SpecialFunctions.jl/issues/19

"""
# 1. The parameter Eu is calculated as described in Malkova et al., and the real part of the
# dielectric function is assumed to be the same as described in Ferlauto, considering that Eu is a contant.
function cody_lorentz(
    x::Vector{Vector{T0}}, ħω::AbstractArray{T1};
    cme::T2=:dipole,
) where {T0<:Real, T1<:Real, T2<:Symbol}
    isequal(cme, :dipole) && return _cody_lorentz_dipole(float.(x), vec(float.(ħω)))
    isequal(cme, :momentum) && return _cody_lorentz_momentum(float.(x), vec(float.(ħω)))
end

# Lorentz oscillator function _L and variable band edge function _G
_G(ħω::T1, Eg::T1) where {T1<:Float64} = ((ħω - Eg)/ħω)^2
_G(ħω::T1, Eg::T1, Ep::T1) where {T1<:Float64} = (ħω - Eg)^2/((ħω - Eg)^2 + Ep^2)
_L(A::T1, E0::T1, Γ::T1, ħω::T1) where {T1<:Float64} = A*E0*Γ*ħω/((ħω^2 - E0^2)^2 + (Γ*ħω)^2)
_Ld(E0::T1, Γ::T1, ħω::T1) where {T1<:Float64} = (ħω^2 - E0^2)^2 + (Γ*ħω)^2

# Determination of IcL assuming constant dipole matrix element
function _cody_lorentz_dipole(x::Vector{Vector{T1}}, ħω::Vector{T1}) where {T1<:Float64}
    ϵinf, Eg, Et, Ep, Eu = x[1]
    G = _G.(ħω, Eg, Ep)
    ϵ = zeros(ComplexF64, length(vec(ħω)))
    Iu = similar(ϵ)
    for i = 2:length(x)
        A, E0, Γ = x[i]
        # Eg<=Et<=E0 || throw("The Cody-Lorentz model requires that Eg<=Et<=E0.")
        D = (E0^2 - Et^4)/((E0^2 - Et^2)^2 + (Γ*Et)^2) + Ep^2*Et/(Et - Eg)/((Et - Eg)^2 + Ep^2)
        Eu = 0.5*Et/D
        L = _L.(A, E0, Γ, ħω)
        E1 = Et*_G.(Et, Eg, Ep)*_L.(A, E0, Γ, Et)
        _urbach_tail!.(Iu, E1, ħω, Eu, Et)
        Icl = _matrix_element_dipole.(A, E0, Eg, Et, Eu, Ep, Γ, ħω, G)
        ϵf1 = @. real(Iu) + real(Icl)
        ϵf2 = @. G*L*(ħω>Et) + E1*exp((ħω - Et)/Eu)/ħω*(ħω<=Et)
        @. ϵf2 *= ϵf2>0.0
        @. ϵ += ϵf1 + im*ϵf2
    end
    @. ϵ += ϵinf
    return sqrt.(ϵ)
end

function _matrix_element_dipole(
    A::T1, E0::T1, Eg::T1, Et::T1, Eu::T1, Ep::T1, Γ::T1, ħω::T1, G::T1,
) where {T1<:Float64}
    I0c = (π/2.0 - atan((Et - Ep)/Ep))/Ep
    c0c = 0.5*ħω*G/_Ld(E0,Γ,ħω)
    d0c = - 0.5*ħω*(ħω + Eg)^2/((ħω + Eg)^2 + Ep^2)/_Ld(E0,Γ,ħω)
    Fsq = Ep^2 + Eg^2
    ζ = E0^2 - Γ^2/2.0
    χ = sqrt(complex(4.0*E0^2 - Γ^2))
    Ksq = 2.0*Fsq + 2.0*ζ - 4.0*Eg^2
    Y4 = E0^4 + Fsq*(Ksq - Fsq) - 4.0*Eg^2*Ksq
    b0c = Y4*Fsq*(_Ld(E0,Γ,ħω)*((c0c - d0c)/ħω + 2.0*Eg*Ksq/Y4*(c0c + d0c)) - 1.0)/
          ((Ksq - Fsq)*Fsq*Y4 + E0^4*Y4 + 4.0*Eg^2*Fsq*Ksq^2)
    b1c = (2.0*Eg*Ksq*b0c - _Ld(E0,Γ,ħω)*(c0c + d0c))/Y4
    a3c = - (b1c + c0c + d0c)
    a2c = - (b0c + 2.0*Eg*b1c + ħω*(c0c - d0c))
    a1c = - (2.0*Eg*b0c - (Ksq - Fsq)*b1c + (ħω^2 - 2.0*ζ)*(c0c + d0c))
    a0c = 1.0 + (Ksq - Fsq)*b0c + 2.0*Eg*Ksq*b1c - ħω*(ħω^2 - 2.0*ζ)*(c0c - d0c)
    _Icl = 2.0*A*E0*Γ/π*(b1c*(Eg*I0c - sqrt(complex(log((Et - Eg)^2 + Ep^2)))) + b0c*I0c)
    I1t = 0.5/χ/Γ*(π - 2.0*atan(2.0*(Et^2 - ζ)/χ/Γ))
    I0at = 0.5/Γ*(π - atan((2.0*Et + χ)/Γ) + atan((-2.0*Et + χ)/Γ))
    I0bt = 0.25/χ*log((Et^2 + E0^2 + χ*Et)/(Et^2 + E0^2 - χ*Et))
    Itl = 2.0*A*E0*Γ/π*(a3c*(ζ*I1t - log(_Ld(E0,Γ,Et))^(0.25)) +
          a2c*(I0at + I0bt) + a1c*I1t + a0c*((I0at - I0bt)/E0^2) -
          c0c*log(abs(ħω - Et)) - d0c*log(ħω + Et))
    Icl = _Icl + Itl
    return Icl
end

# Determination of ItL assuming constant momentum matrix element
function _cody_lorentz_momentum(x::Vector{Vector{T1}}, ħω::Vector{T1}) where {T1<:Float64}
    ϵinf, Eg, Et, Ep, Eu = x[1]
    G = _G.(ħω, Eg)
    ϵ = zeros(ComplexF64, length(vec(ħω)))
    Iu = similar(ϵ)
    for i = 2:length(x)
        A, E0, Γ = x[i]
        # Eg<=Et<=E0 || throw("The Cody-Lorentz model requires that Eg<=Et<=E0.")
        D = (E0^2 - Et^4)/((E0^2 - Et^2)^2 + (Γ*Et)^2) + Ep^2*Et/(Et - Eg)/((Et - Eg)^2 + Ep^2)
        Eu = 0.5*Et/D
        L = _L.(A, E0, Γ, ħω)
        E1 = Et*_G.(Et, Eg)*_L.(A, E0, Γ, Et)
        _urbach_tail!.(Iu, E1, ħω, Eu, Et)
        Itl = _matrix_element_momentum.(A, E0, Eg, Et, Eu, Γ, ħω, G)
        ϵf1 = @. real(Iu) + real(Itl)
        ϵf2 = @. G*L*(ħω>Et) + E1*exp((ħω - Et)/Eu)/ħω*(ħω<=Et)
        @. ϵf2 *= ϵf2>0.0
        @. ϵ += ϵf1 + im*ϵf2
    end
    @. ϵ += ϵinf
    return sqrt.(ϵ)
end

function _matrix_element_momentum(
    A::T1, E0::T1, Eg::T1, Et::T1, Eu::T1, Γ::T1, ħω::T1, G::T1,
) where {T1<:Float64}
    ζ = E0^2 - Γ^2/2.0 # I define ζ as ζ squared by removing the sqrt
    χ = sqrt(complex(4.0*E0^2 - Γ^2))
    c0t = 0.5*ħω*G/_Ld(E0,Γ,ħω)
    d0t = - 0.5*(ħω + Eg)^2/ħω/_Ld(E0,Γ,ħω)
    a3t = - (c0t + d0t)
    a2t = - ħω*(c0t - d0t)
    a0t = 1.0 + (ħω^2 - 2.0*ζ)*a2t
    a1t = (ħω^2 - 2.0*ζ)*a3t
    I1t = 0.5/χ/Γ*(π - 2.0*atan(2.0*(Et^2 - ζ)/χ/Γ))
    I0at = 0.5/Γ*(π - atan((2.0*Et + χ)/Γ) + atan((-2.0*Et + χ)/Γ))
    I0bt = 0.25/χ*log((Et^2 + E0^2 + χ*Et)/(Et^2 + E0^2 - χ*Et))
    Itl = 2.0*A*E0*Γ/π*(a3t*(ζ*I1t - log(_Ld(E0,Γ,Et))^(0.25)) +
          a2t*(I0at + I0bt) + a1t*I1t + a0t*((I0at - I0bt)/E0^2) -
          c0t*log(abs(ħω - Et)) - d0t*log(ħω + Et))
    return Itl
end

function _urbach_tail!(
        Iu::T0, E1::T1, ħω::T1, Eu::T1, Et::T1,
) where {T0<:ComplexF64, T1<:Float64}
    Iu = E1/π/ħω*(
            exp((ħω - Et)/Eu)*(expint(complex((Et - ħω)/Eu)) - expint(complex(-ħω/Eu))) -
            exp(-(ħω + Et)/Eu)*(expint(complex((Et + ħω)/Eu)) - expint(complex(ħω/Eu)))
    )
    return nothing
end

end # module RI
