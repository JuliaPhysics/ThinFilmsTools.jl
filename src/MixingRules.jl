module DF

export bruggemanspheres, looyengacylinders, looyengaspheres, lorentzlorenz, maxwellgarnettspheres, monecke, gedf, gem, Info

function Info()
    tmp1 = "Available mixing rules functions:" *
    "\n   bruggemanspheres(n1, n2, p)" *
    "\n   looyengacylinders(n1, n2, p)" *
    "\n   looyengaspheres(n1, n2, p)" *
    "\n   lorentzlorenz(n1, n2, p)" *
    "\n   maxwellgarnettspheres(n1, n2, p)" *
    "\n   monecke(n1, n2, p)" *
    "\n   gedf(n1, n2, p, β)" *
    "\n   gedf(n1, n2, ϕ, ϕc, tp)"
    return println(tmp1)
end

"""
    Calculates the complex refractive index of a porous material using Bruggeman model, for nonmagnetic and isotropic materials. It considers spheres (n1) hosted in another material (n2).

    Usage: neff = bruggemanspheres(n1,n2,p)

    Input:
       n1:    refractive index of material 1
       n2:    refractive index of material 2
       p:     porosity parameter (real number, 0<p<1), proportion of n1 over total

    Output:
       neff:  effective refractive index

    Source: PHYSICAL REVIEW B VOLUME 61, NUMBER 15 15 APRIL 2000-I.
"""
function bruggemanspheres(n1::AbstractArray{T,N}, n2::AbstractArray{T,N}, p::S) where {T<:Number, N, S<:Float64}
    # dielectric function of each media
    df1 = n1.^2
    df2 = n2.^2
    # effective dielectric function: PHYSICAL REVIEW B VOLUME 61, NUMBER 15 15 APRIL 2000-I. solved with mathematica for dfeff.
    dfeff = 0.25 .* (-df1 + 2 .* df2 + 3 .* df1.*p - 3 .* df2 .* p + sqrt.( 8 .* df1 .* df2 + (-df1 + 2 .* df2 + 3 .* df1.*p - 3 .* df2 .* p).^2 ) )
    # compute effective refractive index
    neff = dfeff.^0.5
    return neff
end # EOF bruggemanspheres(...)

"""
    Calculates the complex refractive index of the porous material using Looyenga-Landau-Lifshitz model, for nonmagnetic and isotropic materials. It considers cylinders (n1) hosted in another material (n2).

    Usage:
       neff = looyengacylinders(n1,n2,p)

    Input:
       n1:    refractive index of material 1
       n2:    refractive index of material 2
       p:     porosity parameter (real number, 0<p<1), proportion of n1 over total

    Output:
       neff:  effective refractive index

    Source: Langmuir 2013, 29, 2784−2789.
"""
function looyengacylinders(n1::AbstractArray{T,N}, n2::AbstractArray{T,N}, p::S) where {T<:Number, N, S<:Float64}
    # compute effective reractive index
    neff = (1-p) .* n2.^2  + n1.^2 .* p
end # function looyengacylinders(...)

"""
    Calculates the complex refractive index of the porous material using Looyenga-Landau-Lifshitz model, for nonmagnetic and isotropic materials. It considers spheres (n1) hosted in another material (n2).

    Usage:
       neff = looyengaspheres(n1,n2,p)

    Input:
       n1:    refractive index of material 1
       n2:    refractive index of material 2
       p:     porosity parameter (real number, 0<p<1), proportion of n1 over total

    Output:
       neff:  effective refractive index

    Source: IEEE TRANSACTIONS ON GEOSCIENCE AND REMOTE SENSING, VOL. 38, NO. 3, MAY 2000.
"""
function looyengaspheres(n1::AbstractArray{T,N}, n2::AbstractArray{T,N}, p::S) where {T<:Number, N, S<:Float64}
    # dielectric function of each media
    df1 = n1.^2
    df2 = n2.^2
    # compute effective reractive index
    neff = ( ( ( (1-p) .* (df2.^(1/3)) ) + ( (df1.^(1/3)) .* p) ).^3 ).^0.5
    return neff
end # function looyengaspheres(...)

"""
    Returns the effective dielectric function of a binary liquid mixture. p is the proportion of component n1.

    Usage:
       neff = lorentzlorenz(n1,n2,p)

    Input:
       n1:    refractive index of material 1
       n2:    refractive index of material 2
       p:     proportion of n1 over total, parameter (real number, 0<p<1)

    Output:
       neff:  effective refractive index
"""
function lorentzlorenz(n1::AbstractArray{T,N}, n2::AbstractArray{T,N}, p::S) where {T<:Number, N, S<:Float64}
    neff = sqrt.( ( n2.^2 .* ( -n1.^2 + 2 * p^2 - 2 - 2 .* p .* n1^2) ) ./ ( n1.^2 .* (p-1) - 2 - p .* n2.^2 )  )
end

"""
    Calculates the complex refractive index of the porous material using Maxwell-Garnett model, for nonmagnetic and isotropic materials. It considers spheres (n1) hosted in another material (n2).

    Usage:
       neff = maxwellgarnettspheres(n1,n2,p)

    Input:
       n1:    refractive index of material 1
       n2:    refractive index of material 2
       p:     porosity parameter (real number, 0<p<1), proportion of n1 over total

    Output:
       neff:  effective refractive index

    Source: IEEE TRANSACTIONS ON GEOSCIENCE AND REMOTE SENSING, VOL. 38, NO. 3, MAY 2000.
"""
function maxwellgarnettspheres(n1::AbstractArray{T,N}, n2::AbstractArray{T,N}, p::S) where {T<:Number, N, S<:Float64}
    # dielectric function of each media
    df1 = n1.^2
    df2 = n2.^2
    # compute effective reractive index
    nfeff = ( df1 + ( (2 .* p .* df1) .* (df2-df1) ./ ( df2 + df1 - p .* (df2-df1) ) ) ).^0.5
    return neff
end # EOF maxwellgarnettspheres(n1,n2,p)

"""
    Calculates the complex refractive index of the porous material using Monecke model, for nonmagnetic and isotropic materials.

    Usage:
        neff = monecke(n1, n2, p)

    Input:
        n1:    refractive index of material 1
        n2:    refractive index of material 2
        p:     porosity parameter (real number, 0<p<1), proportion of n1 over total

    Output:
        neff:  effective refractive index

    Source: Phys. Rev. B, Vol. 61, Num. 15, 15 April, 2000-I.
"""
function monecke(n1::AbstractArray{T,N}, n2::AbstractArray{T,N}, p::S) where {T<:Number, N, S<:Float64}
    # dielectric function of each media
    df1 = n1.^2
    df2 = n2.^2
    # compute effective dielectric function
    dfeff = ( 2 * ( p * df1 + (1-p) .* df2 ).^2 ) + ( df1 .* df2) ./ ( ( 1 + p) .* df1 ) + ( (2 - p) .* df2 )
    # compute effective reractive index
    neff = dfeff.^0.5
    return neff
end # EOF monecke(...)

"""
    Calculates the complex refractive index of a porous material using a general dielectric model, for nonmagnetic and isotropic materials. It considers n1) hosted in another material (n2).

    Usage: neff = gedf(n1, n2, p, β)

    Input:
       n1:    refractive index of filling material
       n2:    refractive index of host material
       p:     porosity parameter (real number, 0<p<1), proportion of n1 over total
       β:     either 1 or 3

    Output:
       neff:  effective refractive index

    Source: Appl. Phys. Lett. 108, 102902 (2016); doi: 10.1063/1.4943639.
"""
function gedf(n1::AbstractArray{T,N}, n2::AbstractArray{T,N}, p::S, β::S) where {T<:Number, N, S<:Number}
    # dielectric function of each media
    df1 = n1.^2
    df2 = n2.^2
    deff = df2 .* ( df1 + (β * df2) + (β * p * df1) - (β * p * df2) ) ./ ( df1 + (β * df2) - (p * df1) + (p * df2) )
    # compute effective refractive index
    neff = dfeff.^0.5
    return neff
end # EOF gedf(...)

"""
    Here we used equation (3) and the parameters ϕ, ϕc and tp as fitted. For a more detailed analysis and relation to the depolarization factors see the reference. It considers (n1) hosted in another material (n2).

    Usage: neff = gedf(n1, n2, ϕ, ϕc, tp)

    Input:
        n1:    refractive index of filling material (lower index)
        n2:    refractive index of host material (higher index)
        ϕ:     porosity parameter (real number, 0<ϕ<1), proportion of n1 over total
        ϕc:    critical volume fraction at which the high-index phase first percolates
        tp:    exponent related both to ϕc and to the shape of the grains

    Output:
        neff:  effective refractive index

    Source: Carbon 40 (2002) 2801–2815.
"""
function gem(n1::AbstractArray{T,N}, n2::AbstractArray{T,N}, ϕ::S, ϕc::S, tp::S) where {T<:Number, N, S<:Number}
    # dielectric function of each media
    α = (1-ϕc)/ϕc
    β = 1/tp
    df1 = (n1.^2).^β
    df2 = (n2.^2).^β
    deff = ((α*df1 - ϕ*df1 - α*ϕ*df1 - df2 + ϕ*df2 + α*ϕ*df2 + sqrt.(4*α*df1*df2 + ((α*(ϕ-1) + ϕ)*df1 + df2 - (1 + α)*ϕ*df2).^2))/(2*α)).^tp
    # compute effective refractive index
    neff = dfeff.^0.5
    return neff
end # EOF gem(...)

end # module DF
