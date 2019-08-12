module RIdb

using Interpolations
using HDF5

export aluminum, air, bk7, chrome, dummy, glass, gold, silicon, silicontemperature, silver, sno2f, h2o, etoh

file = "RefractiveIndicesDB.h5"

"""
    Returns the index of refraction of air in complex format, for a given range of wavelengths in nm.
        N = air(λ)
    Input args:
        λ = wavelength range (nm)
"""
function air(λ)
    N::Array{ComplexF64} = 1.00029 * ones.(length(λ)) + im .* zeros.(length(λ))
end # EOF air()

"""
    Returns the index of refraction for constants values in complex format, for a given range of wavelengths in nm.
        N = dummy(λ, a, b)
    Input args:
        λ = wavelength range (nm)
        a = real part constant value
        b = imaginary part constant value
"""
function dummy(λ, x::S, y::S) where {S<:Number}
    N::Array{ComplexF64} = ones.(length(λ)) .* (x + im*y)
end # EOF dummy()

"""
    Returns the index of refraction of glass in complex format, for a given range of wavelengths in nm.
        N = glass(λ)
    Input args:
        λ = wavelength range (nm)

    Source: http://www.cvilaser.com/PublicPages/Pages/TechnicalTips.aspx
"""
function glass(λ)
    N::Array{ComplexF64} = ones.(length(λ)) .* (1.52 + im * 0.0)
end # EOF glass()

"""
    Returns the index of refraction of liquid ethanol in complex format, for a given range of wavelengths in nm.
        N = etoh(λ)
    Input args:
        λ = wavelength range (nm)

    Source: http://refractiveindex.info/
"""
function etoh(λ)
    λ *= 1e6
    c = [1.35265; 0.00306; 0.00002]
    N::Array{ComplexF64} = c[1] .+ c[2]./(λ.^2) .+ c[3]./(λ.^4) .+ im * zeros(λ)
end # EOF etho()

"""
    Returns the index of refraction of aluminum in complex format, for a given range of wavelengths in nm.
        N = aluminum(λ)
    Input args:
        λ = wavelength range (nm)

    Source: http://www-swiss.ai.mit.edu/~jaffer/FreeSnell/nk.html
"""
function aluminum(λ)
    readf = h5open(file, "r") do file
        read(file, "aluminum")
    end
    knots = (sort(vec(readf["lambda"])).*1e9,)
    spl_n = interpolate(knots, vec(readf["n"]), Gridded(Linear()))
    spl_k = interpolate(knots, vec(readf["k"]), Gridded(Linear()))
    N::Array{ComplexF64} = spl_n(vec(λ)) + im*spl_k(vec(λ))
    return N
end # EOF aluminum()

"""
    Returns the index of refraction of BK7 in complex format, for a given range of wavelengths in nm.
        N = bk7(λ)
    Input args:
        λ = wavelength range (nm)

    Source: http://www-swiss.ai.mit.edu/~jaffer/FreeSnell/nk.html
"""
function bk7(λ)
    readf = h5open(file, "r") do file
        read(file, "bk7")
    end
    knots = (sort(vec(readf["lambda"]).*1e9),)
    spl_n = interpolate(knots, vec(readf["n"]), Gridded(Linear()))
    spl_k = interpolate(knots, vec(readf["k"]), Gridded(Linear()))
    N::Array{ComplexF64} = spl_n(vec(λ)) + im*spl_k(λ)
    return N
end # EOF bk7()

"""
    Returns the index of refraction of chrome in complex format, for a given range of wavelengths in nm.
        N = chrome(λ)
    Input args:
        λ = wavelength range (nm)

    Source: http://www-swiss.ai.mit.edu/~jaffer/FreeSnell/nk.html
"""
function chrome(λ)
    readf = h5open(file, "r") do file
        read(file, "chrome")
    end
    knots = (sort(vec(readf["lambda"])),)
    spl_n = interpolate(knots, vec(readf["n"]), Gridded(Linear()))
    spl_k = interpolate(knots, vec(readf["k"]), Gridded(Linear()))
    N::Array{ComplexF64} = spl_n(vec(λ)) + im*spl_k(vec(λ))
    return N
end # EOF chrome()

"""
    Returns the index of refraction of gold in complex format, for a given range of wavelengths in nm.
        N = gold(λ)
    Input args:
        λ = wavelength range (nm)

    Source: http://www-swiss.ai.mit.edu/~jaffer/FreeSnell/nk.html
"""
function gold(λ)
    readf = h5open(file, "r") do file
        read(file, "gold")
    end
    knots = (sort(vec(readf["lambda"])),)
    spl_n = interpolate(knots, vec(readf["n"]), Gridded(Linear()))
    spl_k = interpolate(knots, vec(readf["k"]), Gridded(Linear()))
    N::Array{ComplexF64} = spl_n(vec(λ)) + im*spl_k(vec(λ))
    return N
end # EOF gold()

"""
    Returns the index of refraction of crystalline silicon in complex format, for a given range of wavelengths in nm.
        N = silicon(λ)
    Input args:
        λ = wavelength range (nm)

    Source: http://www-swiss.ai.mit.edu/~jaffer/FreeSnell/nk.html
"""
function silicon(λ)
    readf = h5open(file, "r") do file
        read(file, "silicon")
    end
    knots = (sort(vec(readf["lambda"]))*1e9,)
    spl_n = interpolate(knots, vec(readf["n"]), Gridded(Linear()))
    spl_k = interpolate(knots, vec(readf["k"]), Gridded(Linear()))
    N::Array{ComplexF64} = spl_n(vec(λ)) + im*spl_k(vec(λ))
    return N
end # EOF silicon(...)

"""
    Returns the index of refraction of crystalline silicon in complex format, for a given range of wavelengths in nm and a one temperature value.
        N = silicontemperature(λ, T)
    Input args:
        λ = wavelength range (nm)
        T = value of temperature (C)
    Note: temperature in C (valid range: 20-450)

    Source: http://refractiveindex.info/
"""
function silicontemperature(λ, t::S) where {S<:Number}

    if ( t > 450) | (t < 20 )
        error("Temperature range is invalid (20-450 C).")
    else
        readf = h5open(file, "r") do file
            read(file, "silicontemperature")
        end
        # λ = vec(λ)
        lambda = sort(vec(readf["lambda"]))
        # https://discourse.julialang.org/t/a-question-on-extrapolation-with-interpolations-jl/3669/4
        lambda = [prevfloat(lambda[1]);lambda[2:end-1];nextfloat(lambda[end])]
        knots = (lambda,)
        spl_n20 = interpolate(knots, vec(readf["n20"]), Gridded(Linear()))
        aux_n20 = spl_n20(λ)
        spl_k20 = interpolate(knots, vec(readf["k20"]), Gridded(Linear()))
        aux_k20 = spl_k20(λ)
        spl_n450 = interpolate(knots, vec(readf["n450"]), Gridded(Linear()))
        aux_n450 = spl_n450(λ)
        spl_k450 = interpolate(knots, vec(readf["k450"]), Gridded(Linear()))
        aux_k450 = spl_k450(λ)
        dndt = (aux_n450-aux_n20)/(450-20)
        dkdt = (aux_k450-aux_k20)/(450-20)
        if t < 215
            spl_n = interpolate(knots, vec(readf["n20"]), Gridded(Linear()))
            aux = spl_n(λ).*(1 .+ dndt*(t-20.))
            spl_k = interpolate(knots, vec(readf["k20"]), Gridded(Linear()))
            aux2 = spl_k(λ).*(1 .+ dkdt*(t-20.))
        else
            spl_n = interpolate(knots, vec(readf["n450"]), Gridded(Linear()))
            aux = spl_n(λ).*(1 .+ dndt*(t-20.))
            spl_k = interpolate(knots, vec(readf["k450"]), Gridded(Linear()))
            aux2 = spl_k(λ).*(1 .+ dkdt*(t-20.))
        end
        N::Array{ComplexF64} = aux + im*aux2
        return N
    end
end # EOF silicontemperature()

"""
    Returns the index of refraction of silver in complex format, for a given range of wavelengths in nm.
        N = silver(λ)
    Input args:
        λ = wavelength range (nm)

    Source: http://www-swiss.ai.mit.edu/~jaffer/FreeSnell/nk.html
"""
function silver(λ)
    readf = h5open(file, "r") do file
        read(file, "silver")
    end
    knots = (sort(vec(readf["lambda"])),)
    spl_n = interpolate(knots, vec(readf["n"]), Gridded(Linear()))
    spl_k = interpolate(knots, vec(readf["k"]), Gridded(Linear()))
    N::Array{ComplexF64} = spl_n(vec(λ)) + im*spl_k(vec(λ))
    return N
end # EOF silver()

"""
    Returns the index of refraction of SnO2:F in complex format, for a given range of wavelengths in nm.
        N = sno2f(λ)
    Input args:
        λ = wavelength range (nm)

    Source: Mater. Res. Soc. Symp. Proc., 426, (1996) 449
"""
function sno2f(λ::AbstractArray{T,M}) where {T<:Number, M}
    readf = h5open(file, "r") do file
        read(file, "sno2f")
    end
    knots = (sort(vec(readf["lambda"])*1e9),)
    spl_n = interpolate(knots, vec(readf["n"]), Gridded(Linear()))
    spl_k = interpolate(knots, vec(readf["k"]), Gridded(Linear()))
    N::Array{ComplexF64} = spl_n(vec(λ)) + im*spl_k(vec(λ))
    return N
end # EOF sno2f()

"""
    Returns the index of refraction of liquid water in complex format, for a given range of wavelengths in nm.
        N = h2o(λ)
    Input args:
        λ = wavelength range (nm)

    Source: http://refractiveindex.info/
"""
function h2o(λ)
    readf = h5open(file, "r") do file
        read(file, "h2o")
    end
    knots = (sort(vec(readf["lambda"])*1000),)
    spl_n = interpolate(knots, vec(readf["n"]), Gridded(Linear()))
    spl_k = interpolate(knots, vec(readf["k"]), Gridded(Linear()))
    N::Array{ComplexF64} = spl_n(vec(λ)) + im*spl_k(vec(λ))
    return N
end # EOF h2o()

end # module RIdb
