module RIdb

using Interpolations
using HDF5

export aluminum, air, bk7, chrome, dummy, glass, gold, silicon, silicontemperature, silver, sno2f, h2o, etoh, Info

function Info()
    tmp1 = "Available functions for materials index of refraction:" *
        "\n   aluminum(λ), λ ∈ [4.15, 31000] (nm)" *
        "\n   air(λ), λ ∈ any (nm)" *
        "\n   bk7(λ), λ ∈ [191, 1239] (nm)" *
        "\n   chrome(λ), λ ∈ [207, 1240] (nm)" *
        "\n   dummy(λ), λ ∈ any (nm)" *
        "\n   glass(λ), λ ∈ [300, 2500] (nm)" *
        "\n   gold(λ), λ ∈ [34.15, 10240] (nm)" *
        "\n   silicon(λ), λ ∈ [163.15, 25000] (nm)" *
        "\n   silicontemperature(λ, T), λ ∈ [264, 826.5], T ∈ [20, 450]" *
        "\n   silver(λ), λ ∈ [0.124, 9919] (nm)" *
        "\n   sno2f(λ), λ ∈ [308.25, 2490.9] (nm), fluor doped!" *
        "\n   h2o(λ), λ ∈ [10.0, 1e10] (nm)" *
        "\n   etoh(λ), λ ∈ [476.5, 830] (nm)"
    return println(tmp1)
end

function read_RIfiles_DB(material)
    datapath = joinpath(@__DIR__, "..", "data/")
    readf = h5open(joinpath(datapath, "RefractiveIndicesDB.h5"), "r") do file
        read(file, material)
    end
    readf
end

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
    Returns the index of refraction of glass in complex format, for a given range of wavelengths in nm, using the Sellmeier equation.
        N = glass(λ)
    Input args:
        λ = wavelength range (nm)

    Source: http://refractiveindex.info/
"""
function glass(λ)
    λsq = λ.^2
    N = complex.(sqrt.(1.0 .+ (1.03961212 .* λsq ./ (λsq .- 0.00600069867)) .+ (0.231792344 .* λsq ./ (λsq .- 0.0200179144)) .+ (1.01046945 .* λsq ./ (λsq .- 103.560653))))
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
    readf = read_RIfiles_DB("aluminum")
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
    readf = read_RIfiles_DB("bk7")
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
    readf = read_RIfiles_DB("chrome")
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
    readf = read_RIfiles_DB("gold")
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
    readf = read_RIfiles_DB("silicon")
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
        readf = read_RIfiles_DB("silicontemperature")
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
    readf = read_RIfiles_DB("silver")
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
    readf = read_RIfiles_DB("sno2f")
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
    readf = read_RIfiles_DB("h2o")
    knots = (sort(vec(readf["lambda"])*1000),)
    spl_n = interpolate(knots, vec(readf["n"]), Gridded(Linear()))
    spl_k = interpolate(knots, vec(readf["k"]), Gridded(Linear()))
    N::Array{ComplexF64} = spl_n(vec(λ)) + im*spl_k(vec(λ))
    return N
end # EOF h2o()

end # module RIdb
