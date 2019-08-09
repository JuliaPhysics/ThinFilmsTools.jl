module RIdb

using Interpolations, HDF5

export aluminum, air, bk7, chrome, dummy, glass, gold, silicon, silicontemperature, silver, sno2f, h2o, etoh

# Workig directory
file = "RefractiveIndexDB.h5"

function air(λ::AbstractArray{T,M}) where {T<:Number, M}
    N::Array{ComplexF64} = 1.00029 * ones.(lastindex(λ)) + im .* zeros.(lastindex(λ))
end

function dummy(λ::AbstractArray{T,N}, x::S, y::S) where {T<:Number, N, S<:Number}
    Nc::Array{ComplexF64} = ones.(length(λ)) .* (x + im*y)
end # EOF dummy(...)

#http://www.cvilaser.com/PublicPages/Pages/TechnicalTips.aspx
function glass(λ::AbstractArray{T,M}) where {T<:Number, M}
    N::Array{ComplexF64} = ones.(lastindex(λ)) .* (1.52 + im * 0)
end # EOF glass(...)

#source: http://refractiveindex.info/
function etoh(λ)
    λ = λ * 1e6
    C1 = 1.35265
    C2 = 0.00306
    C3 = 0.00002
    N::Array{ComplexF64} = C1 + C2./(λ.^2) + C3./(λ.^4) + im * zero(λ)
end

#source: http://www-swiss.ai.mit.edu/~jaffer/FreeSnell/nk.html
function aluminum(λ::AbstractArray{T,M}) where {T<:Number, M}
    readf = h5open(file, "r") do file
        read(file, "aluminum")
    end

    knots = (sort(vec(readf["lambda"])).*1e9,)
    spl_n = interpolate(knots, vec(readf["n"]), Gridded(Linear()))
    spl_k = interpolate(knots, vec(readf["k"]), Gridded(Linear()))
    N::Array{ComplexF64} = spl_n(vec(λ)) + im*spl_k(vec(λ))
    return N
end # EOF aluminum(...)

function bk7(λ::AbstractArray{T,M}) where {T<:Number, M}
    readf = h5open(file, "r") do file
        read(file, "bk7")
    end
    knots = (sort(vec(readf["lambda"]).*1e9),)
    spl_n = interpolate(knots, vec(readf["n"]), Gridded(Linear()))
    spl_k = interpolate(knots, vec(readf["k"]), Gridded(Linear()))
    N::Array{ComplexF64} = spl_n(vec(λ)) + im*spl_k(λ)
    return N
end # EOF bk7(...)

#source: http://www-swiss.ai.mit.edu/~jaffer/FreeSnell/nk.html
function chrome(λ::AbstractArray{T,M}) where {T<:Number, M}
    readf = h5open(file, "r") do file
        read(file, "chrome")
    end
    knots = (sort(vec(readf["lambda"])),)
    spl_n = interpolate(knots, vec(readf["n"]), Gridded(Linear()))
    spl_k = interpolate(knots, vec(readf["k"]), Gridded(Linear()))
    N::Array{ComplexF64} = spl_n(vec(λ)) + im*spl_k(vec(λ))
    return N
end # EOF chrome(...)

#source: people.csail.mit.edu/jaffer/FreeSnell/nk.html#nk
function gold(λ::AbstractArray{T,M}) where {T<:Number, M}
    readf = h5open(file, "r") do file
        read(file, "gold")
    end
    knots = (sort(vec(readf["lambda"])),)
    spl_n = interpolate(knots, vec(readf["n"]), Gridded(Linear()))
    spl_k = interpolate(knots, vec(readf["k"]), Gridded(Linear()))
    N::Array{ComplexF64} = spl_n(vec(λ)) + im*spl_k(vec(λ))
    return N
end # EOF gold(...)

#source:
function silicon(λ::AbstractArray{T,M}) where {T<:Number, M}
    readf = h5open(file, "r") do file
        read(file, "silicon")
    end
    knots = (sort(vec(readf["lambda"]))*1e9,)
    spl_n = interpolate(knots, vec(readf["n"]), Gridded(Linear()))
    spl_k = interpolate(knots, vec(readf["k"]), Gridded(Linear()))
    N::Array{ComplexF64} = spl_n(vec(λ)) + im*spl_k(vec(λ))
    return N
end # EOF silicon(...)

#source: http://refractiveindex.info/, temperature in C (valid range: 20-450)
function silicontemperature(λ::AbstractArray{T,M}, t::S) where {T<:Number, M, S<:Number}

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
end # EOF silicontemperature(...)

#source:
function silver(λ::AbstractArray{T,M}) where {T<:Number, M}
    readf = h5open(file, "r") do file
        read(file, "silver")
    end
    knots = (sort(vec(readf["lambda"])),)
    spl_n = interpolate(knots, vec(readf["n"]), Gridded(Linear()))
    spl_k = interpolate(knots, vec(readf["k"]), Gridded(Linear()))
    N::Array{ComplexF64} = spl_n(vec(λ)) + im*spl_k(vec(λ))
    return N
end # EOF silver(...)

#source: Mater. Res. Soc. Symp. Proc., 426, (1996) 449
function sno2f(λ::AbstractArray{T,M}) where {T<:Number, M}
    readf = h5open(file, "r") do file
        read(file, "sno2f")
    end
    knots = (sort(vec(readf["lambda"])*1e9),)
    spl_n = interpolate(knots, vec(readf["n"]), Gridded(Linear()))
    spl_k = interpolate(knots, vec(readf["k"]), Gridded(Linear()))
    N::Array{ComplexF64} = spl_n(vec(λ)) + im*spl_k(vec(λ))
    return N
end # EOF sno2f(...)

function h2o(λ::AbstractArray{T,M}) where {T<:Number, M}
    readf = h5open(file, "r") do file
        read(file, "h2o")
    end
    knots = (sort(vec(readf["lambda"])*1000),)
    spl_n = interpolate(knots, vec(readf["n"]), Gridded(Linear()))
    spl_k = interpolate(knots, vec(readf["k"]), Gridded(Linear()))
    N::Array{ComplexF64} = spl_n(vec(λ)) + im*spl_k(vec(λ))
    return N
end # EOF h2o(...)

end # module RIdb
