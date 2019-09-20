"""

    All the recipes included in this file works with Plots.jl

"""

"""

    Recipe for producing the ThreeOmegaMethod solution plot.

        plot(TOMPlot(), solution; s=(640,480), customcolors=[black orange])
        gui()

            solution: structure from the ThreeOmegaMethod
                s: size of the figure
                customcolors: array with colors for the two lines

"""
struct TOMPlot end
@recipe function f(::TOMPlot, solution; s=(640,480), customcolors=[:black :orange])
    linecolor --> customcolors
    seriestype  :=  :path
    linestyle --> [:solid :dashdot]
    xlabel --> "Frequency [Hz]"
    ylabel --> "Temperature rise [C]"
    label --> ["Real" "Imag"]
    xscale --> :log10
    # tickfont --> font(12)
    # legendfont --> font(10)
    size --> s
    solution.f, [solution.DThreal, solution.DThimag]
end

"""

    Plot the spectrum of Reflectance, Transmittance and Absorbance as input.

        plot(TMMOPlotSpectra1D(), λ, S; s=(640,480))
        gui()

            λ: wavelength range
            S: spectrum
                s: size of the figure

"""
struct TMMOPlotSpectra1D end
@recipe function f(::TMMOPlotSpectra1D, λ, S; s=(640,480))
    # linecolor   --> customcolors
    seriestype  :=  :path
    linestyle --> :solid
    xlabel --> "Wavelength [nm]"
    # ylabel --> "Reflectance"
    # xlim --> (λ[1], λ[end])
    # ylim --> (0.0, 1.0)
    # tickfont --> font(12)
    # legendfont --> font(10)
    size --> s
    vec(λ), S
end

"""

    Plot the spectrum of Reflectance, Transmittance and Absorbance as input for angle of incidence.

        plot(TMMOPlotSpectraAngle1D(), θ, S, s=(640,480))
        gui()

            θ: angle of incidence
            S: spectrum
                s: size of the figure

"""
struct TMMOPlotSpectraAngle1D end
@recipe function f(::TMMOPlotSpectraAngle1D, θ, S; s=(640,480))
    # linecolor   --> RGB(([0,0,0]/255)...)
    seriestype  :=  :path
    linestyle --> :solid # :dash :dashdot]
    xlabel --> L"Angle of incidence [$\degree$]"
    # ylabel --> "Reflectance"
    # xlim --> (θ[1], θ[end])
    # ylim --> (0.0, 1.0)
    # tickfont --> font(12)
    # legendfont --> font(10)
    size --> s
    vec(θ), S
end

"""

    Plot the spectra input (Reflectance, Transmittance and Absorbance) for λ and θ.

        plot(TMMOPlotSpectra2D(), λ, θ, S; num_levels=80, s=(640,480))
        gui()

            λ: wavelength range
            θ: angle of incidence
            S: spectrum
                num_levels: number of levels for the contour plot
                s: size of the figure

"""
struct TMMOPlotSpectra2D end
@recipe function f(::TMMOPlotSpectra2D, λ, θ, S; num_levels=80, s=(640,480))
    seriestype  :=  :contour
    fill --> true
    levels --> num_levels
    # color --> cgrad(:viridis)
    xlabel --> "Wavelength [nm]"
    ylabel --> L"Angle of incidence [$\degree$]"
    # xlim --> (λ[1], λ[end])
    # ylim --> (θ[1], θ[end])
    # tickfont --> font(12)
    # legendfont --> font(10)
    size --> s
    vec(λ), vec(θ), Matrix(S')
end

"""

    Plot the EMF for depth and λ, at a given angle of incidence θ.

        plot(TMMOPlotEMF2D(), λ, ℓ, emf; num_levels=80, s=(640,480))
        gui()

            λ: wavelength range
            ℓ: multilayer depth
            emf: electromagnetic field
                num_levels: number of levels for the contour plot
                s: size of the figure

"""
struct TMMOPlotEMF2D end
@recipe function f(::TMMOPlotEMF2D, λ, ℓ, emf; num_levels=80, s=(640,480))
    seriestype  :=  :contour
    fill --> true
    levels --> num_levels
    # color --> cgrad(:viridis)
    xlabel --> "Wavelength [nm]"
    ylabel --> "Depth profile [nm]"
    # title --> "EMF Intensity for s-wave"
    # tickfont --> font(12)
    # size --> s
    vec(λ), vec(ℓ), Matrix(emf')
end

"""

    Plot the EMF for depth and θ, at a given λ.

        plot(TMMOPlotEMFAngle2D(), θ, ℓ, emf; num_levels=80, s=(640,480))
        gui()

            θ: angle of incidence
            ℓ: multilayer depth
            emf: electromagnetic field
                num_levels: number of levels for the contour plot
                s: size of the figure

"""
struct TMMOPlotEMFAngle2D end
@recipe function f(::TMMOPlotEMFAngle2D, θ, ℓ, emf; num_levels=80, s=(640,480))
    seriestype  :=  :contour
    fill --> true
    levels --> num_levels
    # color --> cgrad(:viridis)
    xlabel --> L"Angle of incidence [$\degree$]"
    ylabel --> "Depth profile [nm]"
    # title --> "EMF Intensity"
    # tickfont --> font(12)
    # size --> s
    vec(θ), vec(ℓ), Matrix(emf')
end

"""

    Recipe for plotting a comparison of the model and experimental spectra.

        plot(PlotFitSpectrum(), x, Xexp, Xmodel; s=(640,480))
        gui()

            x: range of variable (λ or θ)
            Xexp: experimental spectrum
            Xmodel: model spectrum
                s: size of the figure

"""
struct PlotFitSpectrum end
@recipe function f(::PlotFitSpectrum, x, Xexp, Xmodel; s=(640,480), customcolor=:black)
    linecolor --> customcolor
    seriestype := :path
    linestyle --> :solid
    @series begin
        seriestype := :scatter
        markershape --> :circle
        markersize --> 5
        markeralpha --> 0.5
        # markercolor --> RGBA(0,0,0,0.5)
        markerstrokewidth --> 0.5
        label --> "Data"
        x, Xexp
    end
    label --> "Model"
    # tickfont --> font(12)
    # legendfont --> font(10)
    size --> s
    x, Xmodel
end

"""

    Plot the objective function error obtained from SpaceSolutionRef2D and SpaceSolutionTra2D for the space solution given.

        plot(SpaceSolutionOFplot(), x1, x2, S; lims=[x1[1], x1[end], x2[1], x2[end]], num_levels=80, s=(640,480))
        gui()

            x1: range of first dimension
            x2: range of second dimension
            S: matrix with the objective function values
            Optional:
                lims: limits of the two axes
                num_levels: number of levels of the contour plot
                s: size of the figure

"""
struct SpaceSolutionOFplot end
@recipe function f(::SpaceSolutionOFplot, x1, x2, S; lims=[vec(x1)[1], vec(x1)[end], vec(x2)[1], vec(x2)[end]], num_levels=80, s=(640,480))
    seriestype := :contour
    fill --> true
    levels --> num_levels
    # color --> cgrad(:viridis)
    # tickfont --> font(12)
    size --> s
    xlims --> (lims[1], lims[2])
    ylims --> (lims[3], lims[4])
    vec(x1), vec(x2), Matrix(S')
end

"""

    Recipe for plotting a comparison of the model and experimental spectra of ellipsometry data.

        plot(PlotFitSpectrumEllip(), x, Xexp, Xmodel; s=(640,480))
        gui()

            x: range of variable (λ or θ)
            Xexp: experimental spectrum
            Xmodel: model spectrum
                s: size of the figure

"""
struct PlotFitSpectrumEllip end
@recipe function f(::PlotFitSpectrumEllip, x, Xexp, Xmodel; s=(640,480))
    linecolor --> :black
    seriestype := :path
    linestyle --> :solid
    @series begin
        seriestype := :scatter
        markershape --> [:circle :square]
        markersize --> 5
        markeralpha --> 0.5
        # markercolor --> [RGBA(0.9019607843137255, 0.6235294117647059, 0.0, 0.7) RGBA(0.33725490196078434, 0.7058823529411765, 0.9137254901960784, 0.7)]
        markerstrokewidth --> 0.5
        # markerstrokecolor --> [RGB(0.9019607843137255, 0.6235294117647059, 0.0) RGB(0.33725490196078434, 0.7058823529411765, 0.9137254901960784)]
        label --> ["Psi" "Delta"]
        x, Xexp
    end
    label --> ["Model" ""]
    # tickfont --> font(12)
    # legendfont --> font(10)
    size --> s
    x, Xmodel
end

"""

    Plots the index of refraction ath certain wavelength (usually λ0) of the multilayer structure computed from TMMOptics.

        TMMOPlotNprofile(solution; wave=:b, λ=solution.Misc.λ0, θ=solution.Beam.θ[1], s=(640,480))

            solution: structure solution from TMMO1DIsotropic
                wave = :b (both, default), :p (p-wave), :s (s-wave) of the EMF to overlap
                θ: angle for which to overlap the EMF, by default is taken the first one
                λ: wavelength for which to overlap the EMF, by default is taken the reference one
                s: size of the figure

"""
function TMMOPlotNprofile(solution; wave=:b, λ=[solution.Misc.λ0], θ=[solution.Beam.θ[1]], s=(640,480))
    # Generate colors for different layers
    assigned_cols = colorsUniqueLayers(solution.Misc.d)
    # Define the thickness vector with the incident and substrate medium
    d = solution.Misc.d[2:end-1] # remove incident and emergent media
    doffset = 0.05*sum(d)
    new_d = [-doffset; cumsum([0; d; doffset], dims=1)]
    # diffnew_d = diff(new_d)
    plot(TMMOPlot_rectangles.(diff(new_d), solution.Misc.nλ0, new_d[1:end-1], 0.0), opacity=0.6, c=assigned_cols, line=(0.0), xaxis=("Thickness profile [nm]"), legend=false, tickfont=font(12), size=s, label="", yaxis=(L"Refractive index at $\lambda_0$"))
    gui()
    if !isempty(solution.Field.emfp)
        plot(TMMOPlot_rectangles.(diff(new_d), solution.Misc.nλ0, new_d[1:end-1], 0.0), opacity=0.6, c=assigned_cols, line=(0.0), xaxis=("Thickness profile [nm]"), tickfont=font(12), size=s, label="", yaxis=(L"Refractive index at $\lambda_0$"))
        aux1 = findmin(abs.(solution.Beam.λ .- λ[1]))[2][1]
        aux2 = findmin(abs.(solution.Beam.θ .- θ[1]))[2][1]
        λaux1 = solution.Beam.λ[aux1]
        θaux2 = solution.Beam.θ[aux2]
        if wave == :p
            plot!(solution.Misc.ℓ, solution.Field.emfp[aux1, aux2, :], label="", c=RGB(([0,0,0]/255)...))
            title!(@sprintf("EMF-p resonance at λ = %0.0f nm and θ = %0.0f°", λaux1, θaux2), y=1.02)
        elseif wave == :s
            plot!(solution.Misc.ℓ, solution.Field.emfs[aux1, aux2, :], label="", c=RGB(([0,0,0]/255)...))
            title!(@sprintf("EMF-s resonance at λ = %0.0f nm and θ = %0.0f°", λaux1, θaux2), y=1.02)
        else
            plot!(solution.Misc.ℓ, solution.Field.emfp[aux1, aux2, :], label="EMF-p", c=RGB(([0,0,0]/255)...), line=(:solid))
            plot!(solution.Misc.ℓ, solution.Field.emfs[aux1, aux2, :], label="EMF-s", c=RGB(([230,159,0]/255)...), line=(:dashdot))
            title!(@sprintf("EMF resonances at λ = %0.0f nm and θ = %0.0f°", λaux1, θaux2), y=1.02)
        end
        gui()
    end
    return nothing
end

"""

    Draw rectangles for the indices of refraction profile. (IU only)

"""
TMMOPlot_rectangles(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])

"""

    Check the unique layers in the system.

"""
function colorsUniqueLayers(d)
    # Generate colors for different layers, using the thickness vector
    unique_layers = unique(d)
    # cols = distinguishable_colors((num_layers+1), [RGB(1,1,1)])[2:end]
    cols0 = distinguishable_colors(length(unique_layers)*2, [RGB(1,1,1)])[2:end]
    cols = cols0[1:2:end]
    # pcols = map(col -> (red(col), green(col), blue(col)), cols)
    # Match one color for one layer in the structure
    assigned_cols = Array{RGB}(undef, length(d))
    for i ∈ eachindex(cols)
        assigned_cols[unique_layers[i] .== d] .= cols[i]
    end
    # # reassign incident and substrate to differ from possible same layers
    assigned_cols[1] = cols0[1]
    assigned_cols[end] = cols0[end-1]
    return assigned_cols
end

"""

    Plots the photonic dispersion of the multilayer structure processed.

        TMMOPlotDispersion(solution; s=(640,480))

            solution: structure solution from TMMO1DIsotropic
                s: size of the figure

"""
function TMMOPlotDispersion(solution; s=(640,480))
    Λ = solution.Bloch.Λ # periodicity
    ω = solution.Bloch.ω .* Λ ./ 2.0 ./ π # frequency range normalized
    kblochs = solution.Bloch.κs .* Λ ./ π
    kblochp = solution.Bloch.κp .* Λ ./ π
    cc = [RGB(([0,0,0]/255)...) RGB(([230,159,0]/255)...) RGB(([86,180,233]/255)...)]
    if length(solution.Beam.θ) == 1
        tmp1 = ceil(maximum(abs.(imag.(kblochp))), digits=1) # carefull for small numbers like 0.001
        tmp2 = ceil(maximum(abs.(imag.(kblochs))), digits=1)
        p1 = plot(solution.Spectra.Rp, ω, c=cc[1], yaxis=(L"\omega\Lambda/(2\pi)"), xaxis=("Reflectance", (0.0,1.0)), label="p-wave", xticks=(0:0.25:1, string.(0.0:0.25:1)))
        plot!(p1, solution.Spectra.Rs, ω, c=cc[1], line=(:dashdot), label="s-wave")
        p2 = plot(-real.(kblochp), ω, c=cc[2], yaxis=(""), xaxis=(L"K^{r}_{Bloch}\Lambda/\pi", (-1.0,1.0)), label="", title=("p/TM     s/TE"), xticks=(-1:0.5:1, string.([1.0, 0.5, 0.0, 0.5, 1.0])), titlefontsize=12, yticks=([0.], [" "]))
        plot!(p2, real.(kblochs), ω, c=cc[2], yaxis=("", ), label="", line=(:dashdot))
        vline!([0], c=cc[1], line=(:dash), label="")
        p3 = plot(-abs.(imag.(kblochp)), ω, c=cc[3], yaxis=(""), xaxis=(L"K^{i}_{Bloch}\Lambda/\pi"), label="", title=("p/TM     s/TE"), titlefontsize=12, xticks=([-tmp1,0.,tmp2], string.([tmp1, 0.0, tmp2])), yticks=([0.], [""]), size=s)
        plot!(p3, abs.(imag.(kblochs)), ω, c=cc[3], yaxis=(""), label="", line=(:dashdot))
        vline!([0], c=cc[1], line=(:dash), label="")
        plot(p1, p2, p3, layout=(1,3), tickfont=font(12))
        gui()
    else
        # k normalized for angle-frequency dependence
        # logical matrices, used to select points which belong to the forbidden bands
        κp = cos.(kblochp .* π) # rhs equation 45 in the paper, π comes from previous normalization
        κs = cos.(kblochs .* π) # rhs equation 45 in the paper
        κp[abs.(κp) .< 1.0] .= 1.0 # propagating waves
        κp[abs.(κp) .> 1.0] .= 0.0 # evanescent waves
        κs[abs.(κs) .< 1.0] .= 1.0 # propagating waves
        κs[abs.(κs) .> 1.0] .= 0.0 # evanescent waves
        qz = sin.(deg2rad.(solution.Beam.θ)) .* π ./ 2.0 # parallel wavevector qz
        d1 = solution.Misc.d[2]; d2 = solution.Misc.d[3]
        n0 = solution.Misc.nλ0[1]; n1 = solution.Misc.nλ0[2]; n2 = solution.Misc.nλ0[3]
        ωh = Λ / π / (d1*n1 + d2*n2) * acos(-abs(n1-n2)/(n1+n2))
        ωl = Λ / π / (d2*sqrt(n2^2-n0^2) + d1*sqrt(n1^2-n0^2)) * acos(abs((n1^2*sqrt(n2^2-n0^2) - n2^2*sqrt(n1^2-n0^2))/(n1^2*sqrt(n2^2-n0^2) + n2^2*sqrt(n1^2-n0^2))))
        # choose better colormap
        clibrary(:colorbrewer)
        contourf(-qz, ω, real.(κp), levels=90, c=cgrad(:Blues, alpha=0.6), xaxis=(L"Parallel wavevector, $q_z$ (2$\pi$/$\Lambda$)"), yaxis=(L"\omega\Lambda/(2\pi)"), colorbar=false, tickfont=font(12), line=(2), title="p/TM-wave                        s/TE-wave", titlefontsize=10, grid=true, gridstyle=:dot, gridlinewidth=1.5, leg=false, size=s)
        contourf!(qz, ω, real.(κs), levels=90, c=cgrad(:Reds, alpha=0.6), line=(2))
        plot!([0.0, ω[1]], [0.0, ω[1]], line=(:solid, 1., :black, 0.7))
        plot!([0.0, -ω[1]], [0.0, ω[1]], line=(:solid, 1., :black, 0.7))
        plot!([-qz[end], qz[end]], [ωh, ωh], line=(:dash, 1, :black), lab="")
        plot!([-qz[end], qz[end]], [ωl, ωl], line=(:dash, 1, :black), lab="")
        annotate!([(qz[end]-0.05, ωh+0.03, text(L"$\omega_h$", 10)), (qz[end]-0.05, ωl-0.03, text(L"$\omega_l$",10))])
        xlims!(-qz[end], qz[end])
        ylims!(0.0, ω[1])
        gui()
    end
    return nothing
end
