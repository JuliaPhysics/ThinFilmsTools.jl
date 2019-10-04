module PlottingTools

using Colors
using RecipesBase
using LaTeXStrings
using Printf: @sprintf

export TOMPlot,
       TMMOPlotSpectra1D,
       TMMOPlotSpectraAngle1D,
       TMMOPlotSpectra2D,
       TMMOPlotEMF2D,
       TMMOPlotEMFAngle2D,
       PlotFitSpectrum,
       SpaceSolutionOFplot,
       TMMOPlotNprofile,
       TMMOPlotDispersion1D,
       TMMOPlotDispersion1Dalt,
       TMMOPlotDispersion1Dimre,
       TMMOPlotDispersion2D,
       TMMOPlotDispersion2Dalt

"""

    All the recipes included in this file works with Plots.jl

"""

"""

    Recipe for producing the ThreeOmegaMethod solution plot.

        plot(TOMPlot(), solution; s=(640,480))
        gui()

            solution: structure from the ThreeOmegaMethod
                s: size of the figure

"""
struct TOMPlot end
@recipe function f(::TOMPlot, solution; s=(640,480))
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
:black
"""
struct TMMOPlotSpectra2D end
@recipe function f(::TMMOPlotSpectra2D, λ, θ, S; num_levels=80, s=(640,480))
    seriestype  :=  :contour
    fill --> true
    levels --> num_levels
    # color --> colormap("RdBu")
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
    # color --> colormap("RdBu")
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
    # color --> colormap("RdBu")
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
@recipe function f(::PlotFitSpectrum, x, Xexp, Xmodel; s=(640,480))
    # linecolor --> customcolors[1]
    seriestype := :path
    linestyle --> :solid
    @series begin
        seriestype := :scatter
        markershape --> :circle
        markersize --> 5
        markeralpha --> 0.5
        # markercolor --> customcolors[2]
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
    # color --> colormap("RdBu")
    # tickfont --> font(12)
    size --> s
    xlims --> (lims[1], lims[2])
    ylims --> (lims[3], lims[4])
    vec(x1), vec(x2), Matrix(S')
end

"""

    Plots the index of refraction ath certain wavelength (usually λ0) of the
    multilayer structure.

        plot(TMMOPlotNprofile(), solution;
             plotemf=false, wave=:b, λ=[solution.Misc.λ0],
             θ=[solution.Beam.θ[1]], s=(640,480))
        gui()

            solution: structure solution from TMMO1DIsotropic
                wave = :b (both, default), :p (p-wave), :s (s-wave) of the EMF to overlap
                θ: angle for which to overlap the EMF, by default is taken the first one
                λ: wavelength for which to overlap the EMF, by default is taken the reference one
                s: size of the figure

"""
struct TMMOPlotNprofile end
@recipe function f(::TMMOPlotNprofile, solution; plotemf=false, wave=:b, λ=[solution.Misc.λ0], θ=[solution.Beam.θ[1]], s=(640,480))
    d = solution.Misc.d[2:end-1] # remove incident and emergent media
    doffset = 0.05*sum(d)
    new_d = [-doffset; cumsum([0; d; doffset], dims=1)]
    w = diff(new_d)
    h = solution.Misc.nλ0
    x = new_d[1:end-1]
    y = zeros(length(x))
    block1 = [x[i] .+ vcat(0.0, w[i], w[i], 0.0) for i=1:length(x)]
    block2 = [y[i] .+ vcat(0.0, 0.0, h[i], h[i]) for i=1:length(x)]
    xlabel --> "Thickness profile [nm]"
    ylabel --> L"Refractive index at $\lambda_0$"
    @series begin
        seriestype := :shape
        linewidth --> 0.0
        legend --> false
        size --> s
        label --> ""
        # Generate colors for different layers
        color --> Array(colorsUniqueLayers(solution.Misc.d)')
        alpha --> 0.6
        block1, block2
    end
    if plotemf # overlap EMF
        if ~isempty(solution.Field.emfp) # Check is not empty
            aux1 = findmin(abs.(solution.Beam.λ .- λ[1]))[2][1]
            aux2 = findmin(abs.(solution.Beam.θ .- θ[1]))[2][1]
            λaux1 = solution.Beam.λ[aux1]
            θaux2 = solution.Beam.θ[aux2]
            if isequal(wave, :b)
                @series begin
                    seriestype := :path
                    linestyle --> :solid
                    label --> "EMF-p"
                    solution.Misc.ℓ, solution.Field.emfp[aux1, aux2, :]
                end
                @series begin
                    seriestype := :path
                    linestyle --> :dash
                    label --> "EMF-s"
                    solution.Misc.ℓ, solution.Field.emfs[aux1, aux2, :]
                end
            elseif isequal(wave, :p)
                @series begin
                    seriestype := :path
                    linestyle --> :solid
                    label --> "EMF-p"
                    solution.Misc.ℓ, solution.Field.emfp[aux1, aux2, :]
                end
            elseif isequal(wave, :s)
                @series begin
                    seriestype := :path
                    linestyle --> :solid
                    label --> "EMF-s"
                    solution.Misc.ℓ, solution.Field.emfs[aux1, aux2, :]
                end
            end # if isequal(wave, :b)
        end # if ~isempty(solution.Field.emfp) # Check is not empty
    end # if plotemf
end

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

    Plots the photonic dispersion (Bloch wavevector) of a photonic crystal
    structure, computed for a range of frequencies (wavelengths) and one angle
    of incidence.

    Takes both polarisation types and renders them into one piece. You need to
    pass either the real or imaginary parts of them.

        plot(TMMOPlotDispersion1Dalt(), Bloch; kpart=:real, s=(640,480))
        gui()

            Bloch: solution.Bloch structure of the solution
                kpart: either real (:real, default) or imginary (:imag)
                s: size of the figure

"""
struct TMMOPlotDispersion1Dalt end
@recipe function f(::TMMOPlotDispersion1Dalt, Bloch; kpart=:real, s=(640,480))
    ω = Bloch.ω .* Bloch.Λ ./ 2.0 ./ π # frequency range normalized
    xlabel --> L"K_{Bloch}\Lambda/\pi"
    ylabel --> L"\omega\Lambda/(2\pi)"
    # title --> "p/TM     s/TE"
    xticks --> (-1:0.5:1, string.([1.0, 0.5, 0.0, 0.5, 1.0]))
    yticks --> ([0.], [" "])
    size --> s
    framestyle --> :zerolines
    @series begin
        seriestype := :path
        linestyle --> :solid
        label --> "p/TM"
        κ = Bloch.κp .* Bloch.Λ ./ π
        isequal(kpart, :real) ? -real.(κ) : -imag.(κ), vec(ω)
    end
    @series begin
        seriestype := :path
        linestyle --> :dashdot
        label --> "s/TE"
        κ = Bloch.κs .* Bloch.Λ ./ π
        isequal(kpart, :real) ? real.(κ) : imag.(κ), vec(ω)
    end
end

"""

    Plots the photonic dispersion (Bloch wavevector) of a photonic crystal
    structure, computed for a range of frequencies (wavelengths) and one angle
    of incidence.

    Takes one polarisation type in complex format, and plots on the left the
    imaginary part and on the right the real part.

        plot(TMMOPlotDispersion1Dimre(), Bloch; wave=:p, s=(640,480))
        gui()

            Bloch: solution.Bloch structure of the solution
                wave: either p-wave (:p, default) or s-wave (:s)
                s: size of the figure

"""
struct TMMOPlotDispersion1Dimre end
@recipe function f(::TMMOPlotDispersion1Dimre, Bloch; wave=:p, s=(640,480))
    ω = Bloch.ω .* Bloch.Λ ./ 2.0 ./ π # frequency range normalized
    xlabel --> L"K_{Bloch}\Lambda/\pi"
    ylabel --> L"\omega\Lambda/(2\pi)"
    xticks --> (-1:0.5:1, string.([1.0, 0.5, 0.0, 0.5, 1.0]))
    yticks --> ([0.], [" "])
    size --> s
    framestyle --> :zerolines
    if isequal(wave, :p)
        κ = Bloch.κp .* Bloch.Λ ./ π
    elseif isequal(wave, :s)
    else
        throw("The wave parameter must be equal to :p or :s.")
    end
    @series begin
        seriestype := :path
        linestyle --> :solid
        label --> "Imag"
        -abs.(imag.(κ)), vec(ω)
    end
    @series begin
        seriestype := :path
        linestyle --> :dashdot
        label --> "Real"
        real.(κ), vec(ω)
    end
end

"""

    Plots the photonic dispersion (Bloch wavevector) of a photonic crystal
    structure, computed for a range of frequencies (wavelengths) and one
    angle of incidence.

    Takes only one polarisation. If you want you can pass optionally the
    part of the wavevector to plot. By default, plots the real one (as the
    imaginary does not reresent much).

        plot(TMMOPlotDispersion1D(), Bloch; kpart=:real, wave=:p, s=(640,480))
        gui()

            Bloch: solution.Bloch structure of the solution
                wave: either p-wave (:p, default) or s-wave (:s)
                kpart: either real (:real, default) or imginary (:imag)
                s: size of the figure

"""
struct TMMOPlotDispersion1D end
@recipe function f(::TMMOPlotDispersion1D, Bloch; kpart=:real, wave=:p, s=(640,480))
    ω = Bloch.ω .* Bloch.Λ ./ 2.0 ./ π # frequency range normalized
    if isequal(wave, :p)
        κ = Bloch.κp .* Bloch.Λ ./ π
    elseif isequal(wave, :s)
        κ = Bloch.κs .* Bloch.Λ ./ π
    else
        throw("The wave parameter should be either :p or :s.")
    end
    seriestype := :path
    linestyle --> :solid
    xlabel --> L"K_{Bloch}\Lambda/\pi"
    ylabel --> L"\omega\Lambda/(2\pi)"
    label --> ""
    size --> s
    isequal(kpart, :real) ? real.(κ) : imag.(κ), ω
end

"""

    Plots the photonic dispersion (Bloch wavevector) of a photonic crystal
    structure, computed for a range of frequencies (wavelengths) and a range
    of angle of incidences.

        plot(TMMOPlotDispersion2Dalt(), Bloch; s=(640,480), num_levels=90)
        gui()

            Bloch: solution.Bloch structure from TMMO1DIsotropic
                s: size of the figure

"""
struct TMMOPlotDispersion2Dalt end
@recipe function f(::TMMOPlotDispersion2Dalt, Bloch; s=(640,480), num_levels=90)
    qz = Bloch.qz; ωh = Bloch.ωh; ωl = Bloch.ωl
    ω = Bloch.ω .* Bloch.Λ ./ 2.0 ./ π # frequency range normalized
    # k normalized for angle-frequency dependence
    # logical matrices, used to select points which belong to the forbidden bands
    κp = cos.(Bloch.κp .* Bloch.Λ) # rhs equation 45 in the paper, π comes from previous normalization
    κs = cos.(Bloch.κs .* Bloch.Λ) # rhs equation 45 in the paper
    κp[abs.(κp) .< 1.0] .= 1.0 # propagating p wave
    κp[abs.(κp) .> 1.0] .= 0.0 # evanescent p wave
    κs[abs.(κs) .< 1.0] .= 1.0 # propagating s wave
    κs[abs.(κs) .> 1.0] .= 0.0 # evanescent s wave
    size --> s
    xlabel --> L"Parallel wavevector, $q_z$ (2$\pi$/$\Lambda$)"
    ylabel --> L"$\omega\Lambda/(2\pi)$"
    title --> "p/TM-wave                        s/TE-wave"
    grid --> true
    # legend --> false
    gridstyle --> :dot
    gridlinewidth --> 1.5
    xlims --> (-qz[end], qz[end])
    ylims --> (0.0, ω[1])
    @series begin
        seriestype := :contour
        color --> colormap("Blues")
        # linewidth --> 2.0
        colorbar_entry --> false
        fill --> true
        levels --> num_levels
        # alpha --> 0.3
        label --> ""
        -vec(qz), vec(ω), real.(κp)
    end
    @series begin
        seriestype := :contour
        color --> colormap("Reds")
        # linewidth --> 2.0
        colorbar_entry --> false
        fill --> true
        levels --> num_levels
        # alpha --> 0.3
        label --> ""
        vec(qz), vec(ω), real.(κs)
    end
    @series begin
        seriestype := :path
        linestyle --> :solid
        linewidth --> 1.5
        linecolor --> RGBA(0.0, 0.0, 0.0, 0.7)
        label --> "Light line"
        [0.0, ω[1]], [0.0, ω[1]]
    end
    @series begin
        seriestype := :path
        linestyle --> :solid
        linewidth --> 1.5
        linecolor --> RGBA(0.0, 0.0, 0.0, 0.7)
        label --> ""
        [0.0, -ω[1]], [0.0, ω[1]]
    end
    @series begin
        seriestype := :path
        linestyle --> :dash
        linewidth --> 1.5
        linecolor --> RGBA(0.0, 0.0, 0.0, 0.7)
        label --> L"\omega_h"
        [-qz[end], qz[end]], [ωh, ωh]
    end
    @series begin
        seriestype := :path
        linestyle --> :dashdot
        linewidth --> 1.5
        linecolor --> RGBA(0.0, 0.0, 0.0, 0.7)
        label --> L"\omega_l"
        [-qz[end], qz[end]], [ωl, ωl]
    end
end

"""

    Plots the photonic dispersion (Bloch wavevector) of a photonic crystal
    structure, computed for a range of frequencies (wavelengths) and a range
    of angle of incidences.

    This function plots the Bloch wavevector for only one of the polarisation types.

        plot(TMMOPlotDispersion2D(), Bloch; wave=:p, s=(640,480), num_levels=90)
        gui()

            Bloch: solution.Bloch structure from TMMO1DIsotropic
                wave: either :p (p-wave, default) or :s (s-wave)
                s: size of the figure

"""
struct TMMOPlotDispersion2D end
@recipe function f(::TMMOPlotDispersion2D, Bloch; wave=:p, s=(640,480), num_levels=90)
    qz = Bloch.qz; ωh = Bloch.ωh; ωl = Bloch.ωl
    Λ = Bloch.Λ # periodicity
    ω = Bloch.ω .* Λ ./ 2.0 ./ π # frequency range normalized
    if isequal(wave, :p)
        kbloch = Bloch.κp .* Λ# ./ π
        cmap = colormap("Blues")
    elseif isequal(wave, :s)
        kbloch = Bloch.κs .* Λ# ./ π
        cmap = colormap("Reds")
    else
        throw("The wave type must be either :p or :s.")
    end
    # k normalized for angle-frequency dependence
    # logical matrices, used to select points which belong to the forbidden bands
    κ = cos.(kbloch)#.*π) # rhs equation 45 in the paper, π comes from previous normalization
    κ[abs.(κ) .< 1.0] .= 1.0 # propagating wave
    κ[abs.(κ) .> 1.0] .= 0.0 # evanescent wave
    size --> s
    xlabel --> L"Parallel wavevector, $q_z$ (2$\pi$/$\Lambda$)"
    ylabel --> L"\omega\Lambda/(2\pi)"
    grid --> true
    # legend --> false
    gridstyle --> :dot
    gridlinewidth --> 1.5
    xlims --> (0.0, qz[end])
    ylims --> (0.0, ω[1])
    @series begin
        seriestype := :contour
        fill --> true
        levels --> num_levels
        colorbar_entry --> false
        color --> cmap
        # alpha --> 0.9
        label --> ""
        vec(qz), vec(ω), real.(κ)
    end
    @series begin
        seriestype := :path
        linestyle --> :solid
        linewidth --> 1.5
        linecolor --> RGBA(0.0, 0.0, 0.0, 0.7)
        label --> "Light line"
        [0.0, ω[1]], [0.0, ω[1]]
    end
    @series begin
        seriestype := :path
        linestyle --> :dash
        linewidth --> 1.5
        linecolor --> RGBA(0.0, 0.0, 0.0, 0.7)
        label --> L"$\omega_h$"
        [0.0, qz[end]], [ωh, ωh]
    end
    @series begin
        seriestype := :path
        linestyle --> :dashdot
        linewidth --> 1.0
        linecolor --> RGBA(0.0, 0.0, 0.0, 0.7)
        label --> L"$\omega_l$"
        [0.0, qz[end]], [ωl, ωl]
    end
end

end # module
