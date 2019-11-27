using Plots, LaTeXStrings
pyplot()
closeall()
using ThinFilmsTools

## Sellmeier
# https://en.wikipedia.org/wiki/Sellmeier_equation
let
    lg = ["BK7" "Sapphire (OW)" "Sapphire (EW)" "Fused Silica" "MgF2"]
    plt = plot();
    λ = float.(0.2:0.001:1.6)
    B₁ = [1.03961212, 1.43134930, 1.5039759, 0.696166300, 0.48755108]
    B₂ = [0.231792344, 0.65054713, 0.55069141, 0.407942600, 0.39875031]
    B₃ = [1.01046945, 5.3414021, 6.5927379, 0.897479400, 0.897479400]
    C₁ = [6.00069867e-3, 5.2799261e-3, 5.48041129e-3, 4.67914826e-3, 0.001882178]
    C₂ = [2.00179144e-2, 1.42382647e-2, 1.47994281e-2, 1.35120631e-2, 0.008951888]
    C₃ = [103.560653, 325.017834, 402.89514, 97.9340025, 566.13559]
    for i in eachindex(B₁)
        n = RI.sellmeier([B₁[i], B₂[i], B₃[i], C₁[i], C₂[i], C₃[i]], λ)
        plot!(plt, λ.*1e3, real.(n), label=lg[i]);
    end
    xaxis!("Wavelength [nm]")
    yaxis!("Sellmeier index of refraction")
    gui()
end

## drude_lorentz
# F. Wooten, Optical properties of solids, Academic Press (1972), page 48
# https://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Lorentz_Dispersion_Model.pdf
let
    ħω = 0.6:0.01:4
    x = [
         [1.0], # ϵinf
         [1.0, 13.7, 7.22, 0.127], # f1, ħωₚ1, ħω₀1, Γ1
    ]
    n = RI.drude_lorentz(x, ħω)
    plot(ħω, [real.(n) imag.(n)], label=["n = ℜ{N}" "k = ℑ{N}"], line=([:solid :dashdot]))
    xaxis!("Energy [eV]")
    yaxis!("Drude-Lorentz AlGaN index of refraction")
    gui()
end

let
    ħω = 0.6:0.01:5
    x = [
         [1.0], # ϵinf
         [1.0, 12.7, 12, 0.1], # f1, ħωₚ1, ħω₀1, Γ1
    ]
    n = RI.drude_lorentz(x, ħω)
    plot(ħω, [real.(n) imag.(n)], label=["n = ℜ{N}" "k = ℑ{N}"], line=([:solid :dashdot]))
    xaxis!("Energy [eV]")
    yaxis!("Drude-Lorentz SiO2 index of refraction")
    gui()
end

let
    ħω = 0.6:0.01:5
    x = [
         [1.8], # ϵinf
         [1.34, 1.686, 1.686, 0.31], # f1, ħωₚ1, ħω₀1, Γ1
         [0.152, 3.541, 3.541, 0.292], # f2, ħωₚ2, ħω₀2, Γ2
    ]
    n = RI.drude_lorentz(x, ħω)
    plot(ħω, [real.(n) imag.(n)], label=["n = ℜ{N}" "k = ℑ{N}"], line=([:solid :dashdot]))
    xaxis!("Energy [eV]")
    yaxis!("Drude-Lorentz CuPc index of refraction")
    gui()
end

## tauc_lorentz
# Appl. Phys. Lett. 69 (3), 15 July 1996
let
    ħω = 0.8:0.01:5.9
    x = [
         [1.15, 1.2],
         [122, 3.45, 2.54], # A1, E01, Γ1
    ]
    n = RI.tauc_lorentz(x, ħω)
    plot(ħω, [real.(n) imag.(n)], label=["n = ℜ{N}" "k = ℑ{N}"], line=([:solid :dashdot]))
    xaxis!("Energy [eV]")
    yaxis!("Tauc-Lorentz a-Si(II) index of refraction")
    gui()
end

let
    ħω = 1.5:0.01:6
    x = [
         [3.10, 4.5],
         [59.2, 6.78, 0.49], # A1, E01, Γ1
    ]
    n = RI.tauc_lorentz(x, ħω)
    plot(ħω, [real.(n) imag.(n)], label=["n = ℜ{N}" "k = ℑ{N}"], line=([:solid :dashdot]))
    xaxis!("Energy [eV]")
    yaxis!("Tauc-Lorentz Si3N4 index of refraction")
    gui()
end

## forouhi_bloomer
# Physical Review B, 38, 1865 (1988)
let
    ħω = 0:0.01:15
    x = [
           [1.680, 2.5],
           [0.25926, 14.359, 53.747], # first oscillator
     ]
    n = RI.forouhi_bloomer(x, ħω)
    plot(ħω, [real.(n) imag.(n)], label=["n = ℜ{N}" "k = ℑ{N}"], line=([:solid :dashdot]))
    xaxis!("Energy [eV]")
    yaxis!("Forouhi-Bloomer SiC index of refraction (1 term)")
    gui()
end

let
    ħω = 0:0.01:15
    x = [
           [1.337, 2.5],
           [0.18028, 14.222, 52.148], # first oscillator
           [0.10700, 19.397, 99.605], # second oscillator
     ]
    n = RI.forouhi_bloomer(x, ħω)
    plot(ħω, [real.(n) imag.(n)], label=["n = ℜ{N}" "k = ℑ{N}"], line=([:solid :dashdot]))
    xaxis!("Energy [eV]")
    yaxis!("Forouhi-Bloomer SiC index of refraction (2 terms)")
    gui()
end

let
    ħω = 0:0.01:15
    x = [
           [1.353, 2.5],
           [0.00108, 13.227, 43.798], # first oscillator
           [0.19054, 14.447, 53.860], # second oscillator
           [0.00646, 19.335, 94.105], # third oscillator
           [0.05366, 21.940, 125.443], # fourth oscillator
     ]
    n = RI.forouhi_bloomer(x, ħω)
    plot(ħω, [real.(n) imag.(n)], label=["n = ℜ{N}" "k = ℑ{N}"], line=([:solid :dashdot]))
    xaxis!("Energy [eV]")
    yaxis!("Forouhi-Bloomer SiC index of refraction (4 terms)")
    gui()
end

let
    ħω = 0:0.01:8
    x = [
           [1.691, 0.3],
           [0.18463, 5.277, 7.504], # first oscillator
           [0.00941, 9.130, 20.934], # second oscillator
           [0.05242, 9.865, 25.172], # third oscillator
           [0.03467, 13.956, 50.062], # fourth oscillator
     ]
    n = RI.forouhi_bloomer(x, ħω)
    plot(ħω, [real.(n) imag.(n)], label=["n = ℜ{N}" "k = ℑ{N}"], line=([:solid :dashdot]))
    xaxis!("Energy [eV]")
    yaxis!("Forouhi-Bloomer InAs index of refraction (4 terms)")
    gui()
end

## cody_lorentz
# Ferlauto et al. J. Appl. Phys., Vol. 92, No. 5, 2002
let
    ħω = 1.0:0.01:5
    x = [
    	   [1.0, 1.35, 1.5, 0.7, 0.055], # ϵinf, Eg, Et, Ep, Eu
    	   [60, 3.7, 2.7], # first oscillator, A1, E01, Γ1
    ]
    n = RI.cody_lorentz(x, ħω;)# cme=:momentum)
    ϵ = n.^2
    plot(ħω, [real.(ϵ) imag.(ϵ)], label=["ϵ₁ = ℜ{ϵ}" "ϵ₂ = ℑ{ϵ}"], line=([:solid :dashdot]))
    xaxis!("Energy [eV]")
    yaxis!("Cody-Lorentz a-Si(H) dielectric function")
    gui()
end
