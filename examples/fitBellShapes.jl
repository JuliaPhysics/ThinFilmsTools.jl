# Load modules
#using Plots, LaTeXStrings
#pyplot()
#using ThinFilmsTools
# using Optim

# Load modules
using Plots, LaTeXStrings
using PyCall; pygui(true)
pyplot()
closeall()
include("/home/leniac/JuliaLangDev/ThinFilmsTools/src13/ThinFilmsTools.jl")
using Main.ThinFilmsTools
using Optim

x = -10:0.01:10
len_x = length(x)

options = Optim.Options(g_abstol=1e-8, g_reltol=1e-8, iterations=10^5, show_trace=false, store_trace=true);

let
    # Gaussian 1 peak
    seed1 = [[0.0], [15.0, 0.0, sqrt(0.2)]]
    # We generate some data using the Gaussian model
    y_data = Utils.gaussian(x, seed1) .+ 0.5.*randn(len_x)
    sol = FitCurveModel(:gaussian, x, y_data, seed1)
    plot(x, [y_data, sol.ymodel], lab=["Data" "Model"])
    gui()
    # Gaussian 2 peaks
    seed2 = [[1.0], [15.0, -1.0, sqrt(0.2)], [10.0, 1.5, sqrt(0.1)]]
    # We generate some data using the Gaussian model
    y_data = Utils.gaussian(x, seed2) .+ 0.5.*randn(len_x)
    sol = FitCurveModel(:gaussian, x, y_data, seed2; options=options)
    plot(x, [y_data, sol.ymodel], lab=["Data" "Model"])
    gui()
end

let
    # Lorentzian 1 peak
    seed1 = [[0.0], [15.0, 0.0, 0.5]]
    # We generate some data using the Lorentzian model
    y_data = Utils.lorentzian(x, seed1) .+ 0.5.*randn(len_x)
    sol = FitCurveModel(:lorentzian, x, y_data, seed1; options=options)
    plot(x, [y_data, sol.ymodel], lab=["Data" "Model"])
    gui()
    # Lorentzian 2 peaks
    seed2 = [[0.0], [15.0, -1.0, 0.5], [30.0, 1.0, 0.5]]
    # We generate some data using the Lorentzian model
    y_data = Utils.lorentzian(x, seed2) .+ 0.5.*randn(len_x)
    sol = FitCurveModel(:lorentzian, x, y_data, seed2; options=options)
    plot(x, [y_data, sol.ymodel], lab=["Data" "Model"])
    gui()
end

let
    # Voigtian 1 peak
    seed1 = [[0.0], [18.0, 0.0, 0.0, 1.53]]
    # We generate some data using the Voigtian model
    y_data = Utils.voigtian(x, seed1) .+ randn(len_x)
    sol = FitCurveModel(:voigtian, x, y_data, seed1; options=options)
    plot(x, [y_data, sol.ymodel], lab=["Data" "Model"])
    gui()
    # Voigtian 2 peaks
    seed2 = [[0.0], [10.0, -1.6, 0.0, 1.53], [24.0, 5.0, 0.0, 1.53]]
    # We generate some data using the Voigtian model
    y_data = Utils.voigtian(x, seed2) .+ randn(len_x)
    sol = FitCurveModel(:voigtian, x, y_data, seed2; options=options)
    plot(x, [y_data, sol.ymodel], lab=["Data" "Model"])
    gui()
end
