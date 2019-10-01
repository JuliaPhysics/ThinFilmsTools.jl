# Load modules
using Plots, LaTeXStrings
pyplot()
closeall()
using ThinFilmsTools

x = -10:0.01:10
len_x = length(x)

# Gaussian 1 peak
seed1 = [[0.0], [15.0, 0.0, sqrt(0.2)]]
# We generate some data using the Gaussian model
y_data = Utils.Gaussian(x, seed1) .+ randn(len_x)
sol = FitCurveModel(:Gaussian, x, y_data, seed1)
plot(PlotFitSpectrum(), x, y_data, sol.ymodel)
gui()
# Gaussian 2 peaks
seed2 = [[1.0], [15.0, -1.0, sqrt(0.2)], [10.0, 1.5, sqrt(0.1)]]
# We generate some data using the Gaussian model
y_data = Utils.Gaussian(x, seed2) .+ randn(len_x)
sol = FitCurveModel(:Gaussian, x, y_data, seed2)
plot(PlotFitSpectrum(), x, y_data, sol.ymodel)
gui()

# Lorentzian 1 peak
seed1 = [[0.0], [15.0, 0.0, 0.5]]
# We generate some data using the Lorentzian model
y_data = Utils.Lorentzian(x, seed1) .+ randn(len_x)
sol = FitCurveModel(:Lorentzian, x, y_data, seed1)
plot(PlotFitSpectrum(), x, y_data, sol.ymodel)
gui()
# Lorentzian 2 peaks
seed2 = [[0.0], [15.0, -1.0, 0.5], [30.0, 1.0, 0.5]]
# We generate some data using the Lorentzian model
y_data = Utils.Lorentzian(x, seed2) .+ randn(len_x)
sol = FitCurveModel(:Lorentzian, x, y_data, seed2)
plot(PlotFitSpectrum(), x, y_data, sol.ymodel)
gui()

# Voigtian 1 peak
seed1 = [[0.0], [18.0, 0.0, 0.0, 1.53]]
# We generate some data using the Voigtian model
y_data = Utils.Voigtian(x, seed1) .+ randn(len_x)
sol = FitCurveModel(:Voigtian, x, y_data, seed1)
plot(PlotFitSpectrum(), x, y_data, sol.ymodel)
gui()
# Voigtian 2 peaks
seed2 = [[0.0], [10.0, -1.6, 0.0, 1.53], [24.0, 5.0, 0.0, 1.53]]
# We generate some data using the Voigtian model
y_data = Utils.Voigtian(x, seed2) .+ randn(len_x)
sol = FitCurveModel(:Voigtian, x, y_data, seed2)
plot(PlotFitSpectrum(), x, y_data, sol.ymodel)
gui()
