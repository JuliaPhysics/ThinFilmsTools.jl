# Gold heater on a substrate
# SI units throughout the script

# Top layer down to substrate information arrays

# Load modules
using Plots, LaTeXStrings
pyplot()
using ThinFilmsTools

function main()
    # Half-width of heater line [m]
    b = (12.5/2)*1e-6
    # Length of heater line [m]
    l = 1.0e-3
    # Range of frequencies [Hz]
    f = exp10.(LinRange(0, 9, 1500))
    # Power [W]
    p = 0.030^2*22.11
    # Heater thermal resistance
    ρh = [0. 0.]
    # Interface thermal resistances
    thresistances = [1.0e-8 1.5e-9 1.0e-12]
    # Wrap them into structures
    layers = [
        LayerTOM(310.0, 1.0, 0.2e-6, 2.441e6), # heater
        LayerTOM(1.0, 1.0, 1.0e-6, 2320*700.0), # specimen_1
        LayerTOM(0.1, 1.0, 1.0e-6, 2320*700*0.1), # specimen_2
        LayerTOM(160.0, 1.0, 525.0e-6, 2320*700.), # substrate
    ]
    hgeometry = HeaterGeometry(b, l, ρh)
    source = Source(p, f)
    # Call the model
    sol = three_omega(layers, hgeometry, source, thresistances)
    return sol
end

sol = main()

plot(TOMPlot(), sol)
gui()
