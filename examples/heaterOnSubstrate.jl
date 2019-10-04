# Gold heater on a substrate
# SI units throughout the script

# Top layer down to substrate information arrays

# Load modules
using ThinFilmsTools
using Plots
pyplot()
closeall()

### Input data
# half-width of heater line [m]
b = (12.5/2)*1e-6
# length of heater line [m]
l = 1.0e-3
# range of frequencies [Hz]
f = exp10.(LinRange(0, 9, 1500))
# power [W]
p = 0.030^2*22.11
# heater thermal resistance
ρh = [0. 0.] # I only need two actually, see F
# interface thermal resistances
# thresistances = [0. 0.]
thresistances = [0.]

# Wrap them into structures
layers = [ LayerTOM(310.0, 1.0, 0.2e-6, 2.441e6), # heater
           LayerTOM(160.0, 1.0, 525.0e-6, 2320*700.) ] # substrate]
hgeometry = HeaterGeometry(b, l, ρh)
source = Source(p, f)

# call the model
sol = ThreeOmegaMethod(layers, hgeometry, source, thresistances; int_limit=5.0e6)

plot(TOMPlot(), sol)
gui()
