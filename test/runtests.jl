ENV["MPLBACKEND"]="agg" # no GUI

using ThinFilmsTools

include("test3omega.jl")

include("testTMMO.jl")
