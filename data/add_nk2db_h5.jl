# saving the n dataset into the hdf5

using HDF5
using Plots
pyplot(reuse=false)
# using DelimitedFiles

# Workig directory
path = "/home/leniac/JuliaLangDev/ThinFilmsTools/data/"
cd(path)
file = "RefractiveIndicesDB.h5"

# # check the file
# fid = h5open(string(path, file), "r")
# names(fid)
# obj = fid["Silicon"]
# names(obj)
# # g = g_create(parent, name)
# # attrs(parent)[name] = value
# close(fid)

h5open(, "w") do file
    g = g_create(string(path,file),"silicon") # create a group
    g["lambda"] = Matrix(longda')
    g["n"] = Matrix(nSi')
    g["k"] = Matrix(kSi')
    # attrs(g)["Description"] = "Refractive index of silicon: www.refractiveindex.info" # an attribute
end

h5open(string(path,file), "cw") do file
    g = g_create(file, "silicontemperature") # create a group
    g["lambda"] = longda
    g["n20"] = nSi20
    g["k20"] = kSi20
    g["n450"] = nSi450
    g["k450"] = kSi450
    # attrs(g)["Description"] = "This group contains only a single dataset" # an attribute
end

h5open("RefractiveIndexDB.h5", "cw") do file
    g = g_create(file, "aluminum") # create a group
    g["lambda"] = Matrix(longda')
    g["n"] = Matrix(n')
    g["k"] = Matrix(k')
    # attrs(g)["Description"] = "This group contains only a single dataset" # an attribute
end

h5open("RefractiveIndexDB.h5", "cw") do file
    g = g_create(file, "bk7") # create a group
    g["lambda"] = Matrix(longda')
    g["n"] = Matrix(n')
    g["k"] = zero(n')
    # attrs(g)["Description"] = "This group contains only a single dataset" # an attribute
end

h5open("RefractiveIndexDB.h5", "cw") do file
    g = g_create(file, "chrome") # create a group
    g["lambda"] = Matrix(longda')
    g["n"] = Matrix(n')
    g["k"] = Matrix(k')
    # attrs(g)["Description"] = "This group contains only a single dataset" # an attribute
end

h5open("RefractiveIndexDB.h5", "cw") do file
    g = g_create(file, "gold") # create a group
    g["lambda"] = Matrix(longda')
    g["n"] = Matrix(n')
    g["k"] = Matrix(k')
    # attrs(g)["Description"] = "This group contains only a single dataset" # an attribute
end

h5open("RefractiveIndexDB.h5", "cw") do file
    g = g_create(file, "silver") # create a group
    g["lambda"] = Matrix(longda')
    g["n"] = Matrix(n')
    g["k"] = Matrix(k')
    # attrs(g)["Description"] = "This group contains only a single dataset" # an attribute
end

h5open("RefractiveIndexDB.h5", "cw") do file
    g = g_create(file, "sno2f") # create a group
    g["lambda"] = Matrix(longda')
    g["n"] = Matrix(n')
    g["k"] = Matrix(k')
    # attrs(g)["Description"] = "This group contains only a single dataset" # an attribute
end

h5open(string(path,file), "cw") do file
    g = g_create(file, "h2o") # create a group
    g["lambda"] = a[:,1] # saved in um!
    g["n"] = a[:,2]
    g["k"] = a[:,3]
    # attrs(g)["Description"] = "This group contains only a single dataset" # an attribute
end

readf = h5open(string(path,"RefractiveIndexDB.h5"), "r") do file
    read(file, "h2o")
end

lambda = [0.17;0.185;0.2;0.2144;0.2803;0.3021;0.365;0.4046;0.4358;0.5461;0.5876;0.5893;0.6438;0.6563;0.8621;1.083;1.395;1.7091;2.0581;3.2439].*1000
n = [1.615; 1.575; 1.550; 1.5337; 1.4940; 1.4872; 1.4745; 1.4696; 1.4666; 1.4601; 1.4585; 1.4584; 1.4567; 1.4564; 1.4525; 1.4494; 1.4458; 1.4421; 1.4372; 1.4131]
k = zeros(size(n))

h5open(string(path,file), "cw") do file
    g = g_create(file, "fusedsilicauv") # create a group
    g["lambda"] = lambda # saved in um!
    g["n"] = n
    g["k"] = k
    # attrs(g)["Description"] = "This group contains only a single dataset" # an attribute
end

readf = h5open(string(path,file), "r") do file
    read(file, "fusedsilicauv")
end



function nk()
    readf = h5open("RefractiveIndexDB.h5", "r") do file
        read(file, "gold")
    end
    return readf
end
fn = nk()
figure()
plot(fn["lambda"],fn["n"])
