using Pkg; Pkg.activate(".")
Pkg.add("FileIO")


using Glob 

listdir(".", "plots")

save("output.tiff", load("input.png"))