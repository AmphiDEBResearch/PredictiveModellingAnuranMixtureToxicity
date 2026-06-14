using Pkg; Pkg.activate(".")

scriptsdir(args...) = joinpath(pwd(), "scripts", args...)
datadir(args...) = joinpath(pwd(),  "data", args...)
srcdir(args...) = joinpath(pwd(),  "src", args...)
libdir(args...) = joinpath(pwd(),  "lib", args...)
plotsdir(args...) = joinpath(pwd(), "plots", args...)

Pkg.develop(path = libdir("EcotoxSystems.jl"))
Pkg.develop(path = libdir("EcotoxModelFitting.jl"))
Pkg.develop(path = libdir("AmphiDEB"))

using EcotoxSystems
using EcotoxModelFitting
using AmphiDEB
