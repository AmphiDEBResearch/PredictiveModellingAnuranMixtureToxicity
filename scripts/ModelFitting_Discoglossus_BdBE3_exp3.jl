using Base.Threads

using StatsBase
using DataFrames, DataFramesMeta
using DataStructures, CSV
using Chain
using StatsPlots, Plots.Measures
using Logging

default(leg = false)
theme(:default)

using LaTeXStrings
using Suppressor
using Distributions

using EcotoxSystems, AmphiDEB, EcotoxModelFitting
import EcotoxModelFitting: Hyperdist

# re-imports from AmphiDEB
using AmphiDEB.ComponentArrays
using AmphiDEB.EcotoxSystems
using AmphiDEB.Distributions
using AmphiDEB.DataStructures

using OrdinaryDiffEq

using Revise

# import functions from other scripts
includet(scriptsdir("ModelFitting_Discoglossus_24D_UCLM.jl"))
includet(scriptsdir("ModelFitting_Discoglossus_24D_Ugent_exp2.jl"))
includet(scriptsdir("ModelFitting_Discoglossus_BdJEL423_exp2.jl"))
includet(srcdir("traits.jl"))
includet(srcdir("utils.jl"))


const PMOA_24D = "M" # only relevant to simulate combined effects of 2,4-D and Bd
const TREATMENT_MAP_BE3 = OrderedDict(
    1 => "uninfected", 
    2 => "BE3"
)

function load_data_exp3_Bd()

    return load_data_Bd(
        OrderedDict(
            :juveniles_raw => [datadir("exp_raw", "UGent", "exp3", "juveniles.csv"), 1],
            :bdloads => [datadir("exp_raw", "UGent", "exp3", "bdloads_BE3.csv"), 1],
        ),
        "BE3",
        TREATMENT_MAP_BE3
    )

end

# some things never change

plot_data_exp3_Bd() = plot_data_Bd(load_data_exp3_Bd, TREATMENT_MAP_BE3)

function define_defaultparams_exp3_Bd()::ComponentVector

    p = define_defaultparams_exp2_Bd()
    # the authors reported that all individuals reached metamorphosis within a two-day time span 
    # to reproduce this, we propagate the zoom factor to the metamorphosis threshold (different from the previous studies)

    p.spc.propagate_zoom.H_j1 = 1.

    return p
end

simulator_exp3_Bd = simulator_exp2_Bd

function setup_modelfit_exp3_Bd(pmoa::AbstractString)

    f = setup_modelfit_Bd(
        pmoa;
        data_loadfun = load_data_exp3_Bd,
        data_plotfun = plot_data_exp3_Bd,
        define_defaultparams = define_defaultparams_exp3_Bd, 
        simulator = simulator_exp3_Bd, 
        alt_priors = Dict(
            :sigma1 => Truncated(Normal(3.2, 3.2), 0, Inf)
        )
    )

    return f

end