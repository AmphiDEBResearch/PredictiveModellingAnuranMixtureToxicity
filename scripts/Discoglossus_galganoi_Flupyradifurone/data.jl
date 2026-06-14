# functions used to fit Flupyradifurone models
# the functions defined here perform the model fitting for larval and juvenile life stages

# packages

using Base.Threads

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
import AmphiDEB: ComponentVector


# source files

using Revise
includet(srcdir("traits.jl"))
includet(srcdir("utils.jl"))

# constants

const EGG_DRYMASS_MG = 1 # dry mass of an egg, if not fitted
const WETMASS_AT_REPRO_MEASURE_MG = 16530 # wet mass of adult female at time of clutch size measurement (not including egg mass)
const REPRODUCTION_PERIOD = 365 # reproduction period in days
const PMOAS = ["G", "M", "A", "R", "Hneg", "Hpos", "KAP"]


const TREATMENT_IDS = Dict(
    0. => 1,
    10. => 2,
    100. => 3
)

# functions

function load_data(;
    paths::OrderedDict = OrderedDict(
        :aquatic => [datadir("exp_raw", "UCLM", "discoglossus_galganoi", "disco_flp_24D_aquatic.csv"), 1], # number indicates row where data header is located (omitting metadata)
        :metamorphs => [datadir("exp_raw", "UCLM", "discoglossus_galganoi", "disco_flp_24D_metamorphs.csv"), 1],
        ##:adults => [datadir("exp_raw", "Discoglossus_03_adults.csv"), 1]
    ))

    data = OrderedDict()

    for (key,info) in pairs(paths)
        path, header = info
        data[key] = CSV.read(path, DataFrame, header = header)
    end

    dropmissing!(data[:aquatic])

    data[:aquatic] = @subset(
        data[:aquatic], 
        :D_ppm .== 0, # omit 2,4D treatments
        #:D_ppm .< maximum(:D_ppm) # omit the highest treatment due to 100% mortality
        ) 

    data[:aquatic].num_tadpoles = float.(data[:aquatic].num_tadpoles)
    data[:aquatic].fract_tadpoles = data[:aquatic].num_tadpoles ./ data[:aquatic].survival
    rename!(data[:aquatic], :F_ppm => :C_W_1) 
    data[:aquatic] = EcotoxSystems.relative_response(
        data[:aquatic], 
        [:wetmass_mg, :num_tadpoles], 
        :C_W_1; 
        groupby_vars = [:t_exp, :aquarium]
    )

    # process metamorph data

    dropmissing!(data[:metamorphs])

    data[:metamorphs].wetmass_G42_mg = float.(data[:metamorphs].wetmass_G42_mg)
    data[:metamorphs].wetmass_G46_mg = float.(data[:metamorphs].wetmass_G46_mg)

    data[:metamorphs] = @subset(
        data[:metamorphs], 
        :D_ppm .== 0, # omit Flupyradifurone treatments
        ) 

    data[:metamorphs][!,:treatment_id] = [TREATMENT_IDS[t] for t in data[:metamorphs].F_ppm]

    rename!(data[:metamorphs], :F_ppm => :C_W_1)


    return data
end

function plot_data(;kwargs...)

    plt_aqua = @df f.data[:aquatic] plot(
        groupedlineplot(
            :t_exp, :wetmass_mg, :C_W_1, 
            layout = (1,length(unique(:C_W_1))), title = hcat(["$x mg/L" for x in unique(:C_W_1)]...), 
            xlim = (-5, 30), 
            ylim = (100,350), fillalpha = .2, lw = 2, color = :black, marker = true, ylabel = ["Wetmass (mg)" "" ""],
            leg = [:bottomright false false], label = "Observed mean (P₅-P₉₅)"
        ), 
        groupedlineplot(
            :t_exp, :fract_tadpoles, :C_W_1, layout = (1,length(unique(:C_W_1))), leg = false, 
            fillalpha = .2, lw = 2, color = :black, marker = true, ylabel = ["Fraction of \n tadpoles" "" ""],
            xlim = (-5, 30), xlabel = ["" "Time since start of experiment (d)" ""],  
        ), 
        layout = (2,1), size = (1000,600), leftmargin = 7.5mm
    )

    c = 0
    num_concs = length(unique(f.data[:aquatic].C_W_1))

    for (i,C_W) in enumerate(unique(f.data[:aquatic].C_W_1))
        c += 1
        df = @subset(f.data[:aquatic], :C_W_1 .== 0)
        @df df lineplot!(plt_aqua, :t_exp, :wetmass_mg, (.5, .5), subplot = c, color = :gray, lw = 2, linestyle = :dash, label = "")
        @df df lineplot!(plt_aqua, :t_exp, :fract_tadpoles, (.5, .5), subplot = c+num_concs, color = :gray, lw = 2, linestyle = :dash)
    end

    plt = plot(
        plt_aqua, 
        kwargs...
        )

    return plt
end


function plot_metamorphs(;kwargs...)

    @df f.data[:metamorphs] plot(
        boxplot(string.(:C_W_1), :t_exp_G42),
        boxplot(string.(:C_W_1), :t_exp_G42),
        boxplot(string.(:C_W_1), :wetmass_G42_mg),
        boxplot(string.(:C_W_1), :wetmass_G46_mg)
    )

    @df f.data[:metamorphs] plot(
        boxplot(string.(:C_W_1), :t_exp_G42, ylim = (0,20)),
        boxplot(string.(:C_W_1), :t_exp_G42, ylim = (0,20)),
        boxplot(string.(:C_W_1), :wetmass_G42_mg, ylim = (50,400)),
        boxplot(string.(:C_W_1), :wetmass_G46_mg, ylim = (50,400)), 
        layout = (1,4), size = (1200,350),
        leg = [true false false false], label = "Observed", 
        xlabel = "2,4D (mg/L)", leftmargin = 7.5mm, bottommargin = 5mm, 
        fillcolor = :gray, markercolor = :black, fillalpha = .5,
        ylabel = ["Time since \n start of experiment (d)" "Time since \n start of experiment (d)" "Wet mass (mg)" "Wet mass (mg)"], 
        title = ["Timing of Gosner 42" "Timing of Gosner 46" "Mass at Gosner 42" "Mass at Gonser 46"], titlefontsize = 10, labelfontsize = 10, 
    )

end