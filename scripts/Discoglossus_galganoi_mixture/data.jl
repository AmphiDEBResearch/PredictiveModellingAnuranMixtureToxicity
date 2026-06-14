# functions used to fit 2,4D-models
# the functions defined here perform the model fitting for larval and juvenile life stages

# packages

using Base.Threads

using StatsBase
using DataFrames, DataFramesMeta
using DataStructures, CSV
using Chain
using StatsPlots, Plots.Measures
using Logging

default(leg=false)
theme(:default)

using LaTeXStrings
using Suppressor
using Distributions

using EcotoxSystems, AmphiDEB, EcotoxModelFitting
import AmphiDEB: ComponentVector

# source files

using Revise
includet(scriptsdir("traits.jl"))
includet(scriptsdir("utils.jl"))
include(scriptsdir("Discoglossus_galganoi_24D", "fit.jl"))

# fully factorial exposure matrix - including the combinations which have not been tested (we will simulate them)
const EXPOSURE_MATRIX =   [
                0. 0.;
                0. 10.;
                0. 100.;
                0.03 0.;
                0.03 10.;
                0.03 100.;
                0.3  0.;
                0.3 10.;
                0.3 100.;
                3.0 0.;
                3.0 10.;
                3.0 100.;
                30.0 0.;
                30.0 10.;
                30.0 100.;
                100. 0.;
                100. 10.;
                100. 100.
            ]

const TREATMENT_IDS_MIX = DataFrame(
    EXPOSURE_MATRIX, 
    [:D_ppm, :F_ppm]
    ) |> x->@transform(x, :treatment_id = 1:nrow(x))

function load_data_UCLM_mix(;
    paths::OrderedDict=OrderedDict(
        :aquatic => [datadir("exp_raw", "UCLM", "discoglossus_galganoi", "disco_flp_24D_aquatic.csv"), 1], # number indicates row where data header is located (omitting metadata)
        :metamorphs => [datadir("exp_raw", "UCLM", "discoglossus_galganoi", "disco_flp_24D_metamorphs.csv"), 1],
        ##:adults => [datadir("exp_raw", "Discoglossus_03_adults.csv"), 1]
    ))

    data = OrderedDict()

    for (key, info) in pairs(paths)
        path, header = info
        data[key] = CSV.read(path, DataFrame, header=header)
    end

    # process aquatic data

    dropmissing!(data[:aquatic])

    data[:aquatic].num_tadpoles = float.(data[:aquatic].num_tadpoles)
    data[:aquatic][!,:fract_tadpoles] = data[:aquatic].num_tadpoles ./ data[:aquatic].survival

    leftjoin!(
        data[:aquatic], TREATMENT_IDS_MIX, on = [:D_ppm, :F_ppm]
    )

    rename!(data[:aquatic])
    data[:aquatic] = EcotoxSystems.relative_response(
        data[:aquatic],
        [:wetmass_mg, :num_tadpoles],
        :treatment_id;
        groupby_vars=[:t_exp, :aquarium]
    )

    # process metamorph data

    data[:metamorphs].wetmass_G42_mg = float.(data[:metamorphs].wetmass_G42_mg)
    data[:metamorphs].wetmass_G46_mg = float.(data[:metamorphs].wetmass_G46_mg)

    leftjoin!(
        data[:metamorphs], TREATMENT_IDS_MIX, on=[:D_ppm, :F_ppm]
    )
    select!(data[:metamorphs], [:treatment_id, :D_ppm, :F_ppm, :wetmass_G42_mg, :wetmass_G46_mg, :t_exp_G42, :t_exp_G46])
    rename!(data[:metamorphs])
    dropmissing!(data[:metamorphs])

    return data
end

function plot_data_UCLM_mix_growth(; kwargs...)

    data = load_data_UCLM_mix()

    plt_aquatic = plot(
        layout=(
            length(unique(data[:aquatic].D_ppm)),
            length(unique(data[:aquatic].F_ppm))
        ),
        size=(800, 1000),#, thickness_scaling = 1.15,  
        leg=false,
        xlim = (0,50),
        ylim = (100,400),
    )

    c = 0
    for (i, D) in enumerate(unique(data[:aquatic].D_ppm))
        for (j, F) in enumerate(unique(data[:aquatic].F_ppm))
            c += 1

            df = @subset(data[:aquatic], :F_ppm .== F, :D_ppm .== D)

            @df df scatter!(
                :t_exp, :wetmass_mg;
                subplot=c, color=:black,
                leg=c == 1 ? :topright : false, label="Observed",
                title="$D | $F",
                markersize = 3,
                kwargs...
            )

        end
    end
    plt_aquatic

    return plt_aquatic
end

function plot_data_UCLM_mix_numtadpoles(; kwargs...)

    data = load_data_UCLM_mix()
   
    plt_aquatic = plot(
        layout=(
            length(unique(data[:aquatic].D_ppm)),
            length(unique(data[:aquatic].F_ppm))
        ),
        size = (800, 1000),#, thickness_scaling = 1.15,  
        leg=false,
        xlim = (0,50),
        ylim = (-0.5, 10.5),
    )

    c = 0
    for (i, D) in enumerate(unique(data[:aquatic].D_ppm))
        for (j, F) in enumerate(unique(data[:aquatic].F_ppm))
            c += 1

            df = @subset(data[:aquatic], :F_ppm .== F, :D_ppm .== D)

            @df df scatter!(
                :t_exp, :num_tadpoles;
                subplot=c, color=:black,
                leg=c == 1 ? :topright : false, label="Observed",
                title="$D | $F",
                markersize = 3,
                kwargs...
            )

            

        end
    end
    plt_aquatic

    return plt_aquatic
end

function plot_data_UCLM_mix_fracttadpoles(; kwargs...)

    data = load_data_UCLM_mix()
   
    plt_aquatic = plot(
        layout=(
            length(unique(data[:aquatic].D_ppm)),
            length(unique(data[:aquatic].F_ppm))
        ),
        size = (800, 1000),#, thickness_scaling = 1.15,  
        leg=false,
        xlim = (0,50),
        ylim = (-0.05, 1.05),
    )

    c = 0
    for (i, D) in enumerate(unique(data[:aquatic].D_ppm))
        for (j, F) in enumerate(unique(data[:aquatic].F_ppm))
            c += 1

            df = @subset(data[:aquatic], :F_ppm .== F, :D_ppm .== D)

            @df df scatter!(
                :t_exp, :fract_tadpoles;
                subplot=c, color=:black,
                leg=c == 1 ? :topright : false, label="Observed",
                title="$D | $F",
                markersize = 3,
                kwargs...
            )

            

        end
    end
    plt_aquatic

    return plt_aquatic
end

function plot_data_UCLM_mix_metamorphs()
    
    data = load_data_UCLM_mix()
    data[:metamorphs] = @subset(data[:metamorphs], :D_ppm .<= 100)
    num_flp_treatments = length(unique(data[:metamorphs].F_ppm))

    plt_metamorphs = plot(
        layout = (2,num_flp_treatments), 
        size = (1000,500), 
    )

    for (i,F) in enumerate(unique(data[:metamorphs].F_ppm))

        @df data[:metamorphs] dotplot!(
            string.(unique(:D_ppm)), #string.(unique(EXPOSURE_MATRIX[:,1:end-1])), 
            repeat([0], length(unique(:D_ppm))), 
            markersize = 0, markeralpha = 0,
            label = "", 
            subplot = i
        )

        @df data[:metamorphs] dotplot!(
            string.(unique(:D_ppm)), #string.(unique(EXPOSURE_MATRIX[:,1:end-1])), 
            repeat([0], length(unique(:D_ppm))), 
            markersize = 0, markeralpha = 0,
            label = "", 
            subplot = i+num_flp_treatments
        )


        df = @subset(data[:metamorphs], :F_ppm .== F)
        sort!(df, :D_ppm)

        @df df violin!(
            plt_metamorphs, subplot = i,
            string.(:D_ppm), :t_exp_G46, 
            color = :gray, side = :left,
            #xticks = unique(string.(data[:metamorphs].D_ppm)), 
            fillalpha = .15, 
            ylabel = i == 1 ? "Time to G46 (d)" : "",
            leftmargin = i == 1 ? 5mm : 2.5mm,  
            label = i == 1 ? "Observed" : "", 
            title = "$F mg/L FLP", 
            ylim = (10, 50)
            )

        @df df dotplot!(
            plt_metamorphs, subplot = i,
            string.(:D_ppm), :t_exp_G46, 
            color = :black, side = :left, label = ""
            )

        @df df violin!(
            plt_metamorphs, subplot = i+num_flp_treatments,
            string.(:D_ppm), :wetmass_G46_mg, 
            color = :gray, side = :left,
            #xticks = unique(string.(data[:metamorphs].D_ppm)), 
            fillalpha = .15, 
            ylabel = i == 1 ? "Wet mass \n at G46 (d)" : "",
            leftmargin = i == 1 ? 5mm : 2.5mm,  
            label = "", 
            title = "", 
            ylim = (50,300)
            )

        @df df dotplot!(
            plt_metamorphs, subplot = i+num_flp_treatments,
            string.(:D_ppm), :wetmass_G46_mg, 
            color = :black, side = :left, label = ""
            )

    end

    return plt_metamorphs

end


function leftcol(numcols = 3)
    return map(x->((x+numcols-1)%numcols)==0, unique(sims[1][:aquatic].treatment_id)) |> 
    x -> findall(x->x==true, x)
end

function toprow(numcols = 3)
    return 1:numcols
end