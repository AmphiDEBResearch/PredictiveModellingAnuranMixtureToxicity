# ======================================== #
# Cross-validation with metamorph data
# ======================================== #

using DataFrames, DataFramesMeta, CSV 
using StatsPlots, Plots.Measures
using EcotoxSystems, AmphiDEB
using StatsBase
default(leg = false)
theme(:default)

include(scriptsdir("utils.jl"));

using Distances
using OrdinaryDiffEq

function aggregate_metamorphs(data)

    metamorphs_aggregated = combine(groupby(data[:metamorphs], :treatment_id)) do df
        return DataFrame(
            t_exp_G46 = mean(df.t_exp_G46),
            wetmass_G46_mg = mean(df.wetmass_G46_mg)
        )
    end

    return metamorphs_aggregated
end

function plot_metamorphs(
    ;
    xlabel = "2,4-D (mg/L)",
    leftmargin = 7.5mm, 
    bottommargin = 7.5mm, 
    size = (1500,400)
    )

    plt = @df f.data[:metamorphs] plot(    
        violin(
            string.(:treatment_id), :t_exp_G46, 
            side = :left, 
            xticks = (unique(:treatment_id) .- 0.5, unique(:C_W_1)), 
            xrotation = 45
            ), 
        violin(
            string.(:treatment_id), :wetmass_G46_mg,
            side = :left,
            xticks = (unique(:treatment_id) .- 0.5, unique(:C_W_1)),
            xrotation = 45
        ), 
        layout = (1,4), 
        leg = [true false false false], label = "Observed", 
        fillcolor = :gray, markercolor = :black, fillalpha = .5,
        ylabel = ["Time since \n start of experiment (d)" "Wet mass (mg)"], 
        title = ["Timing of Gosner 46" "Mass at Gonser 46"], 
        titlefontsize = 10, labelfontsize = 10, 
        xlabel = xlabel,
        leftmargin = leftmargin, 
        bottommargin = bottommargin, 
        size = size
    )

    @df f.data[:metamorphs] dotplot!(
        string.(:treatment_id), :t_exp_G46, 
        side = :left, 
        color = :black,
        label = "", 
        subplot = 1
    )

    @df f.data[:metamorphs] dotplot!(
        string.(:treatment_id), :wetmass_G46_mg, 
        side = :left, 
        color = :black,
        label = "", 
        subplot = 2
    )

    return plt

end


mean_relative_error(a,b) = mean(@. (a-b)/b)
mape(a,b) = 100*mean(abs.(a .- b) ./ b)

function quant_eval_metamorphs(
    f::ModelFit, 
    sims::AbstractVector
    )

    sims_df = map(x->x[:metamorphs], sims) |>
    x -> [@transform(df, :num_sim = i) for (i,df) in enumerate(x)] |> 
    x -> vcat(x...) |> 
    x -> @transform(x, :treatment_id = [TREATMENT_IDS[c] for c in x.C_W_1]) |>
    clean
            
    eval_df = leftjoin(
        aggregate_metamorphs(f.data),
        sims_df, 
        on = :treatment_id, 
        makeunique = true, 
        renamecols = :_obs => :_pred
    )

    combine(groupby(eval_df, :num_sim_pred)) do df

        df = @subset(df, :treatment_id .> 1)

        mrevals = [
            mean_relative_error(df.t_exp_G46_pred, df.t_exp_G46_obs),
            mean_relative_error(df.wetmass_G46_mg_pred, df.wetmass_G46_mg_obs)
        ]

        mapevals = [
            mape(df.t_exp_G46_pred, df.t_exp_G46_obs),
            mape(df.wetmass_G46_mg_pred, df.wetmass_G46_mg_obs)
        ]

        nrmsdvals = [
            nrmsd(df.t_exp_G46_pred, df.t_exp_G46_obs),
            nrmsd(df.wetmass_G46_mg_pred, df.wetmass_G46_mg_obs)
        ]

        DataFrame(
            endpoint = ["t_exp_G46", "wetmass_G46_mg"], 
            mre = mrevals, 
            mape = mapevals, 
            nrmsd = nrmsdvals
        )
    end |> 
    x -> combine(groupby(x, :endpoint)) do df
        DataFrame(
            mre = mean(df.mre),
            mre_sd = std(df.mre),
            nrmsd = mean(df.nrmsd),
            nrmsd_sd = std(df.nrmsd),
            mape = mean(df.mape),
            mape_sd = std(df.mape)
        )
    end

end