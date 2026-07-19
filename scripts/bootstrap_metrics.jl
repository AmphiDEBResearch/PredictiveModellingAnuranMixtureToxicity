# ============================================================= #
# Functions to bootstrap quantitative metrics on PMoA selection
# Implementation supported by GPT 5.4-nano and revised manually
# ============================================================= #

using Random
import Base.rand
using StatsBase

"""
Compute quantitative metrics for aquatic data (larvae up to GS 42).

"""
function get_metrics(f::ModelFit, sims::Vector)

    obs = combine(groupby(f.data[:aquatic], [:C_W_1, :t_exp])) do df 
        DataFrame(
            :wetmass_mg = mean(df.wetmass_mg), 
            :fract_tadpoles = mean(df.fract_tadpoles)
        )
    end

    sim = [@transform(x[:aquatic], :num = i) for (i,x) in enumerate(sims)] |> 
    x -> vcat(x...) |> 
    x -> combine(groupby(sim, [:C_W_1, :t_exp])) do df
        df = @subset(df, isfinite.(:wetmass_mg), isfinite.(:fract_tadpoles))
        DataFrame(
            :wetmass_mg = mean(df.wetmass_mg), 
            :fract_tadpoles = mean(df.fract_tadpoles)
        )
    end

    joined = leftjoin(
        obs, sim, 
        on = [:t_exp, :C_W_1], makeunique = true, 
        renamecols = "_obs"=>"_sim"
        )

    joined = joined[
        completecases(
            joined, 
            [:wetmass_mg_sim, :wetmass_mg_obs, :fract_tadpoles_sim, :fract_tadpoles_obs]),
        :]

    joined[ismissing.(joined.fract_tadpoles_sim),:fract_tadpoles_sim] .= 0.

    # calculate metrics for each simulation
    metrics = combine(groupby(joined, [:C_W_1])) do df
        return DataFrame(
            MAPE_wetmass = MAPE(df.wetmass_mg_obs, df.wetmass_mg_sim), 
            MAPE_fracttadpoles =  MAPE(df.fract_tadpoles_obs, df.fract_tadpoles_sim), 
            NSE_wetmass = NSE(df.wetmass_mg_obs, df.wetmass_mg_sim), 
            NSE_fracttadpoles = NSE(df.fract_tadpoles_obs, df.fract_tadpoles_sim), 
            )
    end #|> # aggergate over simulations
    #x -> combine(groupby(x, :C_W_1)) do df
    #    subdf_NSE = @subset(df, isfinite.(:NSE_fracttadpoles))
    #    return DataFrame(
    #        MAPE_wetmass = mean(df.MAPE_wetmass), 
    #        MAPE_fracttadpoles = mean(df.MAPE_fracttadpoles), 
    #        NSE_wetmass = mean(df.NSE_wetmass), 
    #        NSE_fracttadpoles = mean(subdf_NSE.NSE_fracttadpoles)
    #    )
    #end

    return metrics
end

"""
Bootstrap confidence interval for a quantitative metric by resampling `sims` (with replacement),
using the same logic as `get_metrics`.

Returns one row per C_W_1 with:
- MAPE_point: MAPE from the original `sims`
- MAPE_lower, MAPE_upper: percentile CI with coverage (1-α)
"""
function get_MAPE_CI(
    f::ModelFit,
    sims::Vector;
    num_iters::Int = 2000,
    α::Float64 = 0.05,
    metric::Symbol = :MAPE_wetmass
    )

    n = length(sims)

    # Point estimate (no resampling)
    point = get_metrics(f, sims) |> x -> select(x, :C_W_1, metric)
    rename!(point, metric => :MAPE)

    # Bootstrap replicates
    boot_dfs = Vector{DataFrame}(undef, num_iters)
    for b in 1:num_iters
        idx = sample(1:n, n)  # resample sims with replacement
        m_b = get_metrics(f, sims[idx])    # compute metrics on resample
        df_b = select(m_b, :C_W_1, metric)
        rename!(df_b, metric => :mape)
        boot_dfs[b] = df_b
    end

    all = vcat(boot_dfs...)


    # Percentile CI by C_W_1
    ci = combine(groupby(all, :C_W_1)) do df
        df = df[completecases(df, [:mape]),:]
        DataFrame(
            C_W_1 = first(df.C_W_1),
            MAPE_lower = quantile(df.mape, α/2),
            MAPE_upper = quantile(df.mape, 1 - α/2),
        )
    end

    return leftjoin(ci, point, on = :C_W_1)
end

using Random, Statistics, DataFrames
using StatsBase

function get_metric_CI(
    f::ModelFit,
    sims::Vector;
    num_iters::Int = 2000,
    α::Float64 = 0.05,
    metric::Symbol = :MAPE_wetmass,
    point_col::Symbol = :MAPE
)

    lower_col = Symbol(string(point_col), "_lower")
    upper_col = Symbol(string(point_col), "_upper")

    n = length(sims)

    # Point estimate
    point = get_metrics(f, sims) |> x -> select(x, :C_W_1, metric)
    rename!(point, metric => point_col)

    # Bootstrap replicates
    val_col = :val
    boot_dfs = Vector{DataFrame}(undef, num_iters)
    for b in 1:num_iters
        idx = sample(1:n, n; replace=true)
        m_b = get_metrics(f, sims[idx])

        df_b = select(m_b, :C_W_1, metric)
        rename!(df_b, metric => val_col)
        boot_dfs[b] = df_b
    end

    all = vcat(boot_dfs...)

    ci = combine(groupby(all, :C_W_1)) do df
        df = df[completecases(df, [val_col]), :]
        c = first(df.C_W_1)
        l = quantile(df[!, val_col], α/2)
        u = quantile(df[!, val_col], 1 - α/2)

        DataFrame(
            :C_W_1 => [c],
            lower_col => [l],
            upper_col => [u],
        )
    end

    return leftjoin(ci, point, on = :C_W_1)
end
