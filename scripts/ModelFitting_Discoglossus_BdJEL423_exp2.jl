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
includet(scriptsdir("ModelFitting_Discoglossus_24D_UCLM.jl"))
includet(scriptsdir("ModelFitting_Discoglossus_24D_Ugent_exp2.jl"))
includet(srcdir("traits.jl"))
includet(srcdir("utils.jl"))

const PMOA_24D = "M" # only relevant to simulate combined effects of 2,4-D and Bd
const TREATMENT_MAP_JEL423 = OrderedDict(
    1 => "uninfected", 
    2 => "JEL423"
)

function aggregate_juveniles(juv::AbstractDataFrame)

    return combine(groupby(juv, [:t_since_mm, :treatment_id])) do df
        DataFrame(
            wetmass_mg_mean = mean(skipmissing(df.wetmass_mg)),
            wetmass_mg_sd = std(skipmissing(df.wetmass_mg)),
            wetmass_mg_var = var(skipmissing(df.wetmass_mg)),

            y_wetmass_mg_mean = mean(skipmissing(df.y_wetmass_mg)),
            y_wetmass_mg_var = var(skipmissing(df.y_wetmass_mg)),
            y_wetmass_mg_sd = std(skipmissing(df.y_wetmass_mg)),
        )

    end
end

function load_data_Bd(
    paths::OrderedDict,
    bd_treatment_name,
    treatment_map
    )

    data = OrderedDict()

    for (key,info) in pairs(paths)
        path, header = info
        data[key] = CSV.read(path, DataFrame, header = header) |> 
        x -> @subset(x, :pretreatment_24D .== 0.) 
        data[key][!,:treatment_id] = [reverse(treatment_map)[x] for x in data[key].treatment_bd]
    end

    # convert mass measurements to mg
    data[:juveniles_raw][!,:wetmass_mg] = data[:juveniles_raw].weight_g * 1e3
    # drop original column
    select!(data[:juveniles_raw], Not(:weight_g))
    "y_weight_g" in names(data[:juveniles_raw]) ? select!(data[:juveniles_raw], Not(:y_weight_g)) : nothing
    # normalize weights
    
    data[:juveniles_raw] = EcotoxSystems.relative_response(
        data[:juveniles_raw],
        [:wetmass_mg],
        :treatment_id;
        groupby_vars = [:t_since_mm]
    )

    # bdload column is not needed in juvenile data, 
    # has to be removed because missing values can cause issues
    select!(data[:juveniles_raw], Not(:bdload))

    data[:juveniles_aggregated] = aggregate_juveniles(data[:juveniles_raw])

    # for the Bd load data, we don't need the control
    data[:bdloads] = @subset(data[:bdloads], :treatment_bd .== bd_treatment_name)

    # return re-ordered data
    return OrderedDict(
        :juveniles_raw => data[:juveniles_raw],
        :juveniles_aggregated => data[:juveniles_aggregated], 
        :bdloads => data[:bdloads]
    )

end

function load_data_exp2_Bd()

    return load_data_Bd(
        OrderedDict(
            :juveniles_raw => [datadir("exp_raw", "UGent", "exp2", "juveniles.csv"), 1],
            :bdloads => [datadir("exp_raw", "UGent", "exp2", "bdloads_JEL423.csv"), 1],
        ),
        "JEL423",
        TREATMENT_MAP_JEL423
    )

end

function plot_data_Bd(
    data_loadfun,
    treatment_map
    )

    data = data_loadfun()

    cpal = palette([:cyan, :purple], 3)

    plt1 = @df data[:juveniles_raw] groupedviolin(
        string.(:t_since_mm), :wetmass_mg, group = :treatment_id, 
        side = :left,
        legend = :topleft, legendtitle = "Bd treatment", legendtitlefontsize = 8,
        label = hcat([treatment_map[id] for id in unique(:treatment_id)]...),
        palette = cpal, fillalpha = .5,
        xlabel = "Time since G42 (d)", ylabel = "Wet mass (mg)", 
        ylim = (100,500)
    )
    
    plt2 = @df data[:juveniles_raw] groupedviolin(
        string.(:t_since_mm), :y_wetmass_mg, group = :treatment_id, 
        side = :left,
        legend = false, 
        palette = cpal, fillalpha = .5,
        xlabel = "Time since G42 (d)", ylabel = "Control-normalized \n wet mass (-)", 
        ylim = (0,2)
    )

    plt3 = plot(xlabel = "Mean Bd load (#)", ylabel = "Probabilty")
    vline!(
        data[:bdloads].load_mean, 
        color = cpal[2], label = "Observed", lw = 1.5, linestyle = :dash, 
        xlim = (
            0,
            maximum(data[:bdloads].load_mean)*3
            ), 
        ylim = (0,1.1)
        )
    plt4 = plot(xlabel = "Bd load variance", leg = false)
    vline!(
        data[:bdloads].load_var, 
        color = cpal[2], lw = 1.5, linestyle = :dash,
        xlim = (
            0,
            maximum(data[:bdloads].load_var)*10
            ),
        ylim = (0,1.1)
        )

    plt = plot(
        plt1, plt2, plt3, plt4, 
        layout = (2,2), size = (800,500), 
        bottommargin = 5mm, leftmargin = 5mm,
        legend_background_color = :transparent,
        foreground_color_legend = nothing
        )

    return plt
end

plot_data_exp2_Bd() = plot_data_Bd(load_data_exp2_Bd, TREATMENT_MAP_JEL423)

"""
    plot_sims_Bd!(plt, sims::AbstractVector; label = "Simulation")

Plot simulations with Bd exposure.
"""
function plot_sims_Bd!(plt, sims::AbstractVector; label = "Simulation")

    existing_xticks = parse.(Float64, xticks(plt.subplots[1])[2])
    
    juveniles_raw = vcat(map(x->x[:juveniles_raw], sims)...) |> clean |> 
    x -> @subset(x, [t in existing_xticks for t in x.t_since_mm])
    
    @df juveniles_raw groupedviolin!(    
        string.(:t_since_mm),
        :wetmass_mg,
        group = :treatment_id,
        side = :right,
        color = :steelblue, 
        linecolor = :steelblue, 
        label = "", 
        subplot = 1,
        fillalpha = .2
    )

    @df juveniles_raw groupedviolin!(
        plt,
        string.(:t_since_mm),
        :y_wetmass_mg,
        group = :treatment_id,
        side = :right,
        color = :steelblue, 
        linecolor = :steelblue, 
        label = label, 
        subplot = 2,
        fillalpha = .2
    )

    bdloads = vcat(map(x->x[:bdloads], sims)...) |> clean

    @df bdloads histogram!(
        :load_mean, 
        normalize = :probability, 
        color = :steelblue, 
        lw = 0, 
        label = label, 
        subplot = 3,
        fillalpha = .2
    )

    @df bdloads histogram!(
        :load_var, 
        normalize = :probability, 
        color = :steelblue, 
        lw = 0, 
        subplot = 4,
        fillalpha = .2
    )

    return plt
    
end

plot_sims_exp2_Bd! = plot_sims_Bd!

function scenario_definition_Bd!(p::ComponentVector)::Nothing

    p.glb.t_max = 100. # setting simulation time conservatively, for cases where metamorphosis is delayed a lot
    p.glb.pathogen_inoculation_time_since_G42 = 16. # approximate time since abirth at which pathogen is added
    p.glb.pathogen_inoculation_dose = 1e3 # number of zoospores (thousands) added at specified time point
    p.glb.medium_renewals = [31. + (8/24)] # time of medium_renewal (pathogen removal), i.e. pathogen renewal; exposure lasted 8h
    p.glb.dX_in = [1e10, 1e10] # ad libitum feeding conditions

    return nothing

end

scenario_definition_exp2_Bd! = scenario_definition_Bd!

function define_defaultparams_exp2_Bd(
    posterior_summary_larvalfit::String = datadir("sims", SAVETAG_LARVALFIT, "posterior_summary.csv"),
    posterior_summary_juvenilefit::String = datadir("sims", SAVETAG_JUVENILEFIT, "posterior_summary.csv"),
    posterior_summary_correction_factors::String = datadir("sims", SAVETAG_CORRECTION_FACTORS, "posterior_summary.csv")
    )::ComponentVector

    p = construct_paramvector()
    scenario_definition_exp2_Bd!(p)
    set_species_params_exp2!(p, posterior_summary_larvalfit, posterior_summary_juvenilefit, posterior_summary_correction_factors)

    return p
end

function extract_bdload_data(sim::AbstractDataFrame)::DataFrame

    t2 = infer_timepoint(sim; t_post_first_metam = 26)
    sim_t = @subset(sim, :treatment_id .== 2, isapprox.(t2, :t))

    if nrow(sim_t)==0
        return DataFrame(
            load_mean = NaN, 
            load_var = NaN, 
            load_skew = NaN
        )
    else
        return DataFrame(
            load_mean = mean(sim_t.P_S .* 1e3) , # P_S needs to be back-converted from thousands to #
            load_var = var(sim_t.P_S .* 1e3),
            load_skew = skewness(collect(skipmissing(sim_t.P_S .* 1e3))),
        )
    end

end

function postprocess_simulation(sim, p; return_raw = false)

    # calculate total dry mass and wet mass
    sim[!,:drymass_mg] = sim.S .+ sim.E_mt 
    sim[!,:wetmass_mg] = [calc_wetmass(r, p.spc.watercontent_larvae, p.spc.watercontent_juveniles) for r in eachrow(sim)]

    # convert simulation time to time since start of experiment
    sim[!,:t_exp] = sim.t .- p.spc.emb_dev_time
    # convert simulation time to time since first metamorphosis            
    time_of_first_metam = infer_timepoint(sim; t_post_first_metam = 0)
    sim[!,:t_since_mm] = sim.t .- time_of_first_metam
    sim[!,:t_since_hatch] = sim.t .- age_at_birth(sim)


    # optionally, return the raw simulation output
    if return_raw
        return sim 
    end

    # convert raw simulation output to dataset

    sim_data = OrderedDict(
        :aquatic => extract_aquatic_data(sim),
        :metamorphs => extract_metamorph_data(sim, 0.),
        :juveniles_raw => extract_juvenile_data(sim),
        :bdloads => extract_bdload_data(sim)
    )

    sim_data[:juveniles_aggregated] = aggregate_juveniles(sim_data[:juveniles_raw])

    return sim_data

end

"""
Simulation that mimicks UGent "experiment 2".
"""
function simulator_exp2_Bd(
    p::EcotoxSystems.ComponentVector; # parameters and forcings
    return_raw::Bool = false, # return raw simulation output? if false, converts output to format of the data
    param_links::NamedTuple = (ind = link_ind_params!,), 
    kwargs... # additional arguments for ODE_simulator
    )

    preprocess_parameters!(p)
    let inoculation_dose = p.glb.pathogen_inoculation_dose
    # try
        p.glb.pathogen_inoculation_dose = 0.
        sim =  @replicates AmphiDEB.ODE_simulator(
            p, 
            param_links = param_links,
            maxiters = 1e3, # default is 1e5
            statevars_init = initialize_statevars_noexposure,
            callbacks = callbacks_exp2(),
            kwargs...
            ) 10
        
        sim[!,:treatment_id] .= 1
        sim[!,:treatment_bd] .= "uninfected"

        p.glb.pathogen_inoculation_dose = inoculation_dose
        sim_infected =  @replicates AmphiDEB.ODE_simulator(
            p, 
            param_links = param_links,
            maxiters = 1e5, # default is 1e5
            statevars_init = initialize_statevars_noexposure,
            callbacks = callbacks_exp2(),
            kwargs...
            ) 10

        sim_infected[!,:treatment_id] .= 2
        sim_infected[!,:treatment_bd] .= "JEL423"

        append!(sim, sim_infected)

        sim = postprocess_simulation(sim, p; return_raw=return_raw)

        return sim
        #catch e 
        #    # Log a task-specific message
        #    with_logger(logger) do
        #        @info("
        #        Encountered error in simulator: 
        #        $e 
        #        Error ocurred for the following parameter sample: 
        #        $p
        #        ")
        #    end
        #
        #    # write buffered message to logging file
        #    flush(io)
        #end
    end
end

function setup_modelfit_Bd(
    pmoa::AbstractString;
    data_loadfun::Function,
    data_plotfun::Function,
    define_defaultparams::Function,
    simulator::Function,
    loss_functions = EcotoxModelFitting.loss_euclidean_logtransform,
    alt_priors = Dict(
        :sigma1 => Truncated(Normal(1.6, 1.6), 0, Inf), # sporangia killing rate [d^-1 1e3 S^-1]
    )
    )

    if !isdir(datadir("sims", SAVETAG))
        mkdir(datadir("sims", SAVETAG))
    end

    pmoa_idx = findfirst(x -> x == pmoa, PMOAS)

    # set up logger
    global io = open(datadir("sims", SAVETAG, "log.txt"), "w+")
    global logger = SimpleLogger(io)

    data = data_loadfun()

    f = ModelFit( 
        prior = Prior(
            "pth.gamma" => Truncated(Normal(0.1, 0.1), 0, Inf), # zoospore encounter rate [d^-1],
            "pth.eta" => Truncated(Normal(100, 50), 0, 200), # zoospore production rate  [Z S^-1 d^-1],
            "pth.sigma1" => alt_priors[:sigma1],
            "spc.Chi" => Hyperdist( # killing rate modifier [-]
                σ -> LogNormal(σ^2, σ),
                Truncated(Normal(1, 1), 0.5, 2)
            ),
            "spc.E_P[$(pmoa_idx)]" => Truncated(Normal(15, 7.5), 1.5, 45), # median effective load (10^3 sporangia)
            "spc.B_P[$(pmoa_idx)]" => Truncated(Normal(2, 2), 2, 4) # bd effect slope [-]
            ),
        defaultparams = define_defaultparams(), 
        simulator = simulator, 
        data = data, 
        response_vars = [
            Symbol[],
            [:y_wetmass_mg_mean], 
            [:load_mean, :load_var],
        ],
        time_resolved = [false, true, false],
        data_weights = [
            Float64[],
            [1.0],
            [1.0, 1.0]
        ],
        time_var = :t_since_mm,
        grouping_vars = [
            Symbol[],
            [:t_since_mm, :treatment_id],
            Symbol[]
            ], 
        plot_data = data_plotfun,
        loss_functions = loss_functions
    )

    return f

end


function setup_modelfit_exp2_Bd(pmoa::AbstractString)

    f = setup_modelfit_Bd(
        pmoa;
        data_loadfun = load_data_exp2_Bd,
        data_plotfun = plot_data_exp2_Bd,
        define_defaultparams = define_defaultparams_exp2_Bd, 
        simulator = simulator_exp2_Bd
    )

    return f

end

paramlabels["pth.gamma"] = L"\gamma_{Bd}"
paramlabels["pth.eta"] = L"\eta_{Bd}"
paramlabels["pth.sigma1"] = L"\sigma_{1,Bd}"
paramlabels["spc.Chi"] = L"\chi_{Bd}"
paramlabels["spc.E_P[2]"] = L"e_{Bd,M}"
paramlabels["spc.B_P[2]"] = L"b_{Bd,M}"

function fit_model!(
    f; 
    pmcsettings = (
        :n => 100,
        :q_dist => 0.1, 
        :t_max => 0,
        :evals_per_sample => 1
    ),
    savetag = nothing, 
    paramlabels = paramlabels, 
    continue_from = nothing,
    n_posterior_check = 100,
    kwargs...
    )

    @info "Using savetag $(savetag)"

    let pmchist, posterior_check

            
        # run the calibration
        @suppress pmchist = run_PMC!(
            f; 
            pmcsettings...,
            savedir = datadir("sims"),
            savetag = savetag, 
            paramlabels = paramlabels, 
            continue_from = continue_from
        )

        closeall() # reset plots

        @info "#### ---- Best fit ---- ####"

        let plt = f.plot_data()

            p_opt = f.accepted[:,argmin(vec(f.losses))]
            sim_opt = [f.simulator(p_opt) for _ in 1:100]

            plot_sims_exp2_Bd!(plt, sim_opt, label = "Best fit")
            display(plt)
            
            if !isnothing(savetag)
                savefig(plot(plt, dpi = 400), datadir("sims", "$(savetag)", "VPC_posterior_bestfit.png"))
            end
        end
        
        @info "#### ---- Posterior retrodictions ---- ####"

        let plt = f.plot_data()
            
            @suppress posterior_check = posterior_predictions(f, n_posterior_check)

            plot_sims_exp2_Bd!(plt, posterior_check.predictions, label = "Retrodictions")

            if !isnothing(savetag)
                savefig(plot(plt, dpi = 400), datadir("sims", "$(savetag)", "VPC_posterior_samples.png"))
            end

            display(plt)
        end

        @info "#### ---- Marginal posteriors ---- ####"

        let plt, num_params = length(f.prior.dists), num_cols = 4, num_rows = Int(ceil(num_params / num_cols))

            plt = plot(
                plot.(f.prior.dists, color = :black)..., layout = (num_rows, num_cols), 
                leg = hcat(vcat(true, repeat([false], num_params-1))...), 
                label = "Prior", 
                size = (1200,200*num_rows), bottommargin = 5mm
                )
        
                
            for (i,param) in enumerate(f.prior.labels)
        
                histogram!(
                    plt, subplot = i, 
                    f.accepted[i,:], weights = Weights(vec(f.weights)), 
                    xlabel = param in keys(paramlabels) ? paramlabels[param] : param, 
                    normalize = :pdf, label = "Posterior", color = :gray, lw = 0.5, fillalpha = .5
                    )
            end
        
            display(plt)
            
            if !isnothing(savetag)
                savefig(plot(plt, dpi = 400), datadir("sims", "$(savetag)", "marginal_posteriors.png"))
            end
        end

        @info "#### ---- Convergence behaviour ---- ####"

        let plt
            plt = plot(eachindex(pmchist.dists) .- 1, map(median, pmchist.dists), marker = true, lw = 1.2, xlabel = "PMC step", ylabel = "Loss", label = "Median", xticks = eachindex(pmchist.dists) .- 1)
            plot!(plt, eachindex(pmchist.dists) .- 1, map(minimum, pmchist.dists), marker = true, lw = 1.2, label = "Minimum")
            display(plt)
            savefig(plot(plt, dpi = 400), datadir("sims", savetag, "loss.png"))
        end

        begin
            posterior_medians = [mapslices(x -> median(x, Weights(pmchist.weights[i])), pmchist.particles[i], dims = 2) for i in eachindex(pmchist.particles)] |>
            x -> hcat(x...)

            q25 = [mapslices(x -> quantile(0.25, x, Weights(pmchist.weights[i])), pmchist.particles[i], dims = 2) for i in eachindex(pmchist.particles)] |>
            x -> hcat(x...)
            
            q75 = [mapslices(x -> quantile(0.25, x, Weights(pmchist.weights[i])), pmchist.particles[i], dims = 2) for i in eachindex(pmchist.particles)] |>
            x -> hcat(x...)
            
            num_params = length(f.prior.dists) 
            num_cols = 4 
            num_rows = Int(ceil(num_params / num_cols))

            plt = plot(
                eachindex(pmchist.particles) .-1,
                posterior_medians', layout = size(posterior_estimates)[1], 
                size = (1200,200*num_rows), marker = true,
                leg = hcat(vcat(:topleft, repeat([false], length(f.prior.dists)-1))...),
                label = "Median",
                ylabel = hcat(f.prior.labels...), titlefontsize = 12,
                bottommargin = 5mm, leftmargin = 5mm, 
                xlabel = "PMC step"
                )

            plot!(
                eachindex(pmchist.particles) .-1, 
                q25', fillrange = q75', 
                lw = 0, fillcolor = :gray, fillalpha = .25,
                marker = :diamond, label = "IQR"
                )


            for (i,dist) in enumerate(f.prior.dists)

                # set ylim based in prior limits

                q1 = quantile(dist, 0.01)
                q2 = quantile(dist, 0.99)

                # indicate IQR of priors 

                l = repeat([quantile(dist, 0.25)], length(pmchist.particles))
                u = repeat([quantile(dist, 0.75)], length(pmchist.particles))

                plot!(plt, subplot = i, ylim = (q1,q2))

                #plot!(
                #    plt, subplot = i, 
                #    eachindex(pmchist.particles) .- 1, 
                #    l,
                #    fillrange = u, 
                #    fillalpha = 0.35, 
                #    color = :lightgray, 
                #    label = "Prior IQR"
                #    )

                hline!(plt, subplot = i, [quantile(dist, 0.25)], color = :gray, linestyle = :dash, label = "")
                hline!(plt, subplot = i, [quantile(dist, 0.75)], color = :gray, linestyle = :dash, label = "")

            end
        end

        display(plt)

        @info "#### ---- Posterior summary ---- ####"

        generate_posterior_summary(
            f.prior.labels,
            f.accepted, 
            f.losses, 
            f.weights; 
            tex = false,
            paramlabels = paramlabels,
            savetag = nothing
        ) |> display

        return pmchist, posterior_check
    end
end