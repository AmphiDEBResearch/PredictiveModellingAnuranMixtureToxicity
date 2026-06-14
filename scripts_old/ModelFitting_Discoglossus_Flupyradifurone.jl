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

function plot_sims!(plt, predictions::AbstractVector; label = "Simulation", color = :steelblue)::Nothing

    df_aqua = sort(vcat([df[:aquatic] for df in predictions]...), :t_exp)

    c = 0
    num_concs = length(unique(df_aqua.C_W_1))

    for (i,C_W) in enumerate(unique(df_aqua.C_W_1))
        c += 1
        df = @subset(df_aqua, :C_W_1 .== C_W)
        @df df lineplot!(plt, :t_exp, :wetmass_mg, subplot = c, color = color, lw = 2, fillalpha = .2, label = label)
        @df df lineplot!(plt, :t_exp, :fract_tadpoles, subplot = c+num_concs, color = color, lw = 2, fillalpha = .2)
    end

    return nothing
end


function define_defaultparams()::AmphiDEB.ComponentVector

    p = AmphiDEB.ComponentVector(
        glb = AmphiDEB.defaultparams.glb, 
        pth = AmphiDEB.defaultparams.pth,
        spc = EcotoxSystems.ComponentVector(
            AmphiDEB.defaultparams.spc; 
            # auxiliary parameters
            log_k_D_G = log(1e-10),
            log_k_D_M = log(1e-10),
            log_k_D_A = log(1e-10),
            log_k_D_KAP = log(1e-10),
            log_e_G = log(1e10),
            log_e_M = log(1e10),
            log_e_A = log(1e10),
            log_e_KAP = log(1e10),
            watercontent_larvae = 0.93, 
            watercontent_juveniles = 0.85,
            time_since_birth = 15.,
            emb_dev_time = 2. 
        ))

    # setting global parameters

    p.glb.t_max = 100. # setting simulation time conservatively, for cases where metamorphosis is delayed a lot
    p.glb.pathogen_inoculation_time = Inf # no pathogen inoculation
    p.glb.dX_in = [1e10, 1e10] # ad libitum feeding conditions

    p.spc.Z = truncated(Normal(1, 0.1), 0, Inf)
    # propagation of zoom factor to H_j1 is turned off => we want variability in the transition to metamorphs
    p.spc.propagate_zoom.H_j1 = 0.

    # adding point estimates as defaults

    posterior_summary_larvalfit = CSV.read(datadir("sims", SAVETAG_LARVALFIT, "posterior_summary.csv"), DataFrame)
    posterior_summary_juvenilefit = CSV.read(datadir("sims", SAVETAG_JUVENILEFIT, "posterior_summary.csv"), DataFrame)

    @info "Overwriting default values of $(posterior_summary_larvalfit.param)"
    for (label,value) in zip(posterior_summary_larvalfit.param, posterior_summary_larvalfit.best_fit)
        if label == "spc.Z"
            p.spc.Z = truncated(Normal(1, value), 0, Inf)
        else
            assign_value_by_label!(p, label, value)
        end
    end

    @info "Overwriting default values of $(posterior_summary_juvenilefit.param)"
    for (label,value) in zip(posterior_summary_juvenilefit.param, posterior_summary_juvenilefit.best_fit)
        if label == "spc.Z"
            p.spc.Z = truncated(Normal(1, value), 0, Inf)
        else
            assign_value_by_label!(p, label, value)
        end
    end

    p.spc.k_M_juv = p.spc.k_M_emb
    p.spc.X_emb_int = 1. # ≈ initial dry mass of an egg (mg)

    p.spc.KD .= 0.
    p.spc.B .= 2.
    p.spc.E .= 1e10

    p.spc.emb_dev_time = estimate_emb_dev_time(p)

    return p
end

"""
    simulate_preexperiment(p::AmphiDEB.ComponentVector)::AmphiDEB.ComponentVector

Simulate individual life history before the start of the experiment. <br>
This is used to simulate different conditions before and during the experiment. 
"""
function simulate_preexperiment(p::AmphiDEB.ComponentVector)::AmphiDEB.ComponentVector
    
    # remember actual values of t_max and C_W
    let t_max = deepcopy(p.glb.t_max), 
        C_W = deepcopy(p.glb.C_W)

        # set simulation time to estimated time at start of the experiment
        p.glb.t_max = p.ind.emb_dev_time + p.ind.time_since_birth

        # turn off chemical exposure
        p.glb.C_W = [0.]

        # simulate the model until the start of the experiment
        sim = AmphiDEB.ODE_simulator(
            p, 
            returntype = EcotoxSystems.odesol, # directly return the ODE solution object - we don't need a DataFrame
            gen_ind_params = x -> x # skip generation of individual-level parameters - we already have them
            )
        
        # re-set global parameters to actual values
        p.glb.t_max = t_max
        p.glb.C_W .= C_W

        # retrieve the final state pf the simulated experiment
        u0 = sim.u[end].ind

        # return the final state as initial state of the actual simulation, 
        # re-setting global states
        return ComponentVector(
            glb = AmphiDEB.initialize_global_statevars(p),
            ind = u0
        )
    end

end

function link_ind_params!(ind::AmphiDEB.ComponentVector)::Nothing
    
    # defining links between parameters 
    # this is mostly relevant for parameters which may be subject to the zoom factor
    # this is handled by EcotoxSystems.jl for each simulated individual 

    ind.dI_max_emb = ind.dI_max_lrv # ingestion rate for embryos assumed to be same as for larave
    
    # FIXME: this does not appear to have the desired effect
    # the issue is currently fixed by adding expression to simulator; works for now because k_M is not linked to Z
    # this should be fixed though and a test added in EcotoxSystems.jl and AmphiDEB.jl
    ind.k_M_juv = ind.k_M_emb # somatic maintenace rate is assumed to remain constant across life stages
    
    ind.k_J_emb = (1-ind.kappa_emb)/ind.kappa_emb * ind.k_M_emb # maturity maintenance is linked to somatic (assuming same cumulative investment in both branches)
    #ind.k_J_juv = (1-ind.kappa_juv)/ind.kappa_juv * ind.k_M_juv
    ind.kappa_juv = ind.kappa_emb

    return nothing
end

"""
    simulator(
        p::EcotoxSystems.ComponentVector; # parameters and forcings
        return_raw::Bool = false, # return raw simulation output? if false, converts output to format of the data
        param_links::NamedTuple = (ind = link_ind_params!,),
        kwargs... # additional arguments for ODE_simulator
        )


Simulate toxicity test with *D. galganoi* exposed to 2,4D.
"""
function simulator(
    p::EcotoxSystems.ComponentVector; # parameters and forcings
    return_raw::Bool = false, # return raw simulation output? if false, converts output to format of the data
    param_links::NamedTuple = (ind = link_ind_params!,), 
    kwargs... # additional arguments for ODE_simulator
    )

    # parameter pre-processing

    let Z = deepcopy(p.spc.Z)

        p.glb.C_W .= 0. # reset C_W => exposure() will take care of this
        p.spc.time_since_birth = ceil(p.spc.time_since_birth) # convert time since birth to whole day
        # add parameter links

        p.spc.dI_max_emb = p.spc.dI_max_lrv
        p.spc.k_M_juv = p.spc.k_M_emb

        # estimate population mean of the embryonic development time

        p.spc.Z = Dirac(1.) # turn off individual variabiilty to get estimate based on popmean 
        
        p.spc.emb_dev_time = estimate_emb_dev_time(p) # assign embryonic development time
        p.spc.Z = Z # re-assign original zoom factor
        
        # reverting log-transformed parameters

        #p.spc.KD[1,1] = exp(p.spc.log_k_D_G)
        #p.spc.KD[1,2] = exp(p.spc.log_k_D_M)
        #p.spc.KD[1,3] = exp(p.spc.log_k_D_A)
        #p.spc.KD[1,6] = exp(p.spc.log_k_D_KAP)
        #
        #p.spc.E[1,1] = exp(p.spc.log_e_G)
        #p.spc.E[1,2] = exp(p.spc.log_e_M)
        #p.spc.E[1,3] = exp(p.spc.log_e_A)
        #p.spc.E[1,6] = exp(p.spc.log_e_KAP)

        # "inner simulator" is the function called for each treatment
        try
            inner_sim(p) = @replicates AmphiDEB.ODE_simulator(
                p, 
                param_links = param_links,
                maxiters = 1e5, # default is 1e5
                statevars_init = simulate_preexperiment,
                #model = AmphiDEB.AmphiDEB_ODE_with_linear_TD!,
                kwargs...
                ) 10

            # epxosure() runs inner_sim for each treatment and collects the results
    
            sim = exposure(
                inner_sim, 
                p,
                Matrix(hcat([0., 10., 100.]...)')
            )

            # calculate total dry mass and wet mass
            sim[!,:drymass_mg] = sim.S .+ sim.E_mt 
            sim[!,:wetmass_mg] = [calc_wetmass(r, p.spc.watercontent_larvae, p.spc.watercontent_juveniles) for r in eachrow(sim)]
            
            # convert simulation time to time since start of experiment
            sim[!,:t_exp] = sim.t #.- age_at_birth(@subset(sim, :C_W_1 .== 0)) .- p.spc.time_since_birth

            # optionally, return the raw simulation output
            if return_raw
                return sim 
            end

            # convert raw simulation output to dataset

            sim_data = OrderedDict(
                :aquatic => extract_aquatic_data(sim),
                :metamorphs => extract_metamorph_data(sim, 0.),
                #:adults => extract_adult_data(sim)
            )

            return sim_data
        catch e 
            # Log a task-specific message
            with_logger(logger) do
                @info("
                Encountered error in simulator: 
                $e 
                Error ocurred for the following parameter sample: 
                $p
                ")
            end

            # write buffered message to logging file
            flush(io)
        end
    end
end

function extract_aquatic_data(sim::AbstractDataFrame)::DataFrame
 
    aquatic = @chain sim begin
        @select(:t_exp, :larva, :wetmass_mg, :replicate, :C_W_1)
        @subset(:larva .> 0.5) # only larvae 
        groupby([:t_exp, :C_W_1]) # for every time point and treatment
        combine(
            [:wetmass_mg,:replicate] => 
            ((m, r) -> 
                (wetmass_mg = mean(m), # calculate average dry mass
                num_tadpoles = length(unique(r))) # count how many tadpoles we have
                ) => AsTable
            )
    end

    aquatic[!,:fract_tadpoles] = aquatic.num_tadpoles ./ maximum(aquatic.num_tadpoles)

    return aquatic
end

function extract_metamorph_data(sim::AbstractDataFrame, time_since_hatch)::DataFrame

    # calculate metamorphosis traits for each replicate and treatment
    metamorphs = combine(groupby(sim, [:replicate, :C_W_1])) do df
        DataFrame(
            t_exp_G42 = metamorphosis_timing(df)[1] .- age_at_birth(df) .- time_since_hatch, # timing of Gosner 42 (days since start of experiment)
            t_exp_G46 = metamorphosis_timing(df)[2] .- age_at_birth(df) .- time_since_hatch, # timing of Gosner 46 (days since start of experiment)
            dt_mt = metamorphosis_duration(df), # duration of metamorphosis (d)
            wetmass_G42_mg = wetmass_at_G42(df), # wet mass at Gosner 42 (mg)
            wetmass_G46_mg = wetmass_at_G46(df) # wet mass at Gosner 46 (mg)
        )
    # calculate averages per treatment
    end |> x -> combine(groupby(x, :C_W_1)) do df
        DataFrame(
            t_exp_G42 = mean(df.t_exp_G42),
            t_exp_G46 = mean(df.t_exp_G46),
            dt_mt = mean(df.dt_mt),
            wetmass_G42_mg = mean(df.wetmass_G42_mg),
            wetmass_G46_mg = mean(df.wetmass_G46_mg)
        )
    end
    
    return metamorphs
end

function setup_modelfit(
    pmoa::AbstractString; 
    loss_function::Function = EcotoxModelFitting.loss_euclidean
    )

    pmoa_idx = findfirst(x -> x == pmoa, PMOAS) # convert pmoa from string to index
    @assert !isnothing(pmoa_idx) "Did not find PMoA $(pmoa) in PMOAS"

    if !isdir(datadir("sims", SAVETAG))
        mkdir(datadir("sims", SAVETAG))
    end
    # set up logger
    # FIXME: this does not save the log file to the intended directory, because the setup_modelfit does not know the composite "savetag" 
    global io = open(datadir("sims", SAVETAG, "log.txt"), "w+")
    global logger = SimpleLogger(io)

    data = load_data()

    f = ModelFit( 
        prior = Prior(
            "spc.KD[1,$(pmoa_idx)]" => Truncated(Normal(0.5, 0.5), 0.001, 1), 
            "spc.E[1,$(pmoa_idx)]" => Truncated(Normal(50, 50), 1, Inf),
            "spc.B[1,$(pmoa_idx)]" => Truncated(Normal(2., 4), 0.5, Inf)
            ),
        defaultparams = define_defaultparams(), 
        simulator = simulator, 
        data = data, 
        response_vars = [
            [:wetmass_mg, :fract_tadpoles], 
            Symbol[] # leaving this empty == metamorphs will be ignored during calibration
            #[:t_exp_G42, :t_exp_G46, :wetmass_G42_mg, :wetmass_G46_mg]
        ],
        data_weights = [
            [1., 1.],
            Float64[]
        ],
        grouping_vars = [[:C_W_1], [:C_W_1]], 
        time_resolved = [true, false],
        time_var = :t_exp,
        plot_data = plot_data,
        loss_functions = loss_function
    )

    return f

end


function fit_model!(
    f; 
    n = 100, 
    n_init = 100, 
    q_dist = .1, 
    t_max = 3, 
    savetag = nothing, 
    paramlabels = paramlabels, 
    continue_from = nothing,
    n_posterior_check = 100,
    kwargs...
    )

    if !isdir(plotsdir(savetag))
        mkdir(plotsdir(savetag))
    end

    pmchist = run_PMC!(
        f; 
        n = n, 
        n_init = n_init, 
        q_dist = q_dist, 
        t_max = t_max, 
        savetag = savetag, 
        paramlabels = paramlabels, 
        evals_per_sample = 3,
        continue_from = continue_from
    )

    @info "#### ---- Best fit ---- ####"

    let plt = f.plot_data()

        p_opt = f.samples[:,argmin(vec(f.losses))]
        sim_opt = [f.simulator(p_opt) for _ in 1:100]

        save_sims(sim_opt, savetag, "VPC_bestfit")
        plot_sims!(plt, sim_opt, label = "Best fit")
        display(plt)

        if !isnothing(savetag)
            savefig(plot(plt, dpi = 300), datadir("sims", savetag, "VPC_bestfit.png"))
        end
 
    end

    @info "#### ---- Maximum a posteriori estimation ---- ####"

    let plt = f.plot_data()

        p_MAP = f.samples[:,argmax(vec(f.weights))]
        sim_MAP = [f.simulator(p_MAP) for _ in 1:100]

        plot_sims!(plt, sim_MAP, label = "MAP")
        display(plt)
        
        if !isnothing(savetag)
            savefig(plot(plt, dpi = 400), datadir("sims", "$(savetag)", "VPC_MAP.png"))
        end

    end

    @info "#### ---- Posterior samples ---- ####"

    let plt = f.plot_data()
        @suppress_err begin
            global posterior_check = posterior_predictions(f, n_posterior_check)
        end

        save_sims(posterior_check.predictions, savetag, "VPC_posterior")
        plot_sims!(plt, posterior_check.predictions, label = "Retrodictions")

        if !isnothing(savetag)
            savefig(plot(plt, dpi = 400), datadir("sims", savetag, "VPC_posterior.png"))
        end

        display(plt)
    end

    @info "#### ---- Posterior summary ---- ####"

    display(pmchist.posterior_summary)

    
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
                f.samples[i,:], weights = Weights(vec(f.weights)), 
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
        posterior_estimates = [mapslices(x -> x[argmax(Weights(pmchist.weights[i]))], pmchist.particles[i], dims = 2) for i in eachindex(pmchist.particles)] |>
        x -> hcat(x...)

        posterior_medians = [mapslices(x -> median(x, Weights(pmchist.weights[i])), pmchist.particles[i], dims = 2) for i in eachindex(pmchist.particles)] |>
        x -> hcat(x...)

        num_params = length(f.prior.dists) 
        num_cols = 4 
        num_rows = Int(ceil(num_params / num_cols))

        plt = plot(
            eachindex(pmchist.particles) .-1,
            posterior_estimates', layout = size(posterior_estimates)[1], 
            size = (1200,200*num_rows), marker = true,
            leg = hcat(vcat(:topleft, repeat([false], length(f.prior.dists)-1))...),
            label = "MAP",
            ylabel = hcat(f.prior.labels...), titlefontsize = 12,
            bottommargin = 5mm, leftmargin = 5mm, 
            xlabel = "PMC step"
            )

        plot!(eachindex(pmchist.particles) .-1, posterior_medians', marker = :diamond, label = "Median")

        for (i,dist) in enumerate(f.prior.dists)

            # set ylim based in prior limits

            q1 = quantile(dist, 0.01)
            q2 = quantile(dist, 0.99)

            # indicate IQR of priors 

            l = repeat([quantile(dist, 0.25)], length(pmchist.particles))
            u = repeat([quantile(dist, 0.75)], length(pmchist.particles))

            plot!(plt, subplot = i, ylim = (q1,q2))

            plot!(
                plt, subplot = i, 
                eachindex(pmchist.particles) .- 1, 
                l,
                fillrange = u, 
                fillalpha = 0.35, 
                color = :lightgray, 
                label = "Prior IQR"
                )

            hline!(plt, subplot = i, [quantile(dist, 0.25)], color = :gray, linestyle = :dash, label = "")
            hline!(plt, subplot = i, [quantile(dist, 0.75)], color = :gray, linestyle = :dash, label = "")

        end
    end

    display(plt)

    @info "#### ---- Posterior summary ---- ####"

    generate_posterior_summary(
        f.prior.labels,
        f.samples, 
        f.losses, 
        f.weights; 
        tex = false,
        paramlabels = paramlabels,
        savetag = nothing
    ) |> display



    return pmchist, posterior_check
end

function plot_pmc_gradient(pmchist)


    plt = plot(layout = (1,2))

    loss_median = [median(x) for x in pmchist.dists]
    loss_minimum = [minimum(x) for x in pmchist.dists]
    loss_median_gradient = []
    loss_minimum_gradient = []
    
    for i in eachindex(loss_median)
        if i==1
            push!(loss_median_gradient, NaN)
            push!(loss_minimum_gradient, NaN)
        else
            push!(loss_median_gradient, -(loss_median[i]-loss_median[i-1])/loss_median[1])
            push!(loss_minimum_gradient, -(loss_minimum[i]-loss_minimum[i-1])/loss_minimum[1])
        end
    end
    
    x = eachindex(loss_median) .- 1 
    plot!(x, loss_median, marker = true, label = "Median", xlabel = "SMC population", ylabel = "Loss", subplot = 1)
    plot!(x, loss_minimum, marker = true, label = "Minimum", subplot = 1)
    
    plot!(x, loss_median_gradient, marker = true, subplot = 2, label = "Median", xlabel = "SMC population", ylabel = "Normalized gradient", leg = false)
    plot!(x, loss_minimum_gradient, marker = true, subplot = 2, label = "Minimum")

    return plt
end


function plot_posteriors(f::ModelFit, pmchist)
    plt = plot(layout = (2, 4), size = (1000,500), bottommargin = 7.5mm)

    for t in eachindex(pmchist.particles)
        for k in eachindex(f.prior.dists)

            if t == 1 
                if typeof(f.prior.dists[k]) == Hyperdist
                    plot!(f.prior.dists[k].dist, subplot = k, color = :black, lw = 2, linestyle = :dash, label = "Prior")
                else
                    plot!(f.prior.dists[k], subplot = k, color = :black, lw = 2, linestyle = :dash, label = "Prior")
                end
            end

            density!(
                plt, subplot = k,
                pmchist.particles[t][k,:], weights = Weights(pmchist.weights[t]), 
                leg = k == 1 ? :topright : false, 
                fill = true, fillalpha = .15, lw = 1, 
                label = "t = $t",
                xlabel = paramlabels[f.prior.labels[k]], xlabelfontsize = 14, 
                yticks = [], yaxis = false
                )
        end
    end

    return plt
end


# misc

paramlabels = OrderedDict(
    "spc.eta_AS_emb" => L"\eta_{AS}^{emb}",
    "spc.eta_AS_juv" => L"\eta_{AS}^{juv}",
    "spc.eta_AR" => L"\eta_{AR}",
    "spc.dI_max_lrv" => L"\{ \dot{I} \}_{max}^{lrv}",
    "spc.dI_max_juv" => L"\{ \dot{I} \}_{max}^{juv}",
    "spc.k_M_emb" => L"k_M^{emb}",
    "spc.Z" => L"\sigma_{Z,M}",
    "spc.H_j1_prime" => L"H^{j'}",
    "spc.H_p_prime" => L"H^{p'}",
    "spc.gamma" => L"\gamma",
    "spc.kappa_emb" => L"\kappa",
    "spc.watercontent_larvae" => L"Larval\ water\ content",
    "spc.watercontent_juveniles" => L"Juvenile\ water\ content"
)


function save_sims(
    predictions::AbstractVector, 
    savetag::AbstractString, 
    prefix::AbstractString
    )::Nothing

    for key in keys(predictions[1])
        df = vcat([@transform(p[key], :num_sample = i) for (i,p) in enumerate(predictions)]...)
        CSV.write(datadir("sims", savetag, "$(prefix)_$(key).csv"), df)
    end

    return nothing
end
