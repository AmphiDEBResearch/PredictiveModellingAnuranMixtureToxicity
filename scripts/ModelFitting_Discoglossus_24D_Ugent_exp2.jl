using OrdinaryDiffEq
import AmphiDEB: ComponentVector
using EcotoxModelFitting
import EcotoxModelFitting: Hyperdist

includet(scriptsdir("ModelFitting_Discoglossus_24D_UCLM.jl"))
includet(srcdir("utils.jl"))
includet(srcdir("loss.jl"))

function load_data_exp2_noBd(;
    paths::OrderedDict = OrderedDict(
        :juveniles => [datadir("exp_raw", "UGent", "exp2", "juveniles.csv"), 1]
    ))

    data = OrderedDict()

    for (key,info) in pairs(paths)
        path, header = info
        data[key] = CSV.read(path, DataFrame, header = header) |> 
        x -> @subset(x, :treatment_bd .== "uninfected") 

        data[key].t_since_mm = float.(data[key].t_since_mm)
    end

    treatment_ids = Dict(zip(
        sort(unique(data[:juveniles].pretreatment_24D)),
        collect(eachindex(sort(unique(data[:juveniles].pretreatment_24D))))
    ))

    # convert mass measurements to mg
    data[:juveniles][!,:wetmass_mg] = data[:juveniles].weight_g * 1e3
    # drop original column
    select!(data[:juveniles], Not(:weight_g))
    rename!(data[:juveniles], :y_weight_g => :y_wetmass_mg)

    data[:juveniles][!,:treatment_id] = [treatment_ids[x] for x in data[:juveniles].pretreatment_24D]
    data[:juveniles_control] = @subset(data[:juveniles], :treatment_id .== 1)[:,[:t_since_mm, :wetmass_mg]] 


    data[:juveniles_agg] = combine(groupby(data[:juveniles], [:t_since_mm, :pretreatment_24D, :treatment_bd, :treatment_id])) do df
        DataFrame(
            y_wetmass_mg = mean(df.y_wetmass_mg), 
            wetmass_mg = mean(df.wetmass_mg)
        )

    end

    # re-ordering keys
    data = OrderedDict(
        :juveniles_control => data[:juveniles_control], 
        :juveniles => data[:juveniles], 
        :juveniles_agg => data[:juveniles_agg]
        )

    return data

end

function plot_data_exp2_noBd()

    exp2 = load_data_exp2_noBd()

    plt1 = @df @subset(exp2[:juveniles], :treatment_bd .== "uninfected") groupedviolin(
        string.(:t_since_mm), :wetmass_mg, group = :pretreatment_24D, 
        side = :left,
        legend = :topleft, legendtitle = "2,4D-pretreatment \n (mg/L)", legendtitlefontsize = 8,
        palette = palette([:cyan, :purple], 3), fillalpha = .5,
        xlabel = "Time since G42 (d)", ylabel = "Wet mass (mg)", 
        ylim = (100,500)
    )
    @df @subset(exp2[:juveniles], :treatment_bd .== "uninfected") groupeddotplot!(
        plt1,
        string.(:t_since_mm), :wetmass_mg, group = :pretreatment_24D, 
        side = :left, color = :black, label = ""
        )

    
    plt2 = @df @subset(exp2[:juveniles], :treatment_bd .== "uninfected") groupedviolin(
        string.(:t_since_mm), :y_wetmass_mg, group = :pretreatment_24D, 
        side = :left,
        legend = false, 
        palette = palette([:cyan, :purple], 3), fillalpha = .5,
        xlabel = "Time since G42 (d)", ylabel = "Control-normalized \n wet mass (-)", 
        ylim = (0,2)
    )

     @df @subset(exp2[:juveniles], :treatment_bd .== "uninfected") groupeddotplot!(
        plt2,
        string.(:t_since_mm), :y_wetmass_mg, group = :pretreatment_24D, 
        side = :left,
        color = :black
        )

    plot(
        plt1, plt2, 
        layout = (1,2), size = (800,350), 
        bottommargin = 5mm, leftmargin = 5mm,
        legend_background_color = :transparent,
        foreground_color_legend = nothing
        )

end

"""
    plot_sims_exp2_noBd!(plt, sims::AbstractVector; label = "Simulation")

For UGent experiment 2, plot simulations on top of data. 
"""
function plot_sims_exp2_noBd!(plt, sims::AbstractVector; label = "Simulation")
    
    juveniles = vcat(map(x->x[:juveniles], sims)...) |> clean
    
    @df juveniles groupedviolin!(    
        string.(:t_since_mm), :wetmass_mg,
        group = :pretreatment_24D,
        side = :right,
        color = :steelblue, 
        linecolor = :steelblue, 
        label = "", 
        subplot = 1, 
        fillalpha = .2
    )

    @df juveniles groupedviolin!(
        plt,
        string.(:t_since_mm), :y_wetmass_mg,
        group = :pretreatment_24D,
        side = :right,
        color = :steelblue, 
        linecolor = :steelblue, 
        label = label, 
        subplot = 2, 
        fillalpha = .2
    )

    return plt
    
end

function construct_paramvector()
    p = ComponentVector(
    glb = ComponentVector(
        AmphiDEB.defaultparams.glb; 
        chemical_addition_time = 26., # time since hatching at which chemical is added 
        chemical_exposure_duration = 5., # duration of chemical exposure 
        pathogen_inoculation_time_since_G42 = 16.,  # inoculation time as time since metamorphosis (Gosner 42)
    ), 
    pth = AmphiDEB.defaultparams.pth,
    spc = ComponentVector(
        AmphiDEB.defaultparams.spc; 
        # auxiliary parameters
        watercontent_larvae = 0.93, 
        watercontent_juveniles = 0.85,
        time_since_birth = 15., # approximate time since birth at the beginning of the experiment
        time_to_G42 = 30., # expected time to metamorphosis (estimated will be refined)
        emb_dev_time = 2.,
        Z_UCLM = truncated(Normal(1, 0.1), 0, Inf), # Z estimated from UCLM data - will be assigned from file later
        Z_mean_UGent = 1., # size correction factor for UGent vs UCLM data
        H_j1_UCLM = 10, # H_j1 estimate from UCLM 
    ))

    return p

end

function scenario_definition_exp2_noBd!(p)

    # setting global parameters

    p.glb.t_max = 100. # [d] setting simulation time conservatively, for cases where metamorphosis is delayed a lot
    p.glb.pathogen_inoculation_time_since_G42 = 16. # [days since birth]
    p.glb.pathogen_inoculation_dose = 0  # [1e3 spores]
    p.glb.dX_in = [1e10, 1e10] # [mg/d] simulating ad libitum feeding conditions

end

function set_species_params_exp2!(
    p, 
    posterior_summary_larvalfit, 
    posterior_summary_juvenilefit,
    posterior_summary_correction_factors,
    )::Nothing
    
    p.spc.Z = truncated(Normal(1, 0.1), 0, Inf)
    # if propagation of zoom factor to H_j1 is included, the variability in time to metamorphosis becomes very small
    p.spc.propagate_zoom.H_j1 = 1.

    # adding point estimates as defaults

    assign_values_from_file!(
        p, 
        posterior_summary_larvalfit,
        exceptions = OrderedDict("spc.Z" => (p,label,value) -> p.spc.Z_UCLM = Truncated(Normal(1, value), 0, 1))
        
    )

    assign_values_from_file!(
        p, 
        posterior_summary_juvenilefit,
        exceptions = OrderedDict()
        )

    assign_values_from_file!(
        p, 
        posterior_summary_correction_factors,
        exceptions = OrderedDict()
        )
    
    p.spc.Z = truncated(Normal(p.spc.Z_mean_UGent, p.spc.Z_mean_UGent*p.spc.Z.untruncated.σ), 0, Inf)
    p.spc.H_j1 = p.spc.H_j1_UCLM * mode(p.spc.Z)
    p.spc.KD .= 0.
    p.spc.B .= 2.
    p.spc.E .= 1e10
    
    p.pth.gamma = 0. # disable pathogen in the default params; only once we start fitting will the pathogen be simulated

    p.spc.k_M_juv = p.spc.k_M_emb
    p.spc.X_emb_int = 1. # ≈ initial dry mass of an egg (mg)
    p.spc.emb_dev_time = estimate_emb_dev_time(p)
    p.spc.time_to_G42 = metamorphosis_timing(AmphiDEB.ODE_simulator(p))[1]
    p.glb.pathogen_inoculation_time = p.glb.pathogen_inoculation_time_since_G42 + p.spc.time_to_G42 # update the actual inoculation time (since initialization)

    return nothing

end

function define_defaultparams_exp2_noBd(
    posterior_summary_larvalfit::String = datadir("sims", SAVETAG_LARVALFIT, "posterior_summary.csv"),
    posterior_summary_juvenilefit::String = datadir("sims", SAVETAG_JUVENILEFIT, "posterior_summary.csv"),
    posterior_summary_correction_factors::String = datadir("sims", SAVETAG_CORRECTION_FACTORS, "posterior_summary.csv"),
    )::ComponentVector

    p = construct_paramvector()
    scenario_definition_exp2_noBd!(p)
    set_species_params_exp2!(
        p, 
        posterior_summary_larvalfit, 
        posterior_summary_juvenilefit, 
        posterior_summary_correction_factors
        )

    return p
end

function link_ind_params!(ind::ComponentVector)::Nothing
    
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

## Callbacks
# these are used to reproduce events specific to the experimental setup

function exposure_condition(u, t, integrator) 
    # chemical is added at a fixed time point after hatching
    time_since_hatch = t - integrator.p.ind[:emb_dev_time]
    timing_condition = integrator.p.glb[:chemical_addition_time] - time_since_hatch

    return timing_condition
end

function exposure_effect!(integrator)
    integrator.u.glb.C_W .= integrator.p.glb[:C_W]
end

function nonexposure_condition(u, t, integrator)
    # chemical is removed at a fixed time point after it has been added

    time_since_hatch = t - integrator.p.ind[:emb_dev_time]
    timing_condition = integrator.p.glb[:chemical_addition_time] + integrator.p.glb[:chemical_exposure_duration] - time_since_hatch

    return timing_condition
end

function nonexposure_effect!(integrator)
    integrator.u.glb.C_W .= 0.
end

# pathogen inoculation occurrs at a fixed time after gosner 42
condition_inoculation(u, t, integrator) = t - integrator.p.glb.pathogen_inoculation_time
function effect_inoculation!(integrator) 
    integrator.u.glb.P_Z = integrator.p.glb.pathogen_inoculation_dose
end

# optional: at a given time, individuals are transferred to clean medium and the external spore abundance is re-set to 0
# afterwards, individuals continue shedding spores, so this does not always have a considerable effect on the dynamics
condition_renewal(u,t,integrator) = prod(integrator.p.glb.medium_renewals  .- (t - integrator.p.ind[:emb_dev_time]))  
function effect_renewal!(integrator)
    integrator.u.glb.P_Z = 0.
end

"""
    callbacks_exp2()

Defines callbacks to mimick UGent-exp2.
"""
function callbacks_exp2()

    cb_inoculation = ContinuousCallback(
        condition_inoculation, 
        effect_inoculation!
        )

    cb_renewal = ContinuousCallback(
        AmphiDEB.condition_renewal, 
        AmphiDEB.effect_renewal!
        )

    cb_exposure = ContinuousCallback(
        exposure_condition, 
        exposure_effect!
        )

    cb_nonexposure = ContinuousCallback(
        nonexposure_condition, 
        nonexposure_effect!
        )

    return CallbackSet(
        cb_inoculation, 
        cb_renewal, 
        cb_exposure,
        cb_nonexposure
        )
end

"""
    initialize_statevars_noexposure(p::ComponentVector)

Initialize state variables with C_W set to 0, assuming that callbacks will handle the rest. 
"""
function initialize_statevars_noexposure(p::ComponentVector)

    global_statevars = AmphiDEB.initialize_global_statevars(p)
    global_statevars.C_W .= 0.
    individual_statevars = AmphiDEB.initialize_individual_statevars(p)

    return ComponentVector(
        glb = global_statevars,
        ind = individual_statevars
    )

end

"""
    infer_t2(sim::AbstractDataFrame)::Float64

Infers the final time point of weight measurements in UGent experiment 2. 
This is 26 days afyter the first individual has reached metamorphosis. 
"""
function infer_timepoint(sim::AbstractDataFrame; t_post_first_metam = 26)::Float64

    t_first_metam = @subset(sim, :metamorph .> 0.9).t_exp |> minimum
    return t_first_metam + t_post_first_metam

end

"""
    extract_juvenile_data(sim::AbstractDataFrame)::DataFrame

Computes juveniles data measured in UGent experiment.
"""
function extract_juvenile_data(sim::AbstractDataFrame)::DataFrame

    return @chain sim begin
        @select(:t_since_mm, :wetmass_mg, :treatment_id)
        EcotoxSystems.relative_response(
            [:wetmass_mg],
            :treatment_id;
            groupby_vars = [:t_since_mm]
        )    
    end
    
end

"""
map treatment_id to pretreatment_24D
"""
const TREATMENT_MAP = Dict(
        1 => 0.,
        2 => 0.03, 
        3 => 0.3
    )


    
"""
Mimicks UGent "experiment 2" with optional 2,4D-preexposure during experiment 1. Does not consider Bd exposure.
"""

function preprocess_parameters!(p)

    p.glb.C_W .= 0. # reset C_W => exposure() will take care of this
    p.spc.time_since_birth = ceil(p.spc.time_since_birth) # convert time since birth to whole day
    
    # add parameter links

    p.spc.dI_max_emb = p.spc.dI_max_lrv
    p.spc.k_M_juv = p.spc.k_M_emb

    # estimate population mean of the embryonic development time

    p.spc.Z = Dirac(1.) # turn off individual variabiilty to get estimate based on popmean 
    p.spc.emb_dev_time = estimate_emb_dev_time(p) # assign embryonic development time
    p.spc.Z = Truncated(Normal(p.spc.Z_mean_UGent, p.spc.Z_mean_UGent * p.spc.Z_UCLM.untruncated.σ), 0, Inf)
    # zoom factor does not affect maturity threshold so that we model individual variability correctly, 
    # but we do want the size correction => applying zoom factor separately
    p.spc.H_j1 = p.spc.H_j1_UCLM * p.spc.Z_mean_UGent 
    
    return nothing
end

function postprocess_simulation(sim, p; return_raw = return_raw) 

    sim[!,:treatment_bd] .= "uninfected"

    # calculate total dry mass and wet mass
    sim[!,:drymass_mg] = sim.S .+ sim.E_mt 
    sim[!,:wetmass_mg] = [calc_wetmass(r, p.spc.watercontent_larvae, p.spc.watercontent_juveniles) for r in eachrow(sim)]
    
    ## convert simulation time to time since start of experiment, using first time of metamorphosis
    sim[!,:t_exp] = sim.t .- p.spc.emb_dev_time
    ## convert simulation time to time since first metamorphosis
    #time_of_first_metam = infer_timepoint(sim; t_post_first_metam = 0)
    sim[!,:t_since_mm] = sim.t .- p.spc.time_to_G42
    sim[!,:t_since_hatch] = sim.t .- age_at_birth(sim)

    # optionally, return the raw simulation output
    if return_raw
        return sim 
    end

    # convert raw simulation output to dataset

    sim_data = OrderedDict(
        :aquatic => extract_aquatic_data(sim),
        :metamorphs => extract_metamorph_data(sim, 0.),
        :juveniles => extract_juvenile_data(sim)
    )

    sim_data[:juveniles_agg] = combine(groupby(sim_data[:juveniles], [:t_since_mm, :treatment_id])) do df
        DataFrame(
            y_wetmass_mg = mean(skipmissing(df.y_wetmass_mg)), 
            wetmass_mg = mean(skipmissing(df.wetmass_mg))
        )
    end

    sim_data[:juveniles][!,:pretreatment_24D] = [TREATMENT_MAP[t] for t in sim_data[:juveniles].treatment_id]
    sim_data[:juveniles_agg][!,:pretreatment_24D] = [TREATMENT_MAP[t] for t in sim_data[:juveniles_agg].treatment_id]
    sim_data[:juveniles_control] = @subset(sim_data[:juveniles], :treatment_id .== 1)

    return sim_data

end

function define_simulator_exp2_noBd(;
    C_Wvec = [0, 0.03, 0.3]
    )


    # "inner simulator" is the function called for each treatment
    function inner_sim(p)

        sim = @replicates AmphiDEB.ODE_simulator(
                p, 
                param_links = (ind=link_ind_params!,),
                maxiters = 1e5, # default is 1e5
                statevars_init = initialize_statevars_noexposure,
                callbacks = callbacks_exp2(),
                ) 10 

        return sim

    end

    function simulator_exp2_noBd(
        p::EcotoxSystems.ComponentVector; # parameters and forcings
        return_raw::Bool = false, # return raw simulation output? if false, converts output to format of the data
        kwargs... # additional arguments for ODE_simulator
        )

        preprocess_parameters!(p)

        try

            # epxosure() runs inner_sim for each treatment and collects the results

            sim = exposure(
                inner_sim,
                p,
                C_Wvec
            )

            return postprocess_simulation(sim, p; return_raw = return_raw)
        catch e 
            error(e)
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

    return simulator_exp2_noBd

end

function fit_control_exp2(;
    loss_functions = EcotoxModelFitting.loss_euclidean_logtransform,
    pmcsettings = (
        :n => 1000,
        :q_dist => 0.1,
        :t_max => 3
    ))

    if !isdir(datadir("sims", SAVETAG))
        mkdir(datadir("sims", SAVETAG))
    end

    data = load_data_exp2_noBd()
    data = OrderedDict(:juveniles_control => data[:juveniles_control])

    function plot_control_exp2()

        plt = @df data[:juveniles_control] violin(
            string.(:t_since_mm), :wetmass_mg, 
            side = :left, color = :gray, fillalpha = .25,
            label = "Observed", xlabel = "Time since G42 (d)", ylabel = "Wet mass"
            )

        @df data[:juveniles_control] dotplot!(
            string.(:t_since_mm), :wetmass_mg, 
            side = :left, color = :black, label = ""
            )

        return plt
    end

    function plot_sims_control_exp2!(plt, sims::AbstractVector; label = "Simulation") 

        juveniles_control = EcotoxModelFitting.extract_simkey(sims, :juveniles_control)

        @df juveniles_control violin!(
            plt,
            string.(:t_since_mm), :wetmass_mg, 
            side = :right, color = :steelblue, fillalpha = .25,
            label = label
        )

        return plt

    end

    f = ModelFit( 
        prior = Prior(
            "spc.Z_mean_UGent" => Truncated(Normal(0.8, 0.8), 0.1, 1), 
            "spc.eta_AS_juv" => Truncated(Normal(0.9, 0.1), 0, 1),
            ),
        defaultparams = define_defaultparams_exp2_noBd(), 
        simulator = define_simulator_exp2_noBd(;C_Wvec = [0.]), 
        data = data, 
        response_vars = [
            [:wetmass_mg],
        ],
        data_weights = [
            [1.], 
        ],
        grouping_vars = [
            [:t_since_mm],
            ], 
        time_resolved = [false],
        time_var = :t,
        plot_data = plot_control_exp2,
        loss_functions = loss_functions
    )

    _ = fit_model!(
        f; 
        pmcsettings = pmcsettings,
        evals_per_sample = 1, 
        savetag = "$(SAVETAG)_control", 
        plot_sims! = plot_sims_control_exp2!
        )

    p_opt = EcotoxModelFitting.bestfit(f)

    return ComponentVector(Dict(zip(Symbol.(f.prior.labels), p_opt)))

end

function setup_modelfit_exp2_noBd(
        pmoa::AbstractString; 
        p_opt_control::ComponentVector,
        loss_functions = EcotoxModelFitting.loss_euclidean_logtransform
    )


    pmoa_idx = findfirst(x -> x == pmoa, PMOAS)

    # set up logger
    # FIXME: this does not save the log file to the intended directory, because the setup_modelfit does not know the composite "savetag" 
    global io = open(datadir("sims", SAVETAG, "log.txt"), "w+")
    global logger = SimpleLogger(io)

    data = load_data_exp2_noBd()
    defaultparams = define_defaultparams_exp2_noBd()
    EcotoxModelFitting.assign!(defaultparams, p_opt_control)

    f = ModelFit( 
        prior = Prior(
            "spc.KD[1,$(pmoa_idx)]" => Truncated(Normal(1., 1), 0, 1),
            "spc.E[1,$(pmoa_idx)]" => Truncated(Normal(3e-3, 0.3), 3e-5, 1),
            "spc.B[1,$(pmoa_idx)]" => Truncated(Normal(2, 8), 1, 10)
            ),
        defaultparams = defaultparams, 
        simulator = define_simulator_exp2_noBd(), 
        data = data, 
        response_vars = [
            [:wetmass_mg],
            [:y_wetmass_mg], 
            [:y_wetmass_mg]
        ],
        data_weights = [
            [0.], 
            [0.],
            [1.]
        ],
        grouping_vars = [
            [:t_since_mm],
            [:t_since_mm, :treatment_id],
            [:t_since_mm, :treatment_id]
            ], 
        time_resolved = [false, false, false],
        time_var = :t,
        plot_data = plot_data_exp2_noBd,
        loss_functions = loss_functions
    )

    return f

end

lowsettings = (
        :n_init => 500,
        :n => 500,
        :q_dist => 0.1,
        :t_max => 3,
        :evals_per_sample => 10,
    )

highsettings = (
        :n_init => 50_000,
        :n => 25_000, 
        :q_dist => 1000/25_000,  
        :t_max => 10, 
        :evals_per_sample => 10,
)

function fit_exposure_exp2(
    pmoa::AbstractString; 
    p_opt_control::ComponentVector,
    loss_functions = EcotoxModelFitting.loss_euclidean_logtransform,
    pmcsettings = lowsettings
    )

    if !isdir(datadir("sims", SAVETAG))
        mkdir(datadir("sims", SAVETAG))
    end
    
    f = setup_modelfit_exp2_noBd(
        pmoa; 
        p_opt_control = p_opt_control,
        loss_functions = loss_functions
    )

    fit_model!( # ETA 20 minutes
        f; 
        pmcsettings = pmcsettings,
        plot_sims! = plot_sims_exp2_noBd!,
        savetag = "$(SAVETAG)_$(pmoa)", 
        paramlabels = paramlabels
    );

    return f

end

function plot_raw_sim(f::ModelFit, p::Vector{Float64}; n = 10)

    sim = [f.simulator(p; return_raw = true) for _ in 1:n] |> x-> vcat(x...)

    plt = @df sim plot(
        groupedlineplot(:t, :S .+ :E_mt, :treatment_id, ylabel = "W", fillalpha = .2, lw = 2, label = [0 0.03 0.3], leg = true), 
        groupedlineplot(:t, :y_j_1_1, :treatment_id, ylabel = "y_1" , fillalpha = .2, lw = 2, leg = false),
        groupedlineplot(:t, :y_j_1_2, :treatment_id, ylabel = "y_2" , fillalpha = .2, lw = 2, leg = false),
        groupedlineplot(:t, :y_j_1_3, :treatment_id, ylabel = "y_3" , fillalpha = .2, lw = 2, leg = false),
        groupedlineplot(:t, :y_j_1_7, :treatment_id, ylabel = "y_7" , fillalpha = .2, lw = 2, leg = false),
        xlabel = "Time (d)", size = (800,350), bottommargin = 5mm
    )

    return plt

end

"""
Examples for plausible parameter estimates for this dataset. 
"""
function reference_param_estims()

    p_opt_control = ComponentVector(OrderedDict(zip(
        Symbol.(["spc.Z_mean_Ugent", "spc.eta_AS_juv"]), 
        [0.51, 0.73]
    )))

    p_opt_M = [
        1.374,
        0.065,
        3.65
    ]

    p_opt_G = [
        0.22,
        0.29,
        4.54
    ]

    return p_opt_control, p_opt_M, p_opt_G
end