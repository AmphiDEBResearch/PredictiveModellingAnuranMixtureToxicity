

includet(scriptsdir("ModelFitting_Discoglossus_BdJEL423_exp2.jl"))

const PMOA_FLP = "G" #"KAP"

function extract_aquatic_data(sim::AbstractDataFrame)::DataFrame
 
    aquatic = @chain sim begin
        @select(:t_exp, :larva, :wetmass_mg, :replicate, :treatment_id, :treatment_bd)
        @subset(:larva .> 0.5) # only larvae 
        groupby([:t_exp, :treatment_id, :treatment_bd]) # for every time point and treatment
        combine(
            [:wetmass_mg,:replicate] => 
            ((m, r) -> 
                (wetmass_mg = mean(m), # calculate average dry mass
                num_tadpoles = length(unique(r))) # count how many tadpoles we have
                ) => AsTable
            )
    end

    return aquatic
end

function extract_metamorph_data(sim::AbstractDataFrame, time_since_hatch)::DataFrame

    # calculate metamorphosis traits for each replicate and treatment
    metamorphs = combine(groupby(sim, [:replicate, :treatment_id, :treatment_bd])) do df
        DataFrame(
            t_exp_G42 = metamorphosis_timing(df)[1] .- age_at_birth(df) .- time_since_hatch, # timing of Gosner 42 (days since start of experiment)
            t_exp_G46 = metamorphosis_timing(df)[2] .- age_at_birth(df) .- time_since_hatch, # timing of Gosner 46 (days since start of experiment)
            dt_mt = metamorphosis_duration(df), # duration of metamorphosis (d)
            wetmass_G42_mg = wetmass_at_G42(df), # wet mass at Gosner 42 (mg)
            wetmass_G46_mg = wetmass_at_G46(df) # wet mass at Gosner 46 (mg)
        )
    # calculate averages per treatment
    end |> x -> combine(groupby(x, [:treatment_id, :treatment_bd])) do df
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

treatment_ids_bd = Dict(
    "uninfected" => 1,
    "JEL423" => 2
)

treatment_ids_chem = Dict(zip(
    [0, 10., 100.],
    collect(eachindex([0., 10., 100.]))
))

"""
    simulator(
        p::EcotoxSystems.ComponentVector; # parameters and forcings
        return_raw::Bool = false, # return raw simulation output? if false, converts output to format of the data
        param_links::NamedTuple = (ind = link_ind_params!,),
        kwargs... # additional arguments for ODE_simulator
        )


Mimicks UGent "experiment 2".
"""
function simulator_Bd_Flp(
    p::EcotoxSystems.ComponentVector; # parameters and forcings
    return_raw::Bool = false, # return raw simulation output? if false, converts output to format of the data
    param_links::NamedTuple = (ind = link_ind_params!,), 
    kwargs... # additional arguments for ODE_simulator
    )

    # parameter pre-processing

    let Z = deepcopy(p.spc.Z), inoculation_dose = deepcopy(p.glb.pathogen_inoculation_dose)

        p.glb.C_W .= 0. # reset C_W => exposure() will take care of this
        p.spc.time_since_birth = ceil(p.spc.time_since_birth) # convert time since birth to whole day
        
        # add parameter links

        p.spc.dI_max_emb = p.spc.dI_max_lrv
        p.spc.k_M_juv = p.spc.k_M_emb

        # estimate population mean of the embryonic development time

        p.spc.Z = Dirac(1.) # turn off individual variabiilty to get estimate based on popmean 
        p.spc.emb_dev_time = estimate_emb_dev_time(p) # assign embryonic development time
        p.spc.Z = Z # re-assign original zoom factor
        
         
        #try
            inner_sim(p) = @replicates AmphiDEB.ODE_simulator(
                p, 
                param_links = param_links,
                maxiters = 1e5, # default is 1e5
                statevars_init = initialize_statevars_noexposure,
                callbacks = callbacks_exp2(),
                #model = AmphiDEB.AmphiDEB_ODE_with_linear_TD!,
                kwargs...
                ) 10

            p.glb.pathogen_inoculation_dose = 0.
            sim = exposure(
                inner_sim, 
                p,
                Matrix(hcat([0, 10., 100.]...)')
            )
            
            sim[!,:treatment_id_bd] .= 1
            sim[!,:treatment_bd] .= "uninfected"

            p.glb.pathogen_inoculation_dose = inoculation_dose
            sim_infected = exposure(
                inner_sim, 
                p,
                Matrix(hcat([0, 10., 100.]...)')
            )
            sim_infected[!,:treatment_id_bd] .= 2
            sim_infected[!,:treatment_bd] .= "JEL423"

            append!(sim, sim_infected)

            # overwrite treatment id defined by exposure function,
            # by inferring a combined id for chemical and pathogen exposure
            sim[!,:treatment_id] = ["$(a)-$(b)" for (a,b) in zip(sim.treatment_id_bd,sim.treatment_id)]

            # calculate total dry mass and wet mass
            sim[!,:drymass_mg] = sim.S .+ sim.E_mt 
            sim[!,:wetmass_mg] = [calc_wetmass(r, p.spc.watercontent_larvae, p.spc.watercontent_juveniles) for r in eachrow(sim)]

            # convert simulation time to time since start of experiment
            sim[!,:t_exp] = sim.t .- p.spc.emb_dev_time

            # optionally, return the raw simulation output
            if return_raw
                return sim 
            end

            # convert raw simulation output to dataset

            sim_data = OrderedDict(
                :aquatic => extract_aquatic_data(sim),
                :metamorphs => extract_metamorph_data(sim, 0.),
                :juveniles => extract_juvenile_data(sim),
                :bdloads => extract_bdload_data(sim)
            )

            return sim_data
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


function define_defaultparams_combinedeffects()::ComponentVector

    p = ComponentVector(
        glb = ComponentVector(
            AmphiDEB.defaultparams.glb; 
            chemical_addition_time = 20., # time since hatching at which chemical is added 
            chemical_exposure_duration = 20., # duration of chemical exposure 
        ), 
        pth = AmphiDEB.defaultparams.pth,
        spc = ComponentVector(
            AmphiDEB.defaultparams.spc; 
            # auxiliary parameters
            watercontent_larvae = 0.93, 
            watercontent_juveniles = 0.85,
            time_since_birth = 15.,
            emb_dev_time = 2. 
        ))

    # setting global parameters

    p.glb.t_max = 100. # setting simulation time conservatively, for cases where metamorphosis is delayed a lot
    p.glb.pathogen_inoculation_time = 47. # time since hatching at which pathogen is added
    p.glb.pathogen_inoculation_dose = 1e6 # number of zoospores added at specified time point
    p.glb.medium_renewals = [31. + (8/24)] # time of medium_renewal (pathogen removal), i.e. pathogen renewal; exposure lasted 8h
    p.glb.dX_in = [1e10, 1e10] # ad libitum feeding conditions

    p.spc.Z = truncated(Normal(1, 0.1), 0, Inf)
    # propagation of zoom factor to H_j1 is turned off => we want variability in the transition to metamorphs
    p.spc.propagate_zoom.H_j1 = 0.

    # adding point estimates as defaults

    posterior_summary_larvalfit = CSV.read(datadir("sims", SAVETAG_LARVALFIT, "posterior_summary.csv"), DataFrame)
    posterior_summary_juvenilefit = CSV.read(datadir("sims", SAVETAG_JUVENILEFIT, "posterior_summary.csv"), DataFrame)
    posterior_summary_TKTD = CSV.read(datadir("sims", "$(SAVETAG_TKTDFIT)_$(PMOA_FLP)", "posterior_summary.csv"), DataFrame)

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

    p.spc.KD .= 0.
    p.spc.B .= 2.
    p.spc.E .= 1e10

    @info "Overwriting default values of $(posterior_summary_TKTD.param)"
    for (label,value) in zip(posterior_summary_TKTD.param, posterior_summary_TKTD.best_fit)
        assign_value_by_label!(p, label, value)
    end

    p.spc.k_M_juv = p.spc.k_M_emb
    p.spc.X_emb_int = 1. # ≈ initial dry mass of an egg (mg)
    p.spc.emb_dev_time = estimate_emb_dev_time(p)

    p.spc.E_P .= 1e10
    p.spc.B_P .= Inf

    return p
end


function setup_simulations_combinedeffects(
    pmoa::AbstractString; 
    sigma_factor = 1., 
    loss_functions = loss_mse
    )

    if !isdir(datadir("sims", SAVETAG))
        mkdir(datadir("sims", SAVETAG))
    end

    pmoa_idx = findfirst(x -> x == pmoa, PMOAS)

    # set up logger
    # FIXME: this does not save the log file to the intended directory, because the setup_modelfit does not know the composite "savetag" 
    global io = open(datadir("sims", SAVETAG, "log.txt"), "w+")
    global logger = SimpleLogger(io)

    data = load_data_exp2_Bd()

    f = ModelFit( 
        prior = Prior(
            #"spc.eta_AS_juv" => Truncated(Normal(0.5, 0.5), 0, 1),
            "pth.gamma" => Truncated(Normal(0.1, 0.1), 0, Inf),
            "pth.eta" => Truncated(Normal(100, 10), 0, Inf),
            "pth.sigma1" => Truncated(Normal(0.016, 0.016), 0, Inf),
            "spc.Chi" => Hyperdist(
                σ -> LogNormal(σ^2, σ),
                Truncated(Normal(1, 1), 0, Inf)
            ),
            "spc.E_P[$(pmoa_idx)]" => Truncated(Normal(1e3, 1e3), 1e1, 1e5),
            "spc.B_P[$(pmoa_idx)]" => Truncated(Normal(2, 2), 2, 4)
            ),
        defaultparams = define_defaultparams_combinedeffects(), 
        simulator = simulator_Bd_Flp, 
        data = data, 
        response_vars = [
            [:load_mean, :load_var, :load_skew],
            [:y_W2_mg], 
        ],
        data_weights = [
            [2.0, 1.0, 0.5],
            [1.0]
        ],
        grouping_vars = [
            Symbol[], 
            [:treatment_id]
            ], 
        time_resolved = [false, false],
        time_var = :t_exp,
        plot_data = plot_data_exp2_Bd,
        loss_functions = loss_functions
    )

    return f

end