# functions used to fit 2,4D-models
# the functions defined here perform the model fitting for larval and juvenile life stages

# use fitting setup from M1, only overwrite parts which are different
using Revise
includet(scriptsdir("ModelFitting_Discoglossus_24D_UCLM.jl"))


"""
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
            gen_ind_params = x -> x, # skip generation of individual-level parameters - we already have them
            model = AmphiDEB.M2_complete_ODE_with_loglogistic_TD! # use AmphiDEB M2
            )
        
        # re-set global parameters to actual values
        p.glb.t_max = t_max
        p.glb.C_W .= C_W

        # retrieve the final state of the simulated experiment
        u0 = sim.u[end].ind

        # return the final state as initial state of the actual simulation, 
        # re-setting global states
        return ComponentVector(
            glb = AmphiDEB.initialize_global_statevars(p),
            ind = u0
        )
    end

end

"""
Simulate toxicity test with *D. galganoi* exposed to 2,4D.
"""
function simulator_M2(
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
                model = AmphiDEB.M2_complete_ODE_with_loglogistic_TD!, # use AmphiDEB M2
                kwargs...
                ) 10

            # epxosure() runs inner_sim for each treatment and collects the results
    
            sim = exposure(
                inner_sim, 
                p,
                Matrix(hcat([0, 0.03, 0.3, 3., 30.]...)')
            )

            # calculate total dry mass and wet mass
            sim[!,:drymass_mg] = sim.S .+ sim.E_mt 
            sim[!,:wetmass_mg] = [calc_wetmass(r, p.spc.watercontent_larvae, p.spc.watercontent_juveniles) for r in eachrow(sim)]
            
            # convert simulation time to time since start of experiment
            sim[!,:t_exp] = sim.t

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


"""
Set up AmphiDEB-M2 fitting to UCLM discoglossus 2,4-D data.  
"""
function setup_modelfit_M2(pmoa::AbstractString)

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

    # use M2 larval parameters
    defparams = define_defaultparams()
    assign_values_from_file!(
        defparams, 
        datadir("sims", "input", "Discoglossus_larvae_M2_var_delta_k_M_mt", "posterior_summary.csv"), 
        exceptions = OrderedDict(
            "spc.Z" => (p,label,value) -> p.spc.Z = truncated(Normal(1,value), 0, Inf)
        )
    )

    f = ModelFit( 
        prior = Prior(
            "spc.KD[1,$(pmoa_idx)]" => Truncated(Normal(1, 0.1), 0.001, 1), 
            "spc.E[1,$(pmoa_idx)]" => Truncated(Normal(15, 15), 0.03, Inf),
            "spc.B[1,$(pmoa_idx)]" => Truncated(Normal(2., 4), 0.5, Inf)
            ),
        defaultparams = defparams, 
        simulator = simulator_M2, 
        data = data, 
        response_vars = [
            [:wetmass_mg, :num_tadpoles], 
            [:t_exp_G42, :t_exp_G46, :wetmass_G42_mg, :wetmass_G46_mg]
        ],
        data_weights = [
            [1., 1.],
            Float64[]
        ],
        grouping_vars = [[:C_W_1], [:C_W_1]], 
        time_resolved = [true, false],
        time_var = :t_exp,
        plot_data = plot_data,
        loss_functions = loss_mse_logtransform
    )

    return f

end