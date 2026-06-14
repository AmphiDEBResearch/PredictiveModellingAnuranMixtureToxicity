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