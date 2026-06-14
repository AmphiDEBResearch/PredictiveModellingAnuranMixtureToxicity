
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
        #try 
            p.glb.C_W .= 0.
        #catch e 
        #    println(p.glb.C_W)
        #    error(e)
        #end

        # simulate the model until the start of the experiment
            sim = AmphiDEB.ODE_simulator(
                p,
                returntype = EcotoxSystems.odesol, # directly return the ODE solution object - we don't need a DataFrame
                gen_ind_params = x -> x # skip generation of individual-level parameters - we already have them
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

function link_ind_params!(ind::AmphiDEB.ComponentVector)::Nothing

    # defining links between parameters 
    # this is mostly relevant for parameters which may be subject to the zoom factor
    # this is handled by EcotoxSystems.jl for each simulated individual 

    ind.dI_max_emb = ind.dI_max_lrv # ingestion rate for embryos assumed to be same as for larave

    # FIXME: this does not appear to have the desired effect
    # the issue is currently fixed by adding expression to simulator; works for now because k_M is not linked to Z
    # this should be fixed though and a test added in EcotoxSystems.jl and AmphiDEB.jl
    ind.k_M_juv = ind.k_M_emb # somatic maintenace rate is assumed to remain constant across life stages

    ind.k_J_emb = (1 - ind.kappa_emb) / ind.kappa_emb * ind.k_M_emb # maturity maintenance is linked to somatic (assuming same cumulative investment in both branches)
    #ind.k_J_juv = (1-ind.kappa_juv)/ind.kappa_juv * ind.k_M_juv
    ind.kappa_juv = ind.kappa_emb

    return nothing
end

"""
Simulate toxicity test with *D. galganoi* exposed to mixture of 2,4-D and FLP.
"""
function simulator_UCLM_mixture(
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

        # "inner simulator" is the function called for each treatment
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
            EXPOSURE_MATRIX
        )

        # calculate total dry mass and wet mass
        sim[!,:drymass_mg] = sim.S .+ sim.E_mt 
        sim[!,:wetmass_mg] = [calc_wetmass(r, p.spc.watercontent_larvae, p.spc.watercontent_juveniles) for r in eachrow(sim)]
        
        # convert simulation time to time since start of experiment
        sim[!,:t_exp] = sim.t

        # add Bd treatment column - only needed to have a generically applicable extract_aquatic_data(), extract_metamorph_data()
        sim[!,:treatment_bd] .= "uninfected"

        # optionally, return the raw simulation output
        if return_raw
            return sim 
        end

        # convert raw simulation output to dataset

        sim_data = OrderedDict(
            :aquatic => extract_aquatic_data(sim),
            :metamorphs => extract_metamorph_data(sim, 0.),
        )

        return sim_data
    
    end
end


function plot_sims_UCLM_mix_growth!(plt, sims::AbstractVector)
        
    sims_aquatic = filter(x->!isnothing(x), sims) |> 
    x -> map(d -> d[:aquatic], x) |> 
    x -> vcat(x...) |> 
    clean

    lcol = leftcol()
    trow = toprow()
    single_substance = vcat(lcol, trow)

    for i in sort(unique(sims_aquatic.treatment_id))

        df = @subset(sims_aquatic, :treatment_id .== i)

        @df df lineplot!(
            plt, 
            :t_exp, :wetmass_mg, 
            subplot = i,  lw = 2, 
            label = in(i, single_substance) ? "Retrodicted" : "Predicted", 
            color = in(i, single_substance) ? :gray : :steelblue, 
            ylabel = i in lcol ? "Wet mass \n (mg)" : "",
            xlabel = i == maximum(sims_aquatic.treatment_id)-1 ? "Time since start of experiment (d)" : "",
            fillalpha = .2, 
            leftmargin = 5mm
            )

    end
    lineplot!([], [], color = :steelblue, subplot = 1, label = "Predicted", lw = 2)
    return plt
end

# this is outdated -- we will use fract_tadpoles to correct for variable number of survivors in the data
function plot_sims_UCLM_mix_numtadpoles!(plt, sims::AbstractVector)
        
    sims_aquatic = filter(x->!isnothing(x), sims) |> 
    x -> map(d -> d[:aquatic], x) |> 
    x -> vcat(x...) |> 
    clean

    lcol = leftcol()
    trow = toprow()
    single_substance = vcat(lcol, trow)

    for i in sort(unique(sims_aquatic.treatment_id))

        df = @subset(sims_aquatic, :treatment_id .== i)

        @df df lineplot!(
            plt, 
            :t_exp, :num_tadpoles, 
            subplot = i,  lw = 2, 
            label = in(i, single_substance) ? "Retrodicted" : "Predicted", 
            color = in(i, single_substance) ? :gray : :steelblue, 
            ylabel = i in lcol ? "Wet mass \n (mg)" : "",
            xlabel = i == maximum(sims_aquatic.treatment_id)-1 ? "Time since start of experiment (d)" : "",
            fillalpha = .2, 
            leftmargin = 5mm, 
            fillstyle = ://
            )

    end
    lineplot!([], [], color = :steelblue, subplot = 1, label = "Predicted", lw = 2)
    return plt
end

function plot_sims_UCLM_mix_fracttadpoles!(plt, sims::AbstractVector)
        
    sims_aquatic = filter(x->!isnothing(x), sims) |> 
    x -> map(d -> d[:aquatic], x) |> 
    x -> vcat(x...) |> 
    clean

    lcol = leftcol()
    trow = toprow()
    single_substance = vcat(lcol, trow)

    for i in sort(unique(sims_aquatic.treatment_id))

        df = @subset(sims_aquatic, :treatment_id .== i)

        @df df lineplot!(
            plt, 
            :t_exp, :fract_tadpoles, 
            subplot = i,  lw = 2, 
            label = in(i, single_substance) ? "Retrodicted" : "Predicted", 
            color = in(i, single_substance) ? :gray : :steelblue, 
            fillalpha = .2, 
            leftmargin = 5mm,
            ylabel = i in lcol ? "Fraction of \n tadpoles" : "",
            xlabel = i == maximum(sims_aquatic.treatment_id)-1 ? "Time since start of experiment (d)" : "",
            )

    end
    lineplot!([], [], color = :steelblue, subplot = 1, label = "Predicted", lw = 2)
    return plt
end

function plot_sims_UCLM_mix_metamorphs!(plt, sims::AbstractVector) 

    sim_metamorphs = EcotoxModelFitting.extract_simkey(sims, :metamorphs)
    leftjoin!(sim_metamorphs, TREATMENT_IDS_MIX, on = :treatment_id)

    num_flp_treatments = length(unique(sim_metamorphs.F_ppm))

    for (i,F) in enumerate(unique(sim_metamorphs.F_ppm))

        df = @subset(sim_metamorphs, :F_ppm .== F)
        sort!(df, :D_ppm)

        df_retro = @subset(df, (:D_ppm .== 0) .& (:F_ppm .== 0))
        df_pred = @subset(df, (:D_ppm .> 0) .| (:F_ppm .> 0))

        is_mixture_treatment = [(F>0) && (D>0) for D in unique(df.D_ppm)]
        c = [ismix ? :steelblue : :lightgray for ismix in is_mixture_treatment]
      
        @df df_retro violin!(
            plt, subplot = i, 
            string.(:D_ppm), :t_exp_G46, 
            side = :right, 
            color = :gray, 
            fillalpha = .25, 
            label = i == 1 ? "Retrodicted" : "", 
            fillstyle = ://, 
            trim = true
        )

        @df df_pred violin!(
            plt, subplot = i, 
            string.(:D_ppm), :t_exp_G46, 
            side = :right, 
            color = :steelblue, 
            fillalpha = .25, 
            label = i == 1 ? "Predicted" : "",
            trim = true
        )


        @df df_retro violin!(
            plt, subplot = i+num_flp_treatments, 
            string.(:D_ppm), :wetmass_G46_mg,
            side = :right, 
            color = :gray, 
            fillalpha = .25, 
            label = "",
            fillstyle = ://, 
            xlabel = "2,4-D (mg/L)", 
            bottommargin = 5mm
        )

        @df df_pred violin!(
            plt, subplot = i+num_flp_treatments, 
            string.(:D_ppm), :wetmass_G46_mg, 
            side = :right, 
            color = :steelblue, 
            fillalpha = .25, 
            label = ""
        )


 
    end

    return plt

end
