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
includet(srcdir("traits.jl"))
includet(srcdir("utils.jl"))
include(scriptsdir("ModelFitting_Discoglossus_24D_UCLM.jl"))

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

function define_defaultparams_UCLM_mix(;
    posterior_summary_larvalfit=datadir("sims", SAVETAG_LARVALFIT, "posterior_summary.csv"),
    posterior_summary_juvenilefit=datadir("sims", SAVETAG_JUVENILEFIT, "posterior_summary.csv"),
    posterior_summary_24dfit=datadir("sims", SAVETAG_24DFIT, "posterior_summary.csv"),
    posterior_summary_flpfit=datadir("sims", SAVETAG_FLPFIT, "posterior_summary.csv")
    )

    # we need to reconstruct the parameter vector, because ComponentVector's cannot change their size dynamically
    p = ComponentVector(
        glb = ComponentVector(
            t_max = 56.0,
            N0 = 1.0,
            dX_in = [20.0, 20.0],
            k_V = [0.0, 0.0],
            V_patch = [1.0, 1.0],
            T = 293.15,
            C_W = [0.0, 0.0],
            pathogen_inoculation_dose = 0.0,
            pathogen_inoculation_time = 30.0,
            medium_renewals = [0.0]
            ),
        pth = AmphiDEB.defaultparams.pth,
        spc = ComponentVector(
            Z = Dirac(1.0),
            propagate_zoom = (
                dI_max_emb = 0.3333333333333333,
                dI_max_lrv = 0.3333333333333333,
                dI_max_juv = 0.3333333333333333,
                X_emb_int = 1.0,
                H_j1 = 1.0,
                H_p = 1.0,
                K_X_lrv = 0.3333333333333333,
                K_X_juv = 0.3333333333333333
                ),
            X_emb_int = 1,
            K_X_lrv = 1.0,
            K_X_juv = 1.0,
            dI_max_emb = 1,
            dI_max_lrv = 1,
            dI_max_juv = 1,
            kappa_emb = 0.8,
            kappa_juv = 0.8,
            gamma = 0.5,
            eta_IA = 0.54,
            eta_AS_emb = 0.4,
            eta_AS_juv = 0.4,
            eta_AR = 0.95,
            eta_SA = 0.8,
            k_M_emb = 0.11,
            k_M_juv = 0.11,
            delta_k_M_mt = 1.0,
            k_J_emb = 0.027,
            k_J_juv = 0.027,
            H_j1 = 1,
            H_p = 55.0,
            delta_E = 1.0,
            T_A = 8000.0,
            T_ref = 293.15,
            b_T = 40.0,
            fb_G = 0.0,
            h_b = 0.0,
            
            KD = [0. 0. 0. 0. 0. 0. 0.; 0. 0. 0. 0. 0. 0. 0.], # k_D - value per PMoA (G,M,A,R,H,kap) and stressor (1 row = 1 stressor)
            B = [2. 2. 2. 2. 2. 2. 2.; 2. 2. 2. 2. 2. 2. 2.], # slope parameters
            E = [1e10 1e10 1e10 1e10 1e10 1e10 1e10; 1e10 1e10 1e10 1e10 1e10 1e10 1e10], # sensitivity parameters (thresholds)
            KD_h = [0.; 0.], # k_D - value for GUTS-Sd module (1 row = 1 stressor)
            E_h = [1e10; 1e10], # sensitivity parameter (threshold) for GUTS-SD module
            B_h = [1.; 1.], # slope parameter for GUTS-SD module 
            C_h = [1.; 1.], # proportionality constant to convert relative response to hazard rate 
            
            S_rel_crit = 0.66,
            h_S = 0.6,
            a_max = Truncated(Normal(5475.0, 547.5), 0.0, Inf),
            tau_R = 365.0,
            Chi = LogNormal(1.0, 1.0),
            E_P = [Inf, Inf, Inf, Inf],
            B_P = [2.0, 2.0, 2.0, 2.0],
                
            # auxiliary parameters
            watercontent_larvae = 0.93,
            watercontent_juveniles = 0.85,
            time_since_birth = 15.,
            emb_dev_time = 2.
            )
    )

    # setting global parameters

    p.glb.t_max = 100. # setting simulation time conservatively, for cases where metamorphosis is delayed a lot
    p.glb.pathogen_inoculation_time = Inf # no pathogen inoculation
    p.glb.dX_in = [1e10, 1e10] # ad libitum feeding conditions

    p.spc.Z = truncated(Normal(1, 0.1), 0, Inf)
    # propagation of zoom factor to H_j1 is turned off => we want variability in the transition to metamorphs
    p.spc.propagate_zoom.H_j1 = 0.

    # adding point estimates from calibrations as defaults

    p.spc.KD .= 0.
    p.spc.B .= 2.
    p.spc.E .= 1e10

    assign_values_from_file!(
        p,
        posterior_summary_larvalfit;
        exceptions = OrderedDict(
            "spc.Z" => (p, label, value) -> p.spc.Z = truncated(Normal(1, value), 0, Inf)
        )
    )

    assign_values_from_file!(
        p, 
        posterior_summary_juvenilefit; 
        exceptions = OrderedDict()
        )
    
    assign_values_from_file!(
        p, 
        posterior_summary_24dfit; 
        exceptions = OrderedDict()
        )
    
    # for flp, we need to bump the stressor idx
    postsum_flp = CSV.read(posterior_summary_flpfit, DataFrame) # get the parameters
    set_stressor_idx!(postsum_flp, 2) # bump stressor index
    CSV.write("postsum_flp.csv", postsum_flp) # save modified version to file, 
    assign_values_from_file!(p, "postsum_flp.csv"; exceptions = OrderedDict()) # so that we can use this function, 
    #rm("postsum_flp.csv") # remove the file

    p.spc.k_M_juv = p.spc.k_M_emb
    p.spc.X_emb_int = 1. # ≈ initial dry mass of an egg (mg)

    #p.spc.emb_dev_time = estimate_emb_dev_time(p)

    return p
end

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
