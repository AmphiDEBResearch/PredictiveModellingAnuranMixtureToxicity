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
import AmphiDEB: ComponentVector

using OrdinaryDiffEq

using Revise
includet(srcdir("traits.jl"))
include(scriptsdir("ModelFitting_Discoglossus_24D_UCLM.jl"))
include(scriptsdir("ModelFitting_Discoglossus_BdJEL423_exp2.jl"))
includet(srcdir("utils.jl"))

const TREATMENT_MAP = OrderedDict(
   1 => "uninfected | 0.0", 
   2 => "uninfected | 0.03", 
   3 => "uninfected | 0.3", 
   4 => "JEL423 | 0.0", 
   5 => "JEL423 | 0.03", 
   6 => "JEL423 | 0.3"
)

function load_data_exp2_all(;
    paths::OrderedDict = OrderedDict(
        :juveniles_raw => [datadir("exp_raw", "UGent", "exp2", "juveniles.csv"), 1],
        :bdloads => [datadir("exp_raw", "UGent", "exp2", "bdloads_JEL423.csv"), 1],
    ))

    data = OrderedDict()

    for (key,info) in pairs(paths)
        path, header = info
        data[key] = CSV.read(path, DataFrame, header = header) 
        data[key][!,:treatment_id] = [reverse(TREATMENT_MAP)[x] for x in data[key].treatment_combination]
    end

    # convert mass measurements to mg
    data[:juveniles_raw][!,:wetmass_mg] = data[:juveniles_raw].weight_g * 1e3
    # drop original column
    select!(data[:juveniles_raw], Not(:weight_g))
    rename!(data[:juveniles_raw], :y_weight_g => :y_wetmass_mg)

    # bdload column is not needed in juvenile data, 
    # has to be removed because missing values can cause issues
    select!(data[:juveniles_raw], Not(:bdload))

    data[:juveniles_aggregated] = aggregate_juveniles(data[:juveniles_raw])

    # for the Bd load data, we don't need the control
    data[:bdloads] = @subset(data[:bdloads], :treatment_bd .== "JEL423")

    # return re-ordered data
    return OrderedDict(
        :juveniles_raw => data[:juveniles_raw],
        :juveniles_aggregated => data[:juveniles_aggregated], 
        :bdloads => data[:bdloads]
    )

end

function plot_data_exp2_all()

    data = load_data_exp2_all()

    plt1 = @df @subset(data[:juveniles_raw], :treatment_bd .== "uninfected") groupedviolin(
        string.(:t_since_mm), :wetmass_mg, group = :pretreatment_24D, 
        side = :left,
        legend = :topleft, legendtitle = "2,4D-pretreatment \n (mg/L)", legendtitlefontsize = 8,
        palette = palette([:cyan, :purple], 3), fillalpha = .5,
        xlabel = "Time since G42 (d)", ylabel = "Wet mass (mg)", 
        ylim = (100,500), 
        title = "Uninfected",
    )
    
    plt2 = @df @subset(data[:juveniles_raw], :treatment_bd .== "uninfected") groupedviolin(
        string.(:t_since_mm), :y_wetmass_mg, group = :pretreatment_24D, 
        side = :left,
        legend = false, 
        palette = palette([:cyan, :purple], 3), fillalpha = .5,
        xlabel = "Time since G42 (d)", ylabel = "Control-normalized \n wet mass (-)", 
        ylim = (0,2),
        title = "Uninfected"
    )

    plt3 = @df @subset(data[:juveniles_raw], :treatment_bd .== "JEL423") groupedviolin(
        string.(:t_since_mm), :wetmass_mg, group = :pretreatment_24D, 
        side = :left,
        legend = :topleft, legendtitle = "2,4D-pretreatment \n (mg/L)", legendtitlefontsize = 8,
        palette = palette([:cyan, :purple], 3), fillalpha = .5,
        xlabel = "Time since G42 (d)", ylabel = "Wet mass (mg)", 
        ylim = (100,500),
        title = "BdJEL423"
    )
    
    plt4 = @df @subset(data[:juveniles_raw], :treatment_bd .== "JEL423") groupedviolin(
        string.(:t_since_mm), :y_wetmass_mg, group = :pretreatment_24D, 
        side = :left,
        legend = false, 
        palette = palette([:cyan, :purple], 3), fillalpha = .5,
        xlabel = "Time since G42 (d)", ylabel = "Control-normalized \n wet mass (-)", 
        ylim = (0,2), 
        title = "BdJEL423"
    )

    plot(
        plt1, plt2, plt3, plt4, 
        layout = (2,2), size = (800,700), 
        bottommargin = 5mm, leftmargin = 5mm,
        legend_background_color = :transparent,
        foreground_color_legend = nothing
        )

end

function plot_sims_exp2_all!(plt, sims::AbstractVector)

    juveniles = vcat(map(x->x[:juveniles_raw], sims)...) |> clean

    @df @subset(juveniles, :treatment_id .< 4) groupedviolin!(    
        plt,
        string.(:t_since_mm), :wetmass_mg,
        group = :treatment_id,
        side = :right,
        color = :steelblue, 
        linecolor = :steelblue, 
        label = "", 
        subplot = 1,
        fillalpha = .2
    )
    @df @subset(juveniles, :treatment_id .< 4) groupedviolin!(    
        plt,
        string.(:t_since_mm), :y_wetmass_mg,
        group = :treatment_id,
        side = :right,
        color = :steelblue, 
        linecolor = :steelblue, 
        label = "", 
        subplot = 2,
        fillalpha = .2
        )

    @df @subset(juveniles, :treatment_id .>= 4) groupedviolin!(    
        plt,
        string.(:t_since_mm), :wetmass_mg,
        group = :treatment_id,
        side = :right,
        color = :steelblue, 
        linecolor = :steelblue, 
        label = "", 
        subplot = 3,
        fillalpha = .2
        )

    @df @subset(juveniles, :treatment_id .>= 4) groupedviolin!(    
        plt,
        string.(:t_since_mm), :y_wetmass_mg,
        group = :treatment_id,
        side = :right,
        color = :steelblue, 
        linecolor = :steelblue, 
        label = "", 
        subplot = 4,
        fillalpha = .2
        )

    return plt
end

"""
Set up defaultparams for model validation against 2,4-D/BdJEL423 data. 

- posterior_summary_larvalfit: larval/metamorph DEB parameters
- posterior_summary_juvenilefit: juvenilme/adult DEB parameters
- posterior_summary_TKTDfit: 2,4-D parameters, incl. correction of DEB parameters (zoom factor, growth efficiency)
"""
function define_defaultparams_exp2_all(;
    posterior_summary_larvalfit::String = datadir("sims", SAVETAG_LARVALFIT, "posterior_summary.csv"),
    posterior_summary_juvenilefit::String = datadir("sims", SAVETAG_JUVENILEFIT, "posterior_summary.csv"),
    posterior_summary_correction_factors::String = datadir("sims", SAVETAG_CORRECTION_FACTORS, "posterior_summary.csv"),
    posterior_summary_TKTDfit = datadir("sims", SAVETAG_TKTDFIT, "posterior_summary.csv"),
    posterior_summary_Bdfit::String = datadir("sims", SAVETAG_BDFIT, "posterior_summary.csv"),
    )::ComponentVector

    p = ComponentVector(
        glb = ComponentVector(
            AmphiDEB.defaultparams.glb; 
            chemical_addition_time = 26., # time since hatching at which chemical is added 
            chemical_exposure_duration = 5., # duration of chemical exposure 
        ), 
        pth = AmphiDEB.defaultparams.pth,
        spc = ComponentVector(
            AmphiDEB.defaultparams.spc; 
            # auxiliary parameters
            watercontent_larvae = 0.93, 
            watercontent_juveniles = 0.85,
            time_since_birth = 15.,
            emb_dev_time = 2.,
            Z_UCLM = truncated(Normal(1, 0.1), 0, Inf), # Z estimated from UCLM data - will be assigned from file later
            Z_mean_UGent = 1., # size correction factor for UGent vs UCLM data
            H_j1_UCLM = 10, # H_j1 estimate from UCLM data
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

    p.spc.KD .= 0.
    p.spc.B .= 2.
    p.spc.E .= 1e10

    # adding point estimates as defaults

    assign_values_from_file!(
        p, 
        posterior_summary_larvalfit,
        exceptions = OrderedDict(
            "spc.Z" => (p,label,value) -> p.spc.Z = Truncated(Normal(1, value), 0, 1)
            )
    )
    p.spc.H_j1_UCLM = p.spc.H_j1

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
    p.spc.Z = Truncated(Normal(p.spc.Z_mean_UGent, p.spc.Z_mean_UGent * p.spc.Z.untruncated.σ), 0, Inf)
    p.spc.H_j1 = p.spc.H_j1_UCLM * p.spc.Z_mean_UGent
    
    assign_values_from_file!(
        p, 
        posterior_summary_TKTDfit,
        exceptions = OrderedDict()
        )

    assign_values_from_file!(
        p, 
        posterior_summary_Bdfit, 
        exceptions = OrderedDict(
            "spc.Chi" => (p,label,value) -> p.spc.Chi = LogNormal(value^2, value),
        )
    )

    p.spc.k_M_juv = p.spc.k_M_emb # link juv/ad somatic maintenance to emb/lrv/met somatic maintenance
    p.spc.X_emb_int = 1. # ≈ initial dry mass of an egg (mg)
    p.spc.emb_dev_time = estimate_emb_dev_time(p) # get an estimate of the population mean of embryonic development toimes

    return p
end

"""
Mimicks UGent "experiment 2", including combined 2,4-D/Bd treatments for purpose of validation.
"""
function simulator_exp2_validation(
    p::EcotoxSystems.ComponentVector; # parameters and forcings
    return_raw::Bool = false, # return raw simulation output? if false, converts output to format of the data
    param_links::NamedTuple = (ind = link_ind_params!,), 
    kwargs... # additional arguments for ODE_simulator
    )

    # parameter pre-processing

    let Z = deepcopy(p.spc.Z), inoculation_dose = deepcopy(p.glb.pathogen_inoculation_dose)
        
        # when we uncomment the code block below, we propagate uncertainties in DEB parameters, 
        # but assume that DEB parameters are uncorrelated with pathogen and TKTD parameters (which might be ok)
        # by default, this is disengaged, which means that we take the point estimates of DEB parameters
        ## assign a posterior sample from larval parameters
        #posterior_sample!(
        #    p, 
        #    accepted_larvalfit, 
        #    exceptions = OrderedDict("spc.Z" => (p,label,value) -> p.spc.Z = Truncated(Normal(1,value),0,Inf))
        #    )
        #
        ## assign a posterior sample from juvenile parameters
        #posterior_sample!(p, accepted_juvenilefit)

        p.glb.C_W .= 0. # reset C_W => exposure() will take care of this
        p.spc.time_since_birth = ceil(p.spc.time_since_birth) # convert time since birth to whole day
        
        # add parameter links

        p.spc.dI_max_emb = p.spc.dI_max_lrv
        p.spc.k_M_juv = p.spc.k_M_emb

        # estimate population mean of the embryonic development time

        p.spc.Z = Dirac(1.) # turn off individual variabiilty to get estimate based on popmean 
        p.spc.emb_dev_time = estimate_emb_dev_time(p) # assign embryonic development time
        p.spc.Z = Z # re-assign original zoom factor
        
        # try

            # "inner simulator" is the function called for each treatment
            inner_sim(p) = @replicates AmphiDEB.ODE_simulator(
                p, 
                param_links = param_links,
                maxiters = 1e5, # default is 1e5
                statevars_init = initialize_statevars_noexposure,
                callbacks = callbacks_exp2(),
                kwargs...
                ) 10

            ## simulate chemical exposure without infection

            p.glb.pathogen_inoculation_dose = 0.
            sim_uninfected = exposure(
                inner_sim, 
                p,
                Matrix(hcat([0, 0.03, 0.3]...)')
            )
            
            sim_uninfected[!,:treatment_bd] .= "uninfected"

            ## simulate chemical exposure + infection

            p.glb.pathogen_inoculation_dose = inoculation_dose
            sim_infected = exposure(
                inner_sim, 
                p,
                Matrix(hcat([0, 0.03, 0.3]...)')
            )

            sim_infected[!,:treatment_id] .+= maximum(sim_uninfected.treatment_id)
            sim_infected[!,:treatment_bd] .= "JEL423"

            sim = vcat(sim_uninfected, sim_infected)

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
