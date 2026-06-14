using AmphiDEB, EcotoxModelFitting, EcotoxSystems
using Logging
using Suppressor
using LaTeXStrings
using StatsBase
using Plots, StatsPlots, Plots.Measures
using CSV

import EcotoxModelFitting: Hyperdist

# re-imports from AmphiDEB
using ComponentArrays
using Distributions
using DataStructures

include(srcdir("traits.jl"))
include(srcdir("utils.jl"))

# relative path to the reference parameters
# expecting project structure ./data/sims/input/...
const SAVETAG_LARVALFIT = "input/Discoglossus_larvae" 
const SAVETAG_JUVENILEFIT = "input/Discoglossus_juveniles"
const SAVETAG = "SparlingFellers"
# initial number of tadpoles per experimental unit (tank/aquarium...)
# this is only relevant if we consider biological variability in the analysis (not the default)
const INIT_NUM_TADPOLES = 6 

"""
Select a species / substance from the data.
"""
function subset_data(
    df::AbstractDataFrame, 
    species::AbstractString, 
    substance::AbstractString
    )

    subdf = @subset(df, :species .== species)
    if substance == "C_cpf"
        subdf = @subset(subdf, :C_end .== 0)
    elseif substance == "C_end"
        subdf = @subset(subdf, :C_cpf .== 0)
    else 
        error("Unknown substance $(substance)")
    end

    return subdf
end

function subset_data(
    data::OrderedDict, 
    species::AbstractString,
    substance::AbstractString
    )

    subdata = OrderedDict()

    for (key,value) in pairs(data)
        subdata[key] = subdf(data[value], species, substance)
    end

    return subdata
end

"""
Load aquatic data.
"""
function load_data_sparlingfellers_aquatic(species, substance)
    
    aquatic = CSV.read(datadir("exp_raw", "SparlingFellers", "growth.csv"), DataFrame, missingstring = ["NA"]) |> x-> subset_data(x, species, substance)
    aquatic[!,:wetmass_mg] = 1e3 * aquatic.wetmass_g
    rename!(aquatic, substance => :C_W_1)

    select!(aquatic, :species, :t_exp, :C_W_1, :wetmass_mg)
    dropmissing!(aquatic)

    return aquatic
end


"""
Load metamorph data
"""
function load_data_sparlingfellers_metamorphs(species, substance)
    metamorphs = CSV.read(datadir("exp_raw", "SparlingFellers", "metamorphs.csv"), DataFrame, missingstring = ["NA"]) |> 
    x -> subset_data(x, species, substance)
    
    metamorphs[!,:wetmass_G42_mg] = 1e3 * metamorphs.wetmass_G42_g
    rename!(metamorphs, substance => :C_W_1)

    select!(metamorphs, :species, :C_W_1, :t_exp_G42, :wetmass_G42_mg)
    dropmissing!(metamorphs)

    return metamorphs
end

"""
Load all endpoints + combine into dict
"""
function load_data_sparlingfellers(;
    species::AbstractString, 
    substance::AbstractString,
    controls_only = false
    )

    aquatic = load_data_sparlingfellers_aquatic(species, substance) 
    metamorphs = load_data_sparlingfellers_metamorphs(species, substance) 
    
    if controls_only
        aquatic = @subset(aquatic, :C_W_1 .== 0)
        metamorphs = @subset(metamorphs, :C_W_1 .== 0)
    end
    
    return OrderedDict(
        :aquatic => aquatic,
        :metamorphs => metamorphs
    )
end

"""
Plot control data
"""
function plot_data_controls(data)

    for (key,val) in pairs(data)
        data[key] = @subset(val, :C_W_1 .== 0)
    end

    plt_controls = @df data[:aquatic] scatter(:t_exp, :wetmass_mg, color = :black, label = "Observed")
    vline!(plt_controls, data[:metamorphs].t_exp_G42, color = :black , linestyle = :dash, label = "")
    annotate!(plt_controls, data[:metamorphs].t_exp_G42 .* 0.95, maximum(data[:aquatic].wetmass_mg), Plots.text("G42", 12))

    return plt_controls

end

"""
Plot aquatic data (tadpole mass over time)
"""
function plot_data_aquatic(
    df::DataFrame;
    kwargs...
    )

    plt = @df df scatter(
        :t_exp, :wetmass_mg, group = :C_W_1;
        color = :black,
        layout = (1,length(unique(:C_W_1))), 
        title = hcat(unique(:C_W_1)...), 
        marker = true, size = (1000,300), leg = false,
        leftmargin = 5mm, bottommargin = 5mm, 
        titlefontsize = 10,
        link = :both,
        kwargs...
        )

    plot!(plt, subplot = 1, ylabel = "Wet mass (mg)")
    plot!(plt, xlabel = "Time (d)")
    
    return plt
    
end


"""
Plot metamorph data (mass/time over concentration)
"""
function plot_data_metamorphs(
    df::DataFrame;
    kwargs...
    )


    plt = @df df plot(
        scatter(
            string.(:C_W_1), :t_exp_G42;
            color = :black, label = "Observed", 
            kwargs...
        ),
        scatter(
            string.(:C_W_1), :wetmass_G46_mg;
            color = :black, label = "", 
            kwargs...
        )
    )
        
    plot!(plt, subplot = 1, ylabel = "Time to G42 (d)")
    plot!(plt, subplot = 2, ylabel = "Wet mass \n at metamorphosis (mg)")
    plot!(plt, xlabel = "Exposure (mg/L)")
    
    return plt
    
end


plot_data_aquatic(data::AbstractDict; kwargs...) = plot_data_aquatic(data[:aquatic]; kwargs...)
plot_data_metamorphs(data::AbstractDict; kwargs...) = plot_data_metamorphs(data[:metamorphs]; kwargs...)

"""
Define a default parameter set. 
This will be used as reference for the entire anaylsis. 
All parameters that are not estimated are fixed to the default values.
"""
function define_defaultparams(
    posterior_summary_larvalfit::String = datadir("sims", SAVETAG_LARVALFIT, "posterior_summary.csv"),
    posterior_summary_juvenilefit::String = datadir("sims", SAVETAG_LARVALFIT, "posterior_summary.csv")
    )
    
    p = AmphiDEB.ComponentVector(
        glb = AmphiDEB.defaultparams.glb, 
        pth = AmphiDEB.defaultparams.pth,
        spc = EcotoxSystems.ComponentVector(
            AmphiDEB.defaultparams.spc;

            # auxiliary parameters
            watercontent_larvae = 0.93, 
            watercontent_juveniles = 0.85,
            time_since_birth = 15.,
            emb_dev_time = 2.,
            H_j1_ref = 10.,
            Z_ref = Dirac(1.), # distriubtion of Z from the reference data (actual value will be assigned from file later)
            Z_corr = 1. # correction factor of the population mean of Z; simulator will calculate distribution of Z with account of Z_corr and the cv of Z_ref
        ))

        
    p.glb.t_max = 365. # setting simulation time conservatively, for cases where metamorphosis is delayed a lot
    p.glb.pathogen_inoculation_time = Inf # no pathogen inoculation
    p.glb.dX_in = [1e10, 1e10] # ad libitum feeding conditions

    p.spc.Z = truncated(Normal(1, 0.1), 0, Inf)
    # propagation of zoom factor to H_j1 is turned off => we want variability in the transition to metamorphs
    p.spc.propagate_zoom.H_j1 = 0.

    # adding point estimates as defaults

    assign_values_from_file!(
        p, 
        posterior_summary_larvalfit,
        exceptions = OrderedDict(
            "spc.Z" => (p,label,value) -> p.spc.Z_ref = Truncated(Normal(1, value), 0, 1)
            )
    )

    assign_values_from_file!(
        p, 
        posterior_summary_juvenilefit,
        exceptions = OrderedDict()
        )
    
    
    p.spc.H_j1_ref = p.spc.H_j1
    p.spc.k_M_juv = p.spc.k_M_emb
    p.spc.X_emb_int = 1. # ≈ initial dry mass of an egg (mg)

    p.spc.KD .= 0.
    p.spc.B .= 2.
    p.spc.E .= 1e10

    p.spc.emb_dev_time = estimate_emb_dev_time(p)

    return p
end

"""
Define how parameters are linked for each individual. 
FIXME: this needs to be tested further in EcotoxSystems.jl; parameter linking is for the moment done in the simulator
"""
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
Retrieve vector of exposure concentrations from data
"""
function get_exposure_vector(
    data::AbstractDict
    )

    return data[:aquatic].C_W_1 |> sort |> unique

end

"""
Define a simulator function, that takes 
"""
function define_simulator(
    data::AbstractDict;
    C_Wvec::Union{Nothing,Vector{Float64}} = nothing,
    ind_var::Bool = false
    )

    W_max_obs = maximum(data[:aquatic].wetmass_mg)
    
    function reject_early(p::ComponentVector)

        Wdry_max_obs = (1-p.spc.watercontent_larvae) * W_max_obs

        Wdry_max = calc_Wdry_max_lrv(p.spc)

        # reject parameters if the maximum structure deviates 
        # by more than a factor of 10 from the maximum observed dry mass
        # this calculation ignores the contribution of reserve to total body mass, 
        # but a plausible value should still be in the right order of magnitude...
        if !(0.1<=(Wdry_max/Wdry_max_obs)<=10)
            return true
        end

        return false

    end
    
    if isnothing(C_Wvec)
        C_Wvec = get_exposure_vector(data)
    end

    simcall(p) = AmphiDEB.ODE_simulator(p; maxiters = 1e3)

    function simulator(
        p::ComponentVector;
        return_raw = false
        )

        if reject_early(p)
            return nothing
        end

        p.spc.Z = Dirac(p.spc.Z_corr)
        p.spc.emb_dev_time = estimate_emb_dev_time(p)


        # here, Z is calculated by taking the correction factor, 
        # we can either ignore individual variability,
        # ...or pass on Z-corr to the normal distribution, assuming that cv is the same as in the reference dataset 
        # if we do this, we also need to set the number of replicates to INIT_NUM_TADPOLES below
        # for now we use the simple way
        
        if ind_var
            p.spc.Z = truncated(
                Normal(
                    p.spc.Z_corr, 
                    p.spc.Z_corr * p.spc.Z_ref.untruncated.σ
                    ), 
                    0, Inf
                )
        end

        try
                
            # simulating constant exposure throughout aquatic phase
            sim = exposure(
                p -> replicates(
                    simcall,
                    p,
                    1,#INIT_NUM_TADPOLES; 
                ), 
                p, 
                C_Wvec
            )

            # calculate total dry mass and wet mass
            sim[!,:drymass_mg] = sim.S .+ sim.E_mt 
            sim[!,:wetmass_mg] = [calc_wetmass(r, p.spc.watercontent_larvae, p.spc.watercontent_juveniles) for r in eachrow(sim)]
                                
            # convert simulation time to time since start of experiment
            sim[!,:t_exp] = sim.t .- p.spc.emb_dev_time

            if return_raw
                return sim
            end

            sim_data = OrderedDict(
                :aquatic => extract_aquatic_data(sim),
                :metamorphs => extract_metamorph_data(sim, 0.)
            )

            return sim_data
        catch e
            return nothing
        end

    end

    return simulator

end

"""
Extract aquatic data from raw simulation output. 
I.e., process ODE solution so that it is comparable with the data.
"""
function extract_aquatic_data(sim::AbstractDataFrame)::DataFrame
    return @chain sim begin
        @select(:t_exp, :larva, :wetmass_mg, :replicate, :C_W_1, :treatment_id) # 
        @subset(:larva .> 0.5) # we want larvae only
        groupby([:t_exp, :C_W_1, :treatment_id]) # for every time point and treatment
        combine(
            [:wetmass_mg,:replicate] => 
            ((m, r) -> 
                (wetmass_mg = mean(m), # calculate average dry mass
                num_tadpoles = length(unique(r))) # count how many tadpoles we have
                ) => AsTable
            )
    end
end

"""
Extract aquatic data from raw simulation output. 
I.e., process ODE solution so that it is comparable with the data.
"""
function extract_metamorph_data(sim::AbstractDataFrame, time_since_hatch)::DataFrame

    # calculate metamorphosis traits for each replicate and treatment
    metamorphs = combine(groupby(sim, [:replicate, :C_W_1, :treatment_id])) do df
        DataFrame(
            t_exp_G42 = metamorphosis_timing(df)[1] .- age_at_birth(df), 
            t_exp_G46 = metamorphosis_timing(df)[2] .- age_at_birth(df), 
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

"""
Fit model to cotrol data by estimating DEB correction factors. 
"""
function fit_control(
    species, 
    substance;
    pmcsettings = (
        :n => 10_000,
        :q_dist => 0.1,
        :t_max => 5
    ))

    data = load_data_sparlingfellers(;species=species, substance=substance, controls_only=true)
    savedir = datadir("sims", SAVETAG, species, "control")

    isdir(savedir) ? nothing : mkdir(savedir)

    defparams = define_defaultparams()
    simcontrol = define_simulator(data; C_Wvec = [0.])

    f = ModelFit(
        prior = Prior(
            # size correction factor, in relation to reference parameters (e.g. Discoglossus)
            "spc.Z_corr" => Truncated(Normal(1,2),0,Inf),
            # setting prior of maturity threshold for metamorphosis
            # the cv matches the cv of Z_corr
            "spc.H_j1" => Truncated(Normal(defparams.spc.H_j1, defparams.spc.H_j1*2), 0.1, Inf)
        ),
        defaultparams = defparams, 
        simulator = simcontrol,
        data = data, 
        response_vars = [
            [:wetmass_mg], 
            [:t_exp_G42]
            ], 
        data_weights = [
            [1.], 
            [1.]
            ], 
        grouping_vars = [
            Symbol[], 
            Symbol[]
        ], 
        time_resolved  = [true, false], 
        time_var = :t_exp, 
        plot_data = function plot_data() plot_data_controls(data) end, 
        loss_functions = EcotoxModelFitting.loss_euclidean
    )

    pmchist = run_PMC!(
        f; 
        pmcsettings...,
        savedir = savedir, 
        savetag = "controls"
    )

    p_opt = EcotoxModelFitting.bestfit(f)
    sim_opt = f.simulator(p_opt)

    plt = f.plot_data()

    @df sim_opt[:aquatic] plot!(:t_exp, :wetmass_mg, color = :steelblue, label = "Retrodicted", xlabel = "Time (d)", ylabel = "Wet mass (mg)")
    @df sim_opt[:metamorphs] vline!(:t_exp_G42, linestyle = :dash, color = :steelblue, label = "", lw = 2)

    display(plt)

    savefig(plot(plt, dpi = 300), joinpath(savedir, "controls", "controls_vpc.png"))

    eval_metamorphs = vcat(
        f.data[:metamorphs][:,[:t_exp_G42, :wetmass_G42_mg]],
        sim_opt[:metamorphs][:,[:t_exp_G42, :wetmass_G42_mg]],
    ) |> x-> @transform(x, :source = ["observed", "retrodicted/predicted"])

    
    @info "Evaluation of metamorph data"
    display(eval_metamorphs)

    CSV.write(joinpath(savedir, "controls", "eval_metamorphs.csv"), eval_metamorphs)

    return f, pmchist, ComponentVector(OrderedDict(zip(Symbol.(f.prior.labels), p_opt)))

end

"""
Plot entire data.
"""
function plot_data_exposures(data)

    plt_aqua = @df data[:aquatic] scatter(
        :t_exp, :wetmass_mg, group = :C_W_1, 
        xlabel = "Time (d)", 
        ylabel = hcat(vcat("Wet mass (mg)", repeat([""], length(unique(:C_W_1))-1))...),
        bottommargin = 5mm, leftmargin = 5mm,
        legendtitle = L"C_W", 
        title = "Effects on growth \n (fitted)"
    )

    plt_tG42 = @df data[:metamorphs] scatter(
        :C_W_1, :t_exp_G42, 
        color = :black, label = "Observed", 
        xlabel = L"\mathsf{C_W}" * " (mg/L)", 
        ylabel = "Time to G42 (d)", 
        title = "Effects on time at Gosner stage 42 \n (fitted)"
    )

    plt_WG42 = @df data[:metamorphs] scatter(
        :C_W_1, :wetmass_G42_mg, 
        color = :black, label = "Observed", 
        xlabel = L"\mathsf{C_W}" * " (mg/L)", 
        ylabel = "Wet mass \n at metamorphosis (mg)", 
        title = "Effects on mass at Gosner stage 42 \n (predicted)"
        )

    plt = plot(
        plt_aqua, plt_tG42, plt_WG42, 
        size = (1000,350), layout = (1,3), 
        titlefontsize = 10
        )

    return plt
end


"""
Configure the calibration to fit the TKTD model to exposures.
"""
function setup_modelfit_exposures(
    species, 
    substance,
    pmoa;
    p_opt_controls
    )

    pmoa_idx = findfirst(x -> x == pmoa, AmphiDEB.PMOAS) # convert pmoa from string to index
    @assert !isnothing(pmoa_idx) "Did not find PMoA $(pmoa) in PMOAS. Choose from $(PMOAS)"

    data = load_data_sparlingfellers(;species=species, substance=substance)

    defparams = define_defaultparams()
    [assign_value_by_label!(defparams, label, value) for (label,value) in zip(ComponentArrays.labels(p_opt_controls), p_opt_controls)]

    simexposures = define_simulator(data)

    median_exposure = median(@subset(data[:aquatic], :C_W_1 .> 0).C_W_1)
    lowest_exposure = minimum(@subset(data[:aquatic], :C_W_1 .> 0).C_W_1)

    f = ModelFit(
        prior = Prior(
            "spc.KD[1,$(pmoa_idx)]" => Truncated(Normal(1, 1), 0.001, 1), 
            "spc.E[1,$(pmoa_idx)]" => Truncated(Normal(median_exposure, median_exposure*10), lowest_exposure*0.1, Inf),
            "spc.B[1,$(pmoa_idx)]" => Truncated(Normal(2., 20), 0.5, Inf)
        ),
        defaultparams = defparams, 
        simulator = simexposures,
        data = data, 
        response_vars = [
            [:wetmass_mg], 
            [:t_exp_G42]
            ], 
        data_weights = [
            [1.], 
            [1.]
            ], 
        grouping_vars = [
            Symbol[:C_W_1], 
            Symbol[:C_W_1]
        ], 
        time_resolved  = [true, false], 
        time_var = :t_exp, 
        plot_data = function plot_data() plot_data_exposures(data) end, 
        loss_functions = EcotoxModelFitting.loss_euclidean
    )

    return f

end


"""
Execute TKTD fit to exposure data for given `species`, ``substance` and `pmoa`. 
"""
function fit_exposures(
    species, 
    substance, 
    pmoa; 
    p_opt_controls,
    pmcsettings = (
        :n_init => 25_000, 
        :n => 10_000,
        :q_dist => 1000/25000, 
        :t_max => 5
    ))

    savedir = datadir("sims", "SparlingFellers", species)
    savetag = "$(substance)_$(pmoa)"


    isdir(savedir) ? nothing : mkdir(savedir)

    f = setup_modelfit_exposures(
        species, substance, pmoa; 
        p_opt_controls = p_opt_controls
    )

    pmchist = run_PMC!(
        f; 
        pmcsettings...,
        savedir = savedir, 
        savetag = savetag
    )

    p_opt_exposures = EcotoxModelFitting.bestfit(f)

    ### visual check

    sim_opt = f.simulator(Vector(p_opt_exposures))

    plt_VPC = f.plot_data()

    @df sim_opt[:aquatic] plot!(:t_exp, :wetmass_mg, group = :C_W_1, subplot = 1, palette = palette(:default)[1:length(unique(:C_W_1))], label = "")
    @df sim_opt[:metamorphs] plot!(:C_W_1, :t_exp_G42, subplot = 2, color = :gray , leg = false, lw = 2)
    @df sim_opt[:metamorphs] plot!(:C_W_1, :wetmass_G42_mg, subplot = 3, color = :gray, lw = 2, leg = true, label = "G42")
    @df sim_opt[:metamorphs] plot!(:C_W_1, :wetmass_G46_mg, subplot = 3, color = :gray, lw = 2, leg = true, label = "G46", linestyle = :dot)

    display(plt_VPC)

    savefig(plot(plt_VPC, dpi = 300), joinpath(savedir, savetag, "VPC_bestfitg.png"))

    ### calculating quantitative metrics

    quant_eval_aqua = @chain sim_opt begin
        leftjoin(f.data[:aquatic], _[:aquatic], makeunique = true, on = [:t_exp, :C_W_1])
        select(:wetmass_mg_1, :wetmass_mg)
        dropmissing
        DataFrame(
            endpoint = "wetmass_mg", 
            nrmsd = nrmsd(_.wetmass_mg_1, _.wetmass_mg),
            mre = mrae(_.wetmass_mg_1, _.wetmass_mg)
        )
    end

    quant_eval_meta = @chain sim_opt begin
        leftjoin(f.data[:metamorphs], _[:metamorphs], on = [:C_W_1], makeunique = true)
        select(:t_exp_G42_1, :t_exp_G42, :wetmass_G42_mg, :wetmass_G42_mg_1)
        dropmissing
        DataFrame(
            endpoint =  ["t_exp_G42", "wetmass_G42_mg"], 
            nrmsd = [
                nrmsd(_.t_exp_G42_1, _.t_exp_G42), 
                nrmsd(_.wetmass_G42_mg_1, _.wetmass_G42_mg)
                ],
            mre = [
                mrae(_.t_exp_G42_1, _.t_exp_G42), 
                mre(_.wetmass_G42_mg_1, _.wetmass_G42_mg)
                ],
        )
    end

    quant_eval = vcat(quant_eval_aqua, quant_eval_meta)
    quant_eval[!,:species] .= species
    quant_eval[!,:substance] .= substance
    quant_eval[!,:pmoa] .= pmoa

    display(quant_eval)

    CSV.write(joinpath(savedir, savetag, "quant_eval.csv"), quant_eval)

    return f, pmchist, ComponentVector(OrderedDict(zip(Symbol.(f.prior.labels), p_opt_exposures))), plt_VPC, quant_eval
end

"""
Maps parameter labels to LaTeX strings for plotting and tables. 
"""
paramlabels = OrderedDict(
    "spc.Z_corr" => L"\overline{Z}_{corr}", 
    "spc.eta_AS_emb" => L"\eta_{AS}^{emb}", 
    
    "spc.log_k_D_G" => L"ln k_{D,G}",
    "spc.KD[1,1]" => L"k_{D,G}",
    "spc.log_e_G" => L"ln e_G",
    "spc.B[1,1]" => L"b_G",
    
    "spc.log_k_D_M" => L"ln k_{D,M}",
    "spc.KD[1,2]" => L"k_{D,M}",
    "spc.log_e_M" => L"ln e_M",
    "spc.B[1,2]" => L"b_M",
    
    "spc.log_k_D_A" => L"ln k_{D,A}",
    "spc.KD[1,3]" => L"k_{D,A}",
    "spc.log_e_A" => L"ln(e_A)",
    "spc.B[1,3]" => L"b_A",

    "spc.log_k_D_R" => L"ln k_{D,R}",
    "spc.KD[1,4]" => L"\k_{D,R}",
    "spc.log_e_R" => L"ln(e_R)",
    "spc.B[1,4]" => L"b_R",

    "spc.log_k_D_Hneg" => L"ln k_{D,H^-}",
    "spc.KD[1,5]" => L"\k_{D,H^-}",
    "spc.log_e_KAP" => L"ln(e_{H^-})",
    "spc.B[1,5]" => L"b_{H^-}",

    "spc.log_k_D_Hpos" => L"ln k_{D,H^+}",
    "spc.KD[1,6]" => L"\k_{D,H^+}",
    "spc.log_e_KAP" => L"ln(e_{H^+})",
    "spc.B[1,6]" => L"b_{H^+}",

    "spc.log_k_D_KAP" => L"ln k_{D,\kappa}",
    "spc.KD[1,7]" => L"\k_{D,\kappa}",
    "spc.log_e_KAP" => L"ln(e_{\kappa})",
    "spc.B[1,7]" => L"b_{\kappa}",

    "spc.E[1,1]" => L"e_G",
    "spc.E[1,2]" => L"e_M",
    "spc.E[1,3]" => L"e_A",
    "spc.E[1,5]" => L"e_{H^{-}}",
    "spc.E[1,6]" => L"e_{H^{+}}",
    "spc.E[1,7]" => L"e_{\kappa^{-}}",
)