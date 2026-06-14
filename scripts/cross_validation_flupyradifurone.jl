include("boilerplate.jl")
using Revise
using ProgressMeter

const SAVETAG_LARVALFIT = "input/Discoglossus_larvae" # directory from which larval/metamorph parameters are loaded
const SAVETAG_JUVENILEFIT = "input/Discoglossus_juveniles" # directory from which juvenile/adult parameters are loaded
const SAVEDIR = "Discoglossus_Flupyradifurone_2025-06-23_numtadpoles" # directory from which TKTD parameters are loaded
const SAVETAG = SAVETAG_TKTD = "Discoglossus_Flupyradifurone" 

using Revise

includet(scriptsdir("Discoglossus_galganoi_Flupyradifurone", "fit.jl")) 
includet(scriptsdir("Discoglossus_galganoi_Flupyradifurone", "cross_validation.jl"))

# ======================================== #
# Simulations per PMoA
# ======================================== #

# ---- PMoA G

pmoa_idx = 1
pmoa = PMOAS[pmoa_idx]
f = setup_modelfit(pmoa);

p_opt = CSV.read(datadir("sims", SAVEDIR, "$(SAVETAG_TKTD)_$(pmoa)", "posterior_summary.csv"), DataFrame).best_fit
sim_opt_G = @showprogress [f.simulator(p_opt) for _ in 1:100];

quant_eval_G = quant_eval_metamorphs(f, sim_opt_G)


# ---- PMoA κ

pmoa_idx = 7
pmoa = PMOAS[pmoa_idx]
f = setup_modelfit(pmoa);

p_opt = CSV.read(datadir("sims", SAVEDIR, "$(SAVETAG_TKTD)_$(pmoa)", "posterior_summary.csv"), DataFrame).best_fit
sim_opt_KAP = @showprogress [f.simulator(p_opt) for _ in 1:1_00];

quant_eval_KAP = quant_eval_metamorphs(f, sim_opt_KAP)

# ======================================== #
# Generate plot
# ======================================== #
plt = plot_metamorphs(
    leftmargin = 10mm, 
    bottommargin = 10mm,
    xlabel = "Flupyradifurone (mg/L)"
    )

sim = EcotoxModelFitting.extract_simkey(sim_opt_G, :metamorphs) |> 
x -> @transform(x, :treatment_id = denserank(:C_W_1))
sim_retro = @subset(sim, :treatment_id .== 1)
sim_pred = @subset(sim, :treatment_id .> 1)

@df sim_retro violin!(
    plt, subplot = 1,
    string.(:treatment_id), :t_exp_G46, 
    side = :right,
    fillalpha = .25,
    color = :gray,
    fillstyle = ://,
    label = "Retrodicted",
    )

@df sim_pred violin!(
    plt, subplot = 1,
    string.(:treatment_id), :t_exp_G46, 
    side = :right,
    fillalpha = .25,
    color = :chocolate2,
    label = "Predicted (G)",
    title = "Timing of G46 \n MAPE = $(round(quant_eval_G.mape[1], sigdigits = 2))%"
    )

@df sim_retro violin!(
    plt, subplot = 2,
    string.(:treatment_id), :wetmass_G46_mg, 
    side = :right,
    fillalpha = .25,
    color = :chocolate2,
    fillstyle = ://,
    label = "Predicted (G)",
    )

@df sim_pred violin!(
    plt, subplot = 2,
    string.(:treatment_id), :wetmass_G46_mg, 
    side = :right,
    fillalpha = .25,
    color = :chocolate2,
    label = "Retrodicted",
    title = "Mass at G46 \n MAPE = $(round(quant_eval_G.mape[2], sigdigits = 2))%"
    )


sim = EcotoxModelFitting.extract_simkey(sim_opt_KAP, :metamorphs) |> 
x -> @transform(x, :treatment_id = denserank(:C_W_1))
sim_retro = @subset(sim, :treatment_id .== 1)
sim_pred = @subset(sim, :treatment_id .> 1)

@df sim_pred violin!(
    plt, subplot = 1,
    string.(:treatment_id), :t_exp_G46, 
    side = :right,
    fillalpha = .25,
    color = :mediumseagreen,
    label = "Predicted (κ)",
    title = "Timing of G46 \n MAPE = $(round(quant_eval_G.mape[1], sigdigits = 2))%"
    )

@df sim_pred violin!(
    plt, subplot = 2,
    string.(:treatment_id), :wetmass_G46_mg, 
    side = :right,
    fillalpha = .25,
    color = :mediumseagreen,
    label = "Retrodicted",
    title = "Mass at G46 \n MAPE = $(round(quant_eval_G.mape[2], sigdigits = 2))%"
    )

savefig(
    plot(plt, dpi = 400), 
    plotsdir("CrossValidation_Discoglossus_Flupyradifurone_PMoA_comparison.png")
    #datadir("sims", "$(SAVETAG)_$(PMOAS[pmoa_idx])", "CrossValidation_metamorphs.png")
    )

display(plt)
