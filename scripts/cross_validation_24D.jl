include("boilerplate.jl")

const SAVETAG_LARVALFIT = "input/Discoglossus_larvae" # directory from which larval/metamorph parameters are loaded
const SAVETAG_JUVENILEFIT = "input/Discoglossus_juveniles" # directory from which juvenile/adult parameters are loaded
const SAVEDIR = "Discoglossus_24D_2025-06-23_numtadpoles" # directory from which TKTD parameters are loaded
const SAVETAG_TKTD = "Discoglossus_24D" 
const SAVETAG = "Discoglossus_24D"

using Revise

includet(scriptsdir("Discoglossus_galganoi_24D", "fit.jl")) 
includet(scriptsdir("Discoglossus_galganoi_24D", "cross_validation.jl"))

# ======================================== #
# Predictions for PMoA M
# ======================================== #

using ProgressMeter 

pmoa_idx = 2
pmoa = PMOAS[pmoa_idx]
f = setup_modelfit(pmoa); # reconstructing ModelFit instance

p_opt = CSV.read(datadir("sims", SAVEDIR, "$(SAVETAG_TKTD)_$(pmoa)", "posterior_summary.csv"), DataFrame).best_fit
sim_opt_M = @showprogress [f.simulator(p_opt) for _ in 1:100];

quant_eval_M = quant_eval_metamorphs(f, sim_opt_M)

# ======================================== #
# Predictions for PMoA A
# ======================================== #

pmoa_idx = 3
pmoa = PMOAS[pmoa_idx]
f = setup_modelfit(pmoa);

p_opt = CSV.read(datadir("sims", SAVEDIR, "$(SAVETAG_TKTD)_$(pmoa)", "posterior_summary.csv"), DataFrame).best_fit
sim_opt_A = [f.simulator(p_opt) for _ in 1:100] 


# ======================================== #
# Plot data + all predictions 
# ======================================== #

plt = plot_metamorphs(
    bottommargin = 10mm, leftmargin = 10mm
)
plot!(subplot = 1, leg = :topleft)
plot!(subplot = 1, ylim = (10, 80)) # adjusting plot manually to include sample sizes
plot!(subplot = 2, ylim = (50, 350))

# ---- sample size annotatons for timing of GS 46
ypos = 55
getn(treatment_id) =   nrow(@subset(f.data[:metamorphs], :treatment_id .== treatment_id)) 
[annotate!(
    i-0.75, ypos, Plots.text("n=$(getn(i))", 10), 
    subplot = 1
    ) for i in sort(unique(f.data[:metamorphs].treatment_id))]

# ---- sample size annotations for mass at GS 46
ypos = 300
getn(treatment_id) =   nrow(@subset(f.data[:metamorphs], :treatment_id .== treatment_id)) 
[annotate!(
    i-0.75, ypos, Plots.text("n=$(getn(i))", 10), 
    subplot = 2
    ) for i in sort(unique(f.data[:metamorphs].treatment_id))]
plt

# ---- predictions for M

sim = EcotoxModelFitting.extract_simkey(sim_opt_M, :metamorphs)
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
    title = "Timing of G46 \n MAPE = $(round(quant_eval_M.mape[1], sigdigits = 2))%"
    )

@df sim_pred violin!(
    plt, subplot = 1,
    string.(:treatment_id), :t_exp_G46, 
    side = :right,
    fillalpha = .25,
    color = :steelblue,
    label = "Predicted (M)"
    )

@df sim_retro violin!(
    plt, subplot = 2,
    string.(:treatment_id), :wetmass_G46_mg, 
    side = :right,
    fillalpha = .25,
    color = :gray,
    fillstyle = ://,
    label = "Retrodicted",
    title = "Timing of G46 \n MRE = $(round(quant_eval_M.mape[1], sigdigits = 2))%"
    )

@df sim_pred violin!(
    plt, subplot = 2,
    string.(:treatment_id), :wetmass_G46_mg, 
    side = :right,
    fillalpha = .25,
    color = :steelblue,
    label = "Predicted (M)",
    title = "Timing of G46 \n MAPE = $(round(quant_eval_M.mape[2], sigdigits = 2))%"
    )
  

# ---- predictions for A

sim = EcotoxModelFitting.extract_simkey(sim_opt_A, :metamorphs)
sim_retro = @subset(sim, :treatment_id .== 1)
sim_pred = @subset(sim, :treatment_id .> 1)

@df sim_pred violin!(
    plt, subplot = 1,
    string.(:treatment_id), :t_exp_G46, 
    side = :right,
    fillalpha = .25,
    color = :magenta,
    label = "Predicted (A)",
    )

@df sim_pred violin!(
    plt, subplot = 2,
    string.(:treatment_id), :wetmass_G46_mg, 
    side = :right,
    fillalpha = .25,
    color = :magenta,
    label = "Predicted (A)",
    )    

plot!(plt, size = (850,400))

savefig(
    plot(plt, dpi = 400), 
    plotsdir("CrossValidation_Discoglossus_24D_PMoA_comparison.png")
    )

display(plt)
