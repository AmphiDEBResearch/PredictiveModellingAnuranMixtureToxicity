include("boilerplate.jl")

using Revise

const SAVETAG_LARVALFIT = "input/Discoglossus_larvae" # directory from which larval/metamorph parameters are loaded
const SAVETAG_JUVENILEFIT = "input/Discoglossus_juveniles" # directory from which juvenile/adult parameters are loaded
const SAVETAG_TKTD = "Discoglossus_Flupyradifurone"
const SAVETAG = "Discoglossus_Flupyradifurone"

using Revise

includet(scriptsdir("Discoglossus_galganoi_Flupyradifurone", "fit.jl"))
includet(scriptsdir("Discoglossus_galganoi_Flupyradifurone", "cross_validation.jl")) 

savedir = joinpath("Discoglossus_Flupyradifurone_2025-06-23_numtadpoles")

# ======================================== #
# Simulations for alternative PMoAs
# ======================================== #


@suppress begin
    global fG = setup_modelfit("G")
    global p_opt_G = CSV.read(datadir("sims", savedir, "$(SAVETAG_TKTD)_G", "posterior_summary.csv"), DataFrame).best_fit
    global sim_opt_G = [fG.simulator(p_opt_G) for _ in 1:100]

    global fKAP = setup_modelfit("KAP")
    global p_opt_KAP = CSV.read(datadir("sims", savedir, "$(SAVETAG_TKTD)_KAP", "posterior_summary.csv"), DataFrame).best_fit
    global sim_opt_KAP = [fKAP.simulator(p_opt_KAP) for _ in 1:100];
end;

# ======================================== #
# Generate plot
# ======================================== #

f = fG;
plt = plot_data()
plot_sims!(plt, sim_opt_G, label = "Best fit (G)", color = :chocolate2)
#plot_sims!(plt, sim_opt_M, label = "Best fit (M)", color = :steelblue)
#plot_sims!(plt, sim_opt_A, label = "Best fit (A)", color = :magenta)
plot_sims!(plt, sim_opt_KAP, label = "Best fit (κ)", color = :mediumseagreen)
#plot!(subplot = 1, xlim = (0,20))
#plot!(subplot = 2, xlim = (0,15))

savefig(plot(plt, dpi = 400), plotsdir("Discoglossus_Flupyradifurone_summaries_PMoA_comparison.png"))

plt

