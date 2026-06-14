include("boilerplate.jl")

using Revise

const SAVETAG_LARVALFIT = "input/Discoglossus_larvae" # directory from which larval/metamorph parameters are loaded
const SAVETAG_JUVENILEFIT = "input/Discoglossus_juveniles" # directory from which juvenile/adult parameters are loaded
const SAVETAG = "Discoglossus_24D"
using Revise

includet(scriptsdir("Discoglossus_galganoi_24D", "fit.jl")) 

# ======================================== #
# Read fitted parameters
# ======================================== #

"""
    p_opt_from_tag(savetag::AbstractString)

Read best-fitting parameter combination from a posterior_summary.csv file.
"""
function p_opt_from_tag(savetag::AbstractString)

    posterior_summary = CSV.read(datadir("sims", savetag, "posterior_summary.csv"), DataFrame)
    return posterior_summary.best_fit

end

p_opt_G = p_opt_from_tag(joinpath("Discoglossus_24D_2025-06-23_numtadpoles", "Discoglossus_24D_G"))
p_opt_M = p_opt_from_tag(joinpath("Discoglossus_24D_2025-06-23_numtadpoles", "Discoglossus_24D_M"))
p_opt_A = p_opt_from_tag(joinpath("Discoglossus_24D_2025-06-23_numtadpoles", "Discoglossus_24D_A"));


# ======================================== #
# Run simulations for alternative PMoAs
# ======================================== #

n = 100

fG = setup_modelfit("G")
sim_opt_G = [fG.simulator(p_opt_G) for _ in 1:n]

fM = setup_modelfit("M")
sim_opt_M = [fM.simulator(p_opt_M) for _ in 1:n]

fA = setup_modelfit("A")
sim_opt_A = [fA.simulator(p_opt_A) for _ in 1:n];

# ======================================== #
# Generate plot
# ======================================== #

f = fG; plt = plot_data()

plot_sims!(plt, sim_opt_G, label = "Best fit (G)", color = :chocolate2)
plot_sims!(plt, sim_opt_M, label = "Best fit (M)", color = :steelblue)
plot_sims!(plt, sim_opt_A, label = "Best fit (A)", color = :magenta)
#plot!(subplot = 1, xlim = (0,20))
#plot!(subplot = 2, xlim = (0,15))

plt

savefig(plot(plt, dpi = 400), plotsdir("Discoglossus_24D_summaries_PMoA_comparison.png"))

plt