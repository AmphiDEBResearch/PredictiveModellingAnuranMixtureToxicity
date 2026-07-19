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
p_opt_A = p_opt_from_tag(joinpath("Discoglossus_24D_2025-06-23_numtadpoles", "Discoglossus_24D_A"))


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

# ======================================================== #
# Second case: fG, fM, fA (with their simulations)
# ======================================================== #

using Random, Statistics, DataFrames
include("bootstrap_metrics.jl")

# ------------------ Wet mass ------------------ #

δG = get_metrics(fG, sim_opt_G) |> x -> @transform(x, :pmoa = "G")
δM = get_metrics(fM, sim_opt_M) |> x -> @transform(x, :pmoa = "M")
δA = get_metrics(fA, sim_opt_A) |> x -> @transform(x, :pmoa = "A")

MAPE_CI_G = get_MAPE_CI(fG, sim_opt_G; metric = :MAPE_wetmass, num_iters = 2000) |>
    x -> @transform(x, :pmoa = "G")
MAPE_CI_M = get_MAPE_CI(fM, sim_opt_M; metric = :MAPE_wetmass, num_iters = 2000) |>
    x -> @transform(x, :pmoa = "M")
MAPE_CI_A = get_MAPE_CI(fA, sim_opt_A; metric = :MAPE_wetmass, num_iters = 2000) |>
    x -> @transform(x, :pmoa = "A")

MAPE_CI = vcat(MAPE_CI_G, MAPE_CI_M, MAPE_CI_A)

NSE_CI_G = get_metric_CI(fG, sim_opt_G; metric = :NSE_wetmass, point_col = :NSE, num_iters = 2000) |>
    x -> @transform(x, :pmoa = "G")
NSE_CI_M = get_metric_CI(fM, sim_opt_M; metric = :NSE_wetmass, point_col = :NSE, num_iters = 2000) |>
    x -> @transform(x, :pmoa = "M")
NSE_CI_A = get_metric_CI(fA, sim_opt_A; metric = :NSE_wetmass, point_col = :NSE, num_iters = 2000) |>
    x -> @transform(x, :pmoa = "A")

NSE_CI = vcat(NSE_CI_G, NSE_CI_M, NSE_CI_A)

pMAPE = @df sort(MAPE_CI, :pmoa) groupedbar(
    string.(:C_W_1), :MAPE,
    yerr = :MAPE_upper .- :MAPE_lower,
    group = :pmoa, bar_position = :dodge,
    legendtitle = "PMoA",
    fillalpha = 0.5,
    xlabel = "Treatment level (mg/L)", ylabel = "MAPE"
)

pNSE = @df NSE_CI groupedbar(
    string.(:C_W_1), :NSE,
    yerr = :NSE_upper .- :NSE_lower,
    group = :pmoa, bar_position = :dodge,
    fillalpha = 0.5,
    xlabel = "Treatment level (mg/L)", ylabel = "NSE",
    leg = false
)

pMET = plot(
    pMAPE, pNSE, layout = (2, 1),
    title = ["Larval wet mass" ""],
    bottommargin = 5mm, topmargin = 5mm, leftmargin = 5mm,
    size = (500, 750)
)

# ------------------ Fraction tadpoles ------------------ #

MAPE_CI_G_ft = get_MAPE_CI(fG, sim_opt_G;
    metric = :MAPE_fracttadpoles, num_iters = 2000
) |> x -> @transform(x, :pmoa = "G")

MAPE_CI_M_ft = get_MAPE_CI(fM, sim_opt_M;
    metric = :MAPE_fracttadpoles, num_iters = 2000
) |> x -> @transform(x, :pmoa = "M")

MAPE_CI_A_ft = get_MAPE_CI(fA, sim_opt_A;
    metric = :MAPE_fracttadpoles, num_iters = 2000
) |> x -> @transform(x, :pmoa = "A")

MAPE_CI_ft = vcat(MAPE_CI_G_ft, MAPE_CI_M_ft, MAPE_CI_A_ft)

NSE_CI_G_ft = get_metric_CI(fG, sim_opt_G;
    metric = :NSE_fracttadpoles, point_col = :NSE, num_iters = 2000
) |> x -> @transform(x, :pmoa = "G")

NSE_CI_M_ft = get_metric_CI(fM, sim_opt_M;
    metric = :NSE_fracttadpoles, point_col = :NSE, num_iters = 2000
) |> x -> @transform(x, :pmoa = "M")

NSE_CI_A_ft = get_metric_CI(fA, sim_opt_A;
    metric = :NSE_fracttadpoles, point_col = :NSE, num_iters = 2000
) |> x -> @transform(x, :pmoa = "A")

NSE_CI_ft = vcat(NSE_CI_G_ft, NSE_CI_M_ft, NSE_CI_A_ft)

pMAPE_ft = @df sort(MAPE_CI_ft, :pmoa) groupedbar(
    string.(:C_W_1), :MAPE,
    yerr = :MAPE_upper .- :MAPE_lower,
    group = :pmoa, bar_position = :dodge,
    legendtitle = "PMoA",
    fillalpha = 0.5,
    xlabel = "Treatment level (mg/L)", ylabel = "MAPE"
)

pNSE_ft = @df sort(NSE_CI_ft, :pmoa) groupedbar(
    string.(:C_W_1), :NSE,
    yerr = :NSE_upper .- :NSE_lower,
    group = :pmoa, bar_position = :dodge,
    fillalpha = 0.5,
    xlabel = "Treatment level (mg/L)", ylabel = "NSE",
    leg = false
)

pMET_ft = plot(
    pMAPE_ft, pNSE_ft, layout = (2, 1),
    title = ["Fraction of tadpoles" ""],
    bottommargin = 5mm, topmargin = 5mm, leftmargin = 5mm,
    size = (500, 750)
)

# ------------------ Save/combined ------------------ #

savefig(plot(pMET, pMET_ft, layout = (1, 2), size = (800, 800)),
        plotsdir("24d_metrics_larvae.png"))

