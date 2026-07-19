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
plot_sims!(plt, sim_opt_KAP, label = "Best fit (κ)", color = :mediumseagreen)

plot!(plt, ylim = (80,350), subplot = 3)
plt

savefig(plot(plt, dpi = 400), plotsdir("Discoglossus_Flupyradifurone_summaries_PMoA_comparison.png"))


# ======================================================== #
# Calculate quantitative metrics per PMoA
# ======================================================== #

loss_G = generate_loss_function(fG)
lossvals_G = [loss_G(fG.data, sim) for sim in sim_opt_G]
loss_KAP = generate_loss_function(fKAP)
lossvals_KAP = [loss_KAP(fKAP.data, sim) for sim in sim_opt_KAP]

histogram(
    lossvals_G, normalize = :pdf, 
    fillalpha = .25, label = "G", 
    title = "Distances of accepted particles per PMoA \n (larval calibration data)", 
    titlefontsize = 12, 
    xlabel = "distance", 
    ylabel = "density"
    ) 
histogram!(
    lossvals_KAP, 
    normalize = :pdf, fillalpha = .25, label = "κ"
    )


# ======================================================== #
# Calculation of MAPE and NSE with bootstrapping
# ======================================================== #

MAPE(obs,sim) = mean(@. abs(100 * (sim - obs)/obs))
NSE(obs, sim) = sum((obs .- sim) .^2) ./ sum(obs .- mean(obs).^2)

using Random, Statistics, DataFrames
includet("bootstrap_metrics.jl")

δG = get_metrics(fG, sim_opt_G) |> x->@transform(x, :pmoa = "G")
δκ = get_metrics(fKAP, sim_opt_KAP) |> x->@transform(x, :pmoa = "KAP")
MAPE_CI_G = get_MAPE_CI(fG, sim_opt_G;   metric = :MAPE_wetmass, num_iters = 2000)|> 
x->@transform(x, :pmoa = "G")
MAPE_CI_KAP = get_MAPE_CI(fKAP, sim_opt_KAP;   metric = :MAPE_wetmass, num_iters = 2000) |> 
x->@transform(x, :pmoa = "KAP")

MAPE_CI = vcat(MAPE_CI_G, MAPE_CI_KAP) 

NSE_CI_G = get_metric_CI(
    fG, sim_opt_G;
    metric = :NSE_wetmass,
    point_col = :NSE,
    num_iters = 2000
) |> x -> @transform(x, :pmoa = "G")

NSE_CI_KAP = get_metric_CI(
    fKAP, sim_opt_KAP;
    metric = :NSE_wetmass,
    point_col = :NSE,
    num_iters = 2000
) |> x -> @transform(x, :pmoa = "KAP")

NSE_CI = vcat(NSE_CI_G, NSE_CI_KAP)

# ---- construct bar plot with NSE/MAPE + CI

pMAPE = @df MAPE_CI groupedbar(
    string.(:C_W_1), :MAPE, 
    yerr = :MAPE_upper .- :MAPE_lower,
    group = :pmoa, bar_position = :dodge, 
    label = ["G" "κ"], legendtitle = "PMoA",
    fillalpha = .5, 
    xlabel = "Treatment level (mg/L)", ylabel = "MAPE"
)

pNSE = @df NSE_CI groupedbar(
    string.(:C_W_1), :NSE,
    yerr = :NSE_upper .- :NSE_lower,
    group = :pmoa, bar_position = :dodge,
    fillalpha = .5, xlabel = "Treatment level (mg/L)", ylabel = "NSE", 
    ylim = (-0.05, 0), leg = false
)
pMET = plot(
    pMAPE, pNSE, layout = (2,1), 
    title = ["Larval wet mass" ""], 
    bottommargin=5mm, topmargin=5mm, leftmargin=5mm,
    size = (500,750)
)

savefig(
    plot(
        pMET, background=:transparent, dpi=300), 
        plotsdir("flupy_metrics_wetmass_larvae.png")
        )

# now the same for fraction of tadpoles


metrics = combine(groupby(joined, [:C_W_1, :num_sim])) do df
    return DataFrame(
        MAPE_wetmass = MAPE(df.wetmass_mg_obs, df.wetmass_mg_sim), 
        MAPE_fracttadpoles =  MAPE(df.fract_tadpoles_obs, df.fract_tadpoles_sim), 
        NSE_wetmass = NSE(df.wetmass_mg_obs, df.wetmass_mg_sim), 
        NSE_fracttadpoles = NSE(df.fract_tadpoles_obs, df.fract_tadpoles_sim), 
        )
end

MAPE_CI_G_ft = get_MAPE_CI(fG, sim_opt_G;
    metric = :MAPE_fracttadpoles, num_iters = 2000
) |> x -> @transform(x, :pmoa = "G")

MAPE_CI_KAP_ft = get_MAPE_CI(fKAP, sim_opt_KAP;
    metric = :MAPE_fracttadpoles, num_iters = 2000
) |> x -> @transform(x, :pmoa = "KAP")

MAPE_CI_ft = vcat(MAPE_CI_G_ft, MAPE_CI_KAP_ft)

NSE_CI_G_ft = get_metric_CI(
    fG, sim_opt_G;
    metric = :NSE_fracttadpoles,
    point_col = :NSE,
    num_iters = 2000
) |> x -> @transform(x, :pmoa = "G")

NSE_CI_KAP_ft = get_metric_CI(
    fKAP, sim_opt_KAP;
    metric = :NSE_fracttadpoles,
    point_col = :NSE,
    num_iters = 2000
) |> x -> @transform(x, :pmoa = "KAP")

NSE_CI_ft = vcat(NSE_CI_G_ft, NSE_CI_KAP_ft)

pMAPE_ft = @df MAPE_CI_ft groupedbar(
    string.(:C_W_1), :MAPE,
    yerr = :MAPE_upper .- :MAPE_lower,
    group = :pmoa, bar_position = :dodge,
    label = ["G" "κ"], legendtitle = "PMoA",
    fillalpha = .5,
    xlabel = "Treatment level (mg/L)", ylabel = "MAPE"
)

pNSE_ft = @df NSE_CI_ft groupedbar(
    string.(:C_W_1), :NSE,
    yerr = :NSE_upper .- :NSE_lower,
    group = :pmoa, bar_position = :dodge,
    fillalpha = .5, xlabel = "Treatment level (mg/L)", ylabel = "NSE",
    #ylim = (-0.05, 0),
    leg = false
)

pMET_ft = plot(
    pMAPE_ft, pNSE_ft, layout = (2,1),
    title = ["Larval fraction tadpoles" ""],
    bottommargin=5mm, topmargin=5mm, leftmargin=5mm,
    size = (500,750)
)

savefig(plot(pMET_ft, background=:transparent, dpi=300), plotsdir("flupy_metrics_fraction_tadpoles_larvae.png"))

