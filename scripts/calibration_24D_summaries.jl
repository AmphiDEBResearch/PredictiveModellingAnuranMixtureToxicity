include("boilerplate.jl")

using Revise

const SAVETAG_LARVALFIT = "input/Discoglossus_larvae" # directory from which larval/metamorph parameters are loaded
const SAVETAG_JUVENILEFIT = "input/Discoglossus_juveniles" # directory from which juvenile/adult parameters are loaded
const SAVETAG = "Discoglossus_24D"
using Revise

includet(scriptsdir("Discoglossus_galganoi_24D", "fit.jl")) 

"""
    p_opt_from_tag(savetag::AbstractString)

Read best-fitting parameter combination from a posterior_summary.csv file.
"""
function p_opt_from_tag(savetag::AbstractString)

    posterior_summary = CSV.read(datadir("sims", savetag, "posterior_summary.csv"), DataFrame)
    return posterior_summary.best_fit

end

# read best fits for alternative PMoA
p_opt_G = p_opt_from_tag(joinpath("Discoglossus_24D_2025-06-23_numtadpoles", "Discoglossus_24D_G"))
p_opt_M = p_opt_from_tag(joinpath("Discoglossus_24D_2025-06-23_numtadpoles", "Discoglossus_24D_M"))
p_opt_A = p_opt_from_tag(joinpath("Discoglossus_24D_2025-06-23_numtadpoles", "Discoglossus_24D_A"));

# run simulations for alternative PMoAs

n = 100

fG = setup_modelfit("G")
sim_opt_G = [fG.simulator(p_opt_G) for _ in 1:n]

fM = setup_modelfit("M")
sim_opt_M = [fM.simulator(p_opt_M) for _ in 1:n]

fA = setup_modelfit("A")
sim_opt_A = [fA.simulator(p_opt_A) for _ in 1:n];

f = fG; plt = plot_data()

plot_sims!(plt, sim_opt_G, label = "Best fit (G)", color = :chocolate2)
plot_sims!(plt, sim_opt_M, label = "Best fit (M)", color = :steelblue)
plot_sims!(plt, sim_opt_A, label = "Best fit (A)", color = :magenta)
#plot!(subplot = 1, xlim = (0,20))
#plot!(subplot = 2, xlim = (0,15))

savefig(plot(plt, dpi = 400), plotsdir("Discoglossus_24D_summaries_PMoA_comparison.png"))

plt

"""
    accepted_from_tag(savetag::AbstractString)

Read accepted samples from a file.
"""
function accepted_from_tag(savetag::AbstractString)

    accepted = CSV.read(datadir("sims", savetag, "accepted.csv"), DataFrame)
    
    return accepted

end

accepted_G = accepted_from_tag(joinpath("Discoglossus_24D_2025-06-23_numtadpoles", "Discoglossus_24D_G")) 
accepted_M = accepted_from_tag(joinpath("Discoglossus_24D_2025-06-23_numtadpoles", "Discoglossus_24D_M"))
accepted_A = accepted_from_tag(joinpath("Discoglossus_24D_2025-06-23_numtadpoles", "Discoglossus_24D_A"));

# calculate posterior probabilities based on rejection approach

accepted_G[!,:model] .= "G"
accepted_M[!,:model] .= "M"
accepted_A[!,:model] .= "A"

losses = vcat(
    accepted_G[:,[:model,:loss, :weight]],
    accepted_M[:,[:model,:loss, :weight]],
    accepted_A[:,[:model,:loss, :weight]]
)

ϵ = quantile(losses.loss, 0.25)

losses_accepted = @subset(losses, :loss .<= ϵ)

probs = countmap(losses_accepted.model) |> OrderedDict |>
x -> DataFrame(model = x.keys, prob = x.vals ./ sum(values(x)))

for model in ["G", "M", "A"] # add 0s for entirely rejected models
    if !(model in probs.model)
        append!(probs, DataFrame(
            model = model, 
            prob = 0.
        ))
    end
end

probs

# calculate errors of posterior probabilities using bootstrapping

bootstrap_samplesize = Int(0.5*nrow(losses))
n_samples = 1000

bootstrap_probs_G = []
bootstrap_probs_A = []
bootstrap_probs_M = []

for i in 1:n_samples

    bootstrap_idxs = sample(1:nrow(losses), bootstrap_samplesize)

    losses_i = losses[bootstrap_idxs,:] # get re-sampled subset of losses
    epsilon = quantile(losses_i.loss, 0.25) # compute the threshold for the bootstrap sample

    acc = @subset(losses_i, :loss .<= epsilon, :weight .> 0)

    probs_i = countmap(acc.model) |> OrderedDict |>
    x -> DataFrame(model = x.keys, prob = x.vals ./ sum(values(x)))

    # add 0s for entirely rejected models
    for model in ["G", "M", "A"] 
        if !(model in probs_i.model)
            append!(probs_i, DataFrame(
                model = model, 
                prob = 0.
            ))
        end
    end

    sort!(probs_i, :model)

    push!(bootstrap_probs_A, probs_i.prob[1])
    push!(bootstrap_probs_G, probs_i.prob[2])
    push!(bootstrap_probs_M, probs_i.prob[3])

end

ci(x) = quantile(x, 0.75) - quantile(x, 0.25)

plt_probs = @df probs bar(
    :model, :prob, 
    permute = (:y, :x),
    yerrors = [ci(bootstrap_probs_A), ci(bootstrap_probs_M), ci(bootstrap_probs_G)], 
    xlim = (0,3), ylim = (0,1), 
    xlabel = "PMoA", ylabel = "Posterior probability",
    bar_width = 0.5, leg = false, 
    fillcolor = :gray, fillalpha = .5,
    title = "PMoA selection"
    )
display(plt_probs)
savefig(plot(plt_probs, dpi = 300), plotsdir("Discoglossus_24D_summaries_PMoA_probabilities.png"))

vcat(
    @subset(accepted_G, :loss .== minimum(:loss))[:,[:model,:loss]],
    @subset(accepted_M, :loss .== minimum(:loss))[:,[:model,:loss]],
    @subset(accepted_A, :loss .== minimum(:loss))[:,[:model,:loss]]
) |> x-> @transform(x, :loss = round.(:loss, sigdigits = 2))
