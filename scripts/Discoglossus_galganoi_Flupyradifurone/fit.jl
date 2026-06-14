include("data.jl")
include("parameters.jl")
include("simulation.jl")

function setup_modelfit(
    pmoa::AbstractString; 
    loss_function::Function = EcotoxModelFitting.loss_euclidean
    )

    pmoa_idx = findfirst(x -> x == pmoa, PMOAS) # convert pmoa from string to index
    @assert !isnothing(pmoa_idx) "Did not find PMoA $(pmoa) in PMOAS"

    if !isdir(datadir("sims", SAVETAG))
        mkdir(datadir("sims", SAVETAG))
    end
    
    # set up logger
    global io = open(datadir("sims", SAVETAG, "log.txt"), "w+")
    global logger = SimpleLogger(io)

    data = load_data()

    f = ModelFit( 
        prior = Prior(
            "spc.KD[1,$(pmoa_idx)]" => Truncated(Normal(0.5, 0.5), 0.001, 1), 
            "spc.E[1,$(pmoa_idx)]" => Truncated(Normal(50, 50), 1, Inf),
            "spc.B[1,$(pmoa_idx)]" => Truncated(Normal(2., 4), 0.5, Inf)
            ),
        defaultparams = define_defaultparams(), 
        simulator = simulator, 
        data = data, 
        response_vars = [
            [:wetmass_mg, :fract_tadpoles], 
            Symbol[] # leaving this empty == metamorphs will be ignored during calibration
            #[:t_exp_G42, :t_exp_G46, :wetmass_G42_mg, :wetmass_G46_mg]
        ],
        data_weights = [
            [1., 1.],
            Float64[]
        ],
        grouping_vars = [[:C_W_1], [:C_W_1]], 
        time_resolved = [true, false],
        time_var = :t_exp,
        plot_data = plot_data,
        loss_functions = loss_function
    )

    return f

end


function fit_model!(
    f; 
    n = 100, 
    n_init = 100, 
    q_dist = .1, 
    t_max = 3, 
    savetag = nothing, 
    paramlabels = paramlabels, 
    continue_from = nothing,
    n_posterior_check = 100,
    kwargs...
    )

    if !isdir(plotsdir(savetag))
        mkdir(plotsdir(savetag))
    end

    pmchist = run_PMC!(
        f; 
        n = n, 
        n_init = n_init, 
        q_dist = q_dist, 
        t_max = t_max, 
        savetag = savetag, 
        paramlabels = paramlabels, 
        evals_per_sample = 3,
        continue_from = continue_from
    )

    @info "#### ---- Best fit ---- ####"

    let plt = f.plot_data()

        p_opt = f.samples[:,argmin(vec(f.losses))]
        sim_opt = [f.simulator(p_opt) for _ in 1:100]

        save_sims(sim_opt, savetag, "VPC_bestfit")
        plot_sims!(plt, sim_opt, label = "Best fit")
        display(plt)

        if !isnothing(savetag)
            savefig(plot(plt, dpi = 300), datadir("sims", savetag, "VPC_bestfit.png"))
        end
 
    end

    @info "#### ---- Maximum a posteriori estimation ---- ####"

    let plt = f.plot_data()

        p_MAP = f.samples[:,argmax(vec(f.weights))]
        sim_MAP = [f.simulator(p_MAP) for _ in 1:100]

        plot_sims!(plt, sim_MAP, label = "MAP")
        display(plt)
        
        if !isnothing(savetag)
            savefig(plot(plt, dpi = 400), datadir("sims", "$(savetag)", "VPC_MAP.png"))
        end

    end

    @info "#### ---- Posterior samples ---- ####"

    let plt = f.plot_data()
        @suppress_err begin
            global posterior_check = posterior_predictions(f, n_posterior_check)
        end

        save_sims(posterior_check.predictions, savetag, "VPC_posterior")
        plot_sims!(plt, posterior_check.predictions, label = "Retrodictions")

        if !isnothing(savetag)
            savefig(plot(plt, dpi = 400), datadir("sims", savetag, "VPC_posterior.png"))
        end

        display(plt)
    end

    @info "#### ---- Posterior summary ---- ####"

    display(pmchist.posterior_summary)

    
    @info "#### ---- Marginal posteriors ---- ####"

    let plt, num_params = length(f.prior.dists), num_cols = 4, num_rows = Int(ceil(num_params / num_cols))

        plt = plot(
            plot.(f.prior.dists, color = :black)..., layout = (num_rows, num_cols), 
            leg = hcat(vcat(true, repeat([false], num_params-1))...), 
            label = "Prior", 
            size = (1200,200*num_rows), bottommargin = 5mm
            )
    
            
        for (i,param) in enumerate(f.prior.labels)
    
            histogram!(
                plt, subplot = i, 
                f.samples[i,:], weights = Weights(vec(f.weights)), 
                xlabel = param in keys(paramlabels) ? paramlabels[param] : param, 
                normalize = :pdf, label = "Posterior", color = :gray, lw = 0.5, fillalpha = .5
                )
        end
    
        display(plt)
        
        if !isnothing(savetag)
            savefig(plot(plt, dpi = 400), datadir("sims", "$(savetag)", "marginal_posteriors.png"))
        end
    end

    @info "#### ---- Convergence behaviour ---- ####"

    let plt
        plt = plot(eachindex(pmchist.dists) .- 1, map(median, pmchist.dists), marker = true, lw = 1.2, xlabel = "PMC step", ylabel = "Loss", label = "Median", xticks = eachindex(pmchist.dists) .- 1)
        plot!(plt, eachindex(pmchist.dists) .- 1, map(minimum, pmchist.dists), marker = true, lw = 1.2, label = "Minimum")
        display(plt)
        savefig(plot(plt, dpi = 400), datadir("sims", savetag, "loss.png"))
    end

    begin
        posterior_estimates = [mapslices(x -> x[argmax(Weights(pmchist.weights[i]))], pmchist.particles[i], dims = 2) for i in eachindex(pmchist.particles)] |>
        x -> hcat(x...)

        posterior_medians = [mapslices(x -> median(x, Weights(pmchist.weights[i])), pmchist.particles[i], dims = 2) for i in eachindex(pmchist.particles)] |>
        x -> hcat(x...)

        num_params = length(f.prior.dists) 
        num_cols = 4 
        num_rows = Int(ceil(num_params / num_cols))

        plt = plot(
            eachindex(pmchist.particles) .-1,
            posterior_estimates', layout = size(posterior_estimates)[1], 
            size = (1200,200*num_rows), marker = true,
            leg = hcat(vcat(:topleft, repeat([false], length(f.prior.dists)-1))...),
            label = "MAP",
            ylabel = hcat(f.prior.labels...), titlefontsize = 12,
            bottommargin = 5mm, leftmargin = 5mm, 
            xlabel = "PMC step"
            )

        plot!(eachindex(pmchist.particles) .-1, posterior_medians', marker = :diamond, label = "Median")

        for (i,dist) in enumerate(f.prior.dists)

            # set ylim based in prior limits

            q1 = quantile(dist, 0.01)
            q2 = quantile(dist, 0.99)

            # indicate IQR of priors 

            l = repeat([quantile(dist, 0.25)], length(pmchist.particles))
            u = repeat([quantile(dist, 0.75)], length(pmchist.particles))

            plot!(plt, subplot = i, ylim = (q1,q2))

            plot!(
                plt, subplot = i, 
                eachindex(pmchist.particles) .- 1, 
                l,
                fillrange = u, 
                fillalpha = 0.35, 
                color = :lightgray, 
                label = "Prior IQR"
                )

            hline!(plt, subplot = i, [quantile(dist, 0.25)], color = :gray, linestyle = :dash, label = "")
            hline!(plt, subplot = i, [quantile(dist, 0.75)], color = :gray, linestyle = :dash, label = "")

        end
    end

    display(plt)

    @info "#### ---- Posterior summary ---- ####"

    generate_posterior_summary(
        f.prior.labels,
        f.samples, 
        f.losses, 
        f.weights; 
        tex = false,
        paramlabels = paramlabels,
        savetag = nothing
    ) |> display



    return pmchist, posterior_check
end

function plot_pmc_gradient(pmchist)


    plt = plot(layout = (1,2))

    loss_median = [median(x) for x in pmchist.dists]
    loss_minimum = [minimum(x) for x in pmchist.dists]
    loss_median_gradient = []
    loss_minimum_gradient = []
    
    for i in eachindex(loss_median)
        if i==1
            push!(loss_median_gradient, NaN)
            push!(loss_minimum_gradient, NaN)
        else
            push!(loss_median_gradient, -(loss_median[i]-loss_median[i-1])/loss_median[1])
            push!(loss_minimum_gradient, -(loss_minimum[i]-loss_minimum[i-1])/loss_minimum[1])
        end
    end
    
    x = eachindex(loss_median) .- 1 
    plot!(x, loss_median, marker = true, label = "Median", xlabel = "SMC population", ylabel = "Loss", subplot = 1)
    plot!(x, loss_minimum, marker = true, label = "Minimum", subplot = 1)
    
    plot!(x, loss_median_gradient, marker = true, subplot = 2, label = "Median", xlabel = "SMC population", ylabel = "Normalized gradient", leg = false)
    plot!(x, loss_minimum_gradient, marker = true, subplot = 2, label = "Minimum")

    return plt
end


function plot_posteriors(f::ModelFit, pmchist)
    plt = plot(layout = (2, 4), size = (1000,500), bottommargin = 7.5mm)

    for t in eachindex(pmchist.particles)
        for k in eachindex(f.prior.dists)

            if t == 1 
                if typeof(f.prior.dists[k]) == Hyperdist
                    plot!(f.prior.dists[k].dist, subplot = k, color = :black, lw = 2, linestyle = :dash, label = "Prior")
                else
                    plot!(f.prior.dists[k], subplot = k, color = :black, lw = 2, linestyle = :dash, label = "Prior")
                end
            end

            density!(
                plt, subplot = k,
                pmchist.particles[t][k,:], weights = Weights(pmchist.weights[t]), 
                leg = k == 1 ? :topright : false, 
                fill = true, fillalpha = .15, lw = 1, 
                label = "t = $t",
                xlabel = paramlabels[f.prior.labels[k]], xlabelfontsize = 14, 
                yticks = [], yaxis = false
                )
        end
    end

    return plt
end
