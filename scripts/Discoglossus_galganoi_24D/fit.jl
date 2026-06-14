
include("data.jl")
include("parameters.jl")
include("simulation.jl")

import EcotoxModelFitting: loss_mse

function setup_modelfit(pmoa::AbstractString; sigma_factor = 2.)

    pmoa_idx = findfirst(x -> x == pmoa, PMOAS) # convert pmoa from string to index
    @assert !isnothing(pmoa_idx) "Did not find PMoA $(pmoa) in PMOAS"

    if !isdir(datadir("sims", SAVETAG))
        mkdir(datadir("sims", SAVETAG))
    end
    # set up logger
    # FIXME: this does not save the log file to the intended directory, because the setup_modelfit does not know the composite "savetag" 
    global io = open(datadir("sims", SAVETAG, "log.txt"), "w+")
    global logger = SimpleLogger(io)

    data = load_data()

    f = ModelFit( 
        prior = Prior(
            "spc.KD[1,$(pmoa_idx)]" => Truncated(Normal(1, 0.1), 0.001, 1), 
            "spc.E[1,$(pmoa_idx)]" => Truncated(Normal(15, 15*sigma_factor), 0.03, Inf),
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
        grouping_vars = [
            [:C_W_1], 
            [:C_W_1]
            ], 
        time_resolved = [true, false],
        time_var = :t_exp,
        plot_data = plot_data,
        loss_functions = EcotoxModelFitting.loss_euclidean#loss_mse_logtransform
    )

    return f

end

function save_sims(
    predictions::AbstractVector, 
    savetag::AbstractString, 
    prefix::AbstractString
    )::Nothing

    for key in keys(predictions[1])
        df = vcat([@transform(p[key], :num_sample = i) for (i,p) in enumerate(predictions)]...)
        CSV.write(datadir("sims", savetag, "$(prefix)_$(key).csv"), df)
    end

    return nothing
end


function fit_model!(
    f; 
    pmcsettings = (
        :n => 100,
        :q_dist => 0.1, 
        :t_max => 0,
        :evals_per_sample => 1
    ),
    savedir = datadir("sims"),
    savetag = nothing, 
    paramlabels = paramlabels, 
    continue_from = nothing,
    n_posterior_check = 100,
    plot_sims! = plot_sims!
    )

    if !isdir(plotsdir(savetag))
        mkdir(plotsdir(savetag))
    end

    @info "#### ---- Executing PMC --- ####"

    pmchist = run_PMC!(
        f; 
        pmcsettings...,
        savedir = savedir,
        savetag = savetag, 
        paramlabels = paramlabels, 
        continue_from = continue_from
    )

    @info "#### ---- Best fit ---- ####"

    let plt = f.plot_data()

        p_opt = f.accepted[:,argmin(vec(f.losses))]
        sim_opt = [f.simulator(p_opt) for _ in 1:100]

        save_sims(sim_opt, savetag, "VPC_bestfit")
        plot_sims!(plt, sim_opt, label = "Best fit")
        display(plt)

        if !isnothing(savetag)
            savefig(plot(plt, dpi = 300), datadir("sims", savetag, "VPC_bestfit.png"))
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
    
    @info "#### ---- Marginal posteriors ---- ####"

    let plt, num_params = length(f.prior.dists), num_cols = 4, num_rows = Int(ceil(num_params / num_cols))

        plt = plot(
            plot.(f.prior.dists, color = :black)..., layout = (num_rows, num_cols), 
            leg = hcat(vcat(true, repeat([false], num_params-1))...), 
            label = "Prior", 
            size = (1200,200*num_rows), bottommargin = 7.5mm
            )
    
            
        for (i,param) in enumerate(f.prior.labels)
    
            histogram!(
                plt, subplot = i, 
                f.accepted[i,:], weights = Weights(vec(f.weights)), 
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
        posterior_estimates = [mapslices(x -> x[argmin(vec(pmchist.dists[i]))], pmchist.particles[i], dims = 2) for i in eachindex(pmchist.particles)] |>
        x -> hcat(x...)

        posterior_medians = [mapslices(x -> median(x, Weights(pmchist.weights[i])), pmchist.particles[i], dims = 2) for i in eachindex(pmchist.particles)] |>
        x -> hcat(x...)

        num_params = length(f.prior.dists) 
        num_cols = 4 
        num_rows = Int(ceil(num_params / num_cols))

        plt = plot(
            eachindex(pmchist.particles) .- 1,
            posterior_estimates', layout = (num_rows,num_cols), 
            size = (1200,300*num_rows), marker = true,
            leg = hcat(vcat(:topleft, repeat([false], length(f.prior.dists)-1))...),
            label = "Best fit",
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
        f; 
        tex = false,
        paramlabels = paramlabels,
        savetag = nothing
    ) |> display



    return pmchist, posterior_check
end


function fit_all_pmoas()
    
    @info "### ---- Fitting models for PMoAs: $(PMOAS) ---- ####"

    let fs
        @suppress fs = setup_modelfit.(PMOAS)

        for (i,f) in enumerate(fs)
            @info "#### ---- Running model fit for PMoa $(PMOAS[i]) ---- ####"
            _= fit_model!(
                f; 
                savetag = "$(SAVETAG)_$(PMOAS[i])",
                n_init = 5_000, 
                n = 5_000, 
                t_max = 3, 
                q_dist = 0.1
            )
        
        end
        return fs
    end
end

