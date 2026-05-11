using Distributions
using DataFrames
using ProgressMeter
using DataStructures
using StatsBase
using LaTeXStrings
using Base.Threads
import Base: rand
import Base: getindex
import Base: setindex!
import Base:show


includet(srcdir("loss.jl")) # definitions of basic loss functions
includet(srcdir("utils.jl"))

# reserved column names for the posterior -> cannot be used as parameter names
RESERVED_COLNAMES = ["loss", "weight", "model", "chain"]

abstract type AbstractPrior end

"""
    Hyperdist

A mutable struct for hyper-distributions in multi-level modelling. 
"""
mutable struct Hyperdist
    gendist::Function
    dist::Distribution
end

import Plots:plot
plot(hyper::Hyperdist; kwargs...) = plot(hyper.dist; kwargs...)

import Distributions: mode, mean, median, std, var, minimum, maximum, pdf, quantile
mode(hyper::Hyperdist) = mode(hyper.dist)
mean(hyper::Hyperdist) = mean(hyper.dist)
median(hyper::Hyperdist) = median(hyper.dist)
std(hyper::Hyperdist) = std(hyper.dist)
var(hyper::Hyperdist) = var(hyper.dist)
minimum(hyper::Hyperdist) = minimum(hyper.dist)
maximum(hyper::Hyperdist) = maximum(hyper.dist)
pdf(hyper::Hyperdist, x) = pdf(hyper.dist, x)
quantile(hyper::Hyperdist, q::Float64) = quantile(hyper.dist, q)

rand(hyper::Hyperdist) = rand(hyper.dist)


# transform truncated normal distribution to unit scale
function scaledist(dist::Truncated{Normal{Float64}})

    μ, σ = dist.untruncated.μ, dist.untruncated.σ
    l, u = dist.lower, dist.upper

    l_scaled = (l - μ) / σ
    u_scaled = (u - μ) / σ
    
    dist_scaled = truncated(Normal(), l_scaled, u_scaled)

    return dist_scaled, μ, σ
end

# transform Dirac distribution to unit scale 
function scaledist(dist::Dirac)
    μ = dist.value
    return Dirac(0.), μ, 1
end

# transform Discrete Uniform distribution to unit scale
function scaledist(dist::DiscreteUniform)

    μ = mean(dist)
    σ = (maximum(dist) - minimum(dist))/2 # using half the range leads to correct scaling (using std would be wrong in this case)

    return Uniform(-1, 1), μ, σ
end

scale_param(x::Real, μ::Real, σ::Real) = (x-μ)/σ

function scaledist(hyper::Hyperdist)
    return scaledist(hyper.dist)
end

mutable struct Prior <: AbstractPrior

    labels::Vector{String}
    dists::Vector{Union{Hyperdist,Distribution}}
    gendists::Vector{Function}
    is_hyper::Vector{Bool}
    scaled_dists::Vector{Distribution}
    μs::Vector{Float64}
    σs::Vector{Float64}
    
    """
        Prior(args::Pair...)
        
    Initialize prior instance with a sequence of String/Distribution pairs.

    Distributions can either be a plain distribution as defined in `Distributions.jl` or an instance of `Hyperdist` for multi-level inference 
    (e.g. to estimate the spread of individual variability).

    Example: 

    ```Julia
    prior = Prior(
        :a => truncated(Normal(1, 1), 0, Inf), 
        :b => Beta(1, 1),
    )
    ```

    """
    function Prior(args::Pair...)

        labels = String[]
        dists = Union{Distribution,Hyperdist}[]
        gendists = Function[]
        is_hyper = Bool[]
        scaled_dists = Distribution[]
        μs = Float64[]
        σs = Float64[]

        for (pair) in args

            #@assert !(pair.first in RESERVED_COLNAMES) "The following names are reserved and cannot be used for parameters: $(RESERVED_COLNAMES)"

            scaled_dist, μ, σ = scaledist(pair.second)

            push!(labels, pair.first)
            push!(dists, pair.second)
            push!(scaled_dists, scaled_dist)
            push!(μs, μ)
            push!(σs, σ)
            
            if typeof(pair.second) != Hyperdist
                push!(is_hyper, false)
            else
                push!(gendists, pair.second.gendist)
                push!(is_hyper, true)
            end
        end

        return new(
            labels, 
            dists, 
            gendists, 
            is_hyper, 
            scaled_dists, 
            μs, 
            σs
            )
    end

    """
        Prior(params, prior)

    Initialize prior from a vector of parameter names and prior distributions, respectively.
    """
    function Prior(params, prior)
        return new(params, prior)
    end
end

function show(prior::Prior)
    return OrderedDict(zip(prior.labels, prior.dists))
end

function getindex(prior::Prior, param::Union{String,Symbol})

    index = findfirst(isequal(param), prior.labels)
    @assert index !== nothing "Parameter $param not found in prior object"
    
    return prior.dists[index]
end

function setindex!(prior::Prior, value::Union{Distribution,Hyperdist}, param::Union{String,Symbol})
    
    index = findfirst(isequal(param), prior.labels)
    @assert index !== nothing "Parameter $param not found in prior object"

    prior.dists[index] = value

end


mutable struct ModelFit
    
    prior::Prior
    defaultparams#::EcotoxSystems.ComponentArray
    psim#::EcotoxSystems.ComponentArray
    simulator::Function
    loss::Function
    loss_functions::AbstractVector
    data::OrderedDict
    time_var::Symbol
    response_vars::Vector{Vector{Symbol}}
    grouping_vars::Vector{Vector{Symbol}}
    data_weights::Vector{Vector{Float64}}
    time_resolved::Vector{Bool}
    plot_data::Function
    losses::Matrix{Float64}
    samples::Matrix{Float64}
    thresholds::Vector{Float64}
    weights::Vector{Float64}
    combine_losses::Function

    function ModelFit(;
        prior::Prior, 
        defaultparams::EcotoxSystems.ComponentArray,
        simulator::Function, 
        data::Any, 
        response_vars::Vector{Vector{Symbol}},
        time_resolved::Vector{Bool},
        data_weights::Union{Nothing,Vector{Vector{Float64}}} = nothing,
        time_var::Union{Nothing,Symbol} = nothing,
        grouping_vars::Union{Nothing,Vector{Vector{Symbol}}} = nothing,
        plot_data::Function = function plot_data() plot() end,
        loss_functions::Union{Vector{Vector{Function}},Function} = loss_mse, 
        combine_losses::Function = sum # function that combines losses into single value; can also be x->x to retain all losses
        )

        f = new()

        f.prior = prior

        @info "
        Set up for estimating(length(prior.dists)) parameters: 
        $(f.prior.labels)
        "
        f.defaultparams = deepcopy(defaultparams)
        f.psim = [deepcopy(f.defaultparams) for _ in 1:Threads.nthreads()]
        f.time_var = time_var
        f.response_vars = response_vars

        if isnothing(grouping_vars)
            f.grouping_vars = repeat([Symbol[]], length(data))
        else
            f.grouping_vars = grouping_vars
        end

        if isnothing(data_weights)
            @info "No data weights provided, assuming uniform weights"
            f.data_weights = [ones(size(v)) for v in f.response_vars] |> x-> x ./ sum(vcat(x...))
        else
            @info "Normalizing data weights"
            f.data_weights = data_weights ./ sum(vcat(data_weights...))
        end

        f.time_resolved = time_resolved
        f.simulator = generate_fitting_simulator(defaultparams, prior, simulator)
        f.data = data
        f.plot_data = plot_data
        assign_loss_functions!(f, loss_functions)
        f.combine_losses = combine_losses
        f.loss = generate_loss_function(f)

        return f
    end
end


function update_data_weights!(f::ModelFit, data_weights::Vector{Vector{R}})::Nothing where R <: Real
    
    f.data_weights = data_weights |> x -> x ./ sum(vcat(x...))
    f.loss = generate_loss_function(f)

    return nothing
end

# separate loss function for each response variable
function assign_loss_functions!(f::ModelFit, loss_functions::Vector{Vector{F}}) where F <: Function
    f.loss_functions = loss_functions
end

# single loss function is given => assume same for all response variables
function assign_loss_functions!(f::ModelFit, loss_functions::F) where F <: Function
    f.loss_functions = [repeat([loss_functions], length(vars)) for vars in f.response_vars]
end

"""
    generate_fitting_simulator(defaultparams::ComponentArray, prior::Prior, simulator::Function)::Function

Attempt to define a generic simulator function, 
based on the information given to the ModelFitting object. <br> 
The generated "fitting simulator" 
    - Is a wrapper around the provided `simulator` argument
    - Expects the parameter values as vector of floats
    - Assures that parameters are assigned correctly to a copy of the defaultparams
    - Pre-allocates copies of the defaultparams with account for multithreading
    - Deals with hyperparameters (currently not in a full hierarchical approach, TBC)
    - Defines a second method for the fitting_simulator which dispatches to the original `simulator` function (useful in conjunction with the `ModelFit` struct)

I am sure there are use-cases where this will fail. For the cases tested so far though, it worked fine and was quite helpful.
"""
function generate_fitting_simulator(
    defaultparams, 
    prior::Prior, 
    simulator::Function
    )::Function

    # when using mult-threading, we create a copy of the parameter object for each thread
    pfit = [deepcopy(defaultparams) for _ in 1:Threads.nthreads()]

    # matching parameter labels to indices
    pfit_labels = EcotoxSystems.ComponentArrays.labels(pfit[1])
    idxs = [findfirst(x -> x == l, pfit_labels) for l in prior.labels]

    function fitting_simulator(pvec::Vector{R}; kwargs...) where R <: Real
        
        psim = pfit[threadid()] # pick the parameter copy for the current thread

        psim[idxs[.!prior.is_hyper]] = pvec[.!prior.is_hyper] # assign "normal" parameters directly
        psim[idxs[prior.is_hyper]] = [gendist(h) for (gendist,h) in zip(prior.gendists, pvec[prior.is_hyper])] # assign hyperparameters through the appropriate gendist function

        return simulator(psim; kwargs...)
    end

    # define an additional method that takes the componentarray as argument, and dispatches to the original function
    fitting_simulator(p::CA; kwargs...) where CA <: EcotoxSystems.ComponentArray = simulator(p; kwargs...)

    return fitting_simulator

end

"""
    define_data_scales(f::ModelFit, join_vars::Vector{Vector{Symbol}}; strategy = "maxima")::Vector{Vector{Float64}}

Defines data scales based on strategy "maxima".


Strategy "var2" is still to be implemented.
"var2" uses twice the mean-centered variance of the data, 
where mean-centering is done for each `join_var`. 

(`join_vars` are `grouping_vars`+`time_var`).

In combination with `loss=loss_mse`, the resulting loss function is the simplified normal likelihood, 
corrected for the number of observations per response and grouping variable.
"""
function define_data_scales(
    f::ModelFit, 
    join_vars::Vector{Vector{Symbol}}; 
    strategy = "maxima"
    )::Vector{Vector{Float64}}

    if strategy == "maxima"
        scales = [zeros(size(vars)) for vars in f.response_vars]

        for (i,key) in enumerate(f.data.keys)
            for (j,var) in enumerate(f.response_vars[i])
                scale = maximum(skipmissing(f.data[key][:,var]))
                scales[i][j] = scale
            end
        end

        return scales

    else 
        error("Strategy $strategy not implemented or work in progress")

        #scales = [zeros(size(vars)) for vars in f.response_vars]
        #
        ## for every data table
        #for (i,key) in enumerate(f.data.keys)
        #    # if we have no join vars (i.e. no time variable and no additional grouping vars), 
        #    # then we do the mean-centering over all values for each response variable
        #    if length(join_vars[i]) == 0
        #        for (j,y) in enumerate(f.response_vars[i])
        #            yvals = f.data[key][:,y]
        #            if length(yvals)>1
        #                scales[i][j] = 2*var(yvals .- mean(yvals))
        #            # if we have only one value, then the scale is the squared value itself
        #            else
        #                scales[i][j] = yvals[1]^2
        #            end
        #        end
        #    # if we have join vars (time, grouping vars or both), 
        #    # then we do the mean-centering for each value of each join var
        #    else
        #        for (j,y) in enumerate(f.response_vars[i])
        #            centered_yvals = []
        #            all_yvals = []
        #            for (m,jvar) in enumerate(join_vars[i])
        #                for (l,jvar_val) in enumerate(unique(f.data[key][:,jvar]))
        #                    df = f.data[key]
        #                    yvals = df[df[:,jvar].==jvar_val,y]
        #                    yvals_mean = mean(skipmissing(yvals))
        #                    push!(centered_yvals,  yvals .- mean(yvals))
        #                    push!(all_yvals, yvals)
        #                end
        #            end
        #
        #            if length(unique(centered_yvals))>1
        #                @info """For $key data and repsonse variable $y, found no variance within groups to calculate error variance. Using mean(y)^2 as data scale."""
        #                scales[i][j] = 2*var(centered_yvals)
        #            else
        #                scales[i][j] = mean(all_yvals)^2
        #            end
        #        end
        #    end
        #end

    end

end

"""
    generate_loss_function(f::ModelFit)::Function

Generates a loss function based on some simplifying assumptions: 

    - All data is stored in a dictionary of `DataFrame`s (data tables).
    - The simulation output is give in the same format.
    - Eeach data table is either time-resolved or not.
    - If the data is time-resolved, it has to have a column whose name is indicated by `time_var`.
    - `f.loss_functions` lists the error models applied for each response variable. `loss_mse` by default.
    - Each data table can have multiple response variables, indicated in `f.response_vars`

By default, the individual losses for each response variable are returned separately. <br>
If `run_ABC!` is used, they will be converted to an independent acceptance probability. <br>
Otherwise, they have to be combined into a single loss. 
"""
function generate_loss_function(f::ModelFit)::Function

    # if necessary, add time to grouping variables
    # TODO: this won't work with multinom like or dtw
    # in that case we need to differentiate inside the loss function

    join_vars = similar(f.grouping_vars) # columns used to join data frames
    data_columns = similar(f.grouping_vars) # names of all columns contained in the data; relevant for implementing mutating joins
    select_cols = similar(f.response_vars) # names of eval_df columns which have to be selected before applying dropmissing

    # for every data table
    for (i,key) in enumerate(f.data.keys)
        # check wether we have time-resolved data
        data_columns[i] = Symbol.(names(f.data[key]))
        if f.time_resolved[i]
            # if so, add time to the grouping variables 
            join_vars[i] = vcat(f.time_var, f.grouping_vars[i]) |> unique
        else
            # otherwise, don't
            join_vars[i] = f.grouping_vars[i] |> unique
        end
        
        select_cols[i] = vcat(
            f.grouping_vars[i],
            f.response_vars[i],
            [Symbol("$(x)_1") for x in f.response_vars[i]]
        )
    end

    data_scales = define_data_scales(f, join_vars)

    # get the "nominal length", number of observations in each data table
    nominal_lenghts = [nrow(dropmissing(df)) for df in values(f.data)]
    

    # this function will compute a separate loss for each response variable, 
    # assuming that their acceptance probabilites will be combined later

    function lfun(data::OrderedDict, sim::OrderedDict)::Union{Float64,Vector{Float64}}

        # allocate losses as 1-D Vector across data tables and response vars 
        losses = Vector{Float64}(undef, length(vcat(f.response_vars...)))

        idx = 0

        # for each data table  
        for (i,key) in enumerate(keys(data)) 
            # if we have grouping variables, match data with leftjoin
            if length(join_vars[i])>0
                # merge the corresponding data with prediction 
                eval_df = leftjoin(
                    data[key], sim[key], 
                    on = join_vars[i], 
                    makeunique = true
                    ) |> 
                    x -> select(x, select_cols[i]) |>
                    dropmissing

                # for each response variable in the table   
                for (j,var) in enumerate(f.response_vars[i]) 
                    idx += 1
                    scale = data_scales[i][j]
                    # if data is time-resolved or has grouping variables, but does not have both 
                    # we calculate the loss over all values for the jth response in the ith table
                    if ((f.time_resolved[i]) || (length(f.grouping_vars[i])>0)) && (!((f.time_resolved[i])&&(length(f.grouping_vars[i])>0)))
                        losses[idx] = f.loss_functions[i][j]( # take the ith loss function
                            Vector{Float64}(eval_df[:,Symbol("$(var)")]) ./ scale, 
                            Vector{Float64}(eval_df[:,Symbol("$(var)_1")]) ./ scale,
                            f.data_weights[i][j],
                            nominal_lenghts[i]
                            )

                    # if data is time-resolved and has additional grouping variables, 
                    # we have to calculuate the loss for each grouping variable and then sum up the losses
                    # this is necessary so that we don't mix up different time series when the error model takes temporal dependency into account
                    # (e.g. multinomial likelihood or dynamic time warping)
                    else
                        losses[idx] = begin
                            groupby(eval_df, f.grouping_vars[i]) |> 
                            # this can be accelerated by replacing the do-syntax
                            x -> combine(x) do df
                                DataFrame(
                                    loss = f.loss_functions[i][j](
                                        Vector{Float64}(df[:,Symbol("$(var)")]) ./ scale, 
                                        Vector{Float64}(df[:,Symbol("$(var)_1")]) ./ scale,
                                        f.data_weights[i][j],
                                        nominal_lenghts[i]
                                    )
                                )
                            end.loss |> 
                            sum
                        end
                    end
                end
            # if data is not time-resolved and has no grouping variables, we don't need a join operation
            # NOTE: this works when data.var and sim.var have the same length, 
            # or the length of one of them is 1
            # if we have different lengths and both are > 1, we need a different strategy, 
            # but I'm not sure if that is a case that should ever occur
            else
                for (j,var) in enumerate(f.response_vars[i])
                    idx += 1
                    scale = data_scales[i][j]
                    losses[idx] = f.loss_functions[i][j](
                        data[key][:,var] ./ scale, 
                        sim[key][:,var] ./ scale,
                        f.data_weights[i][j]
                        )
                end
            end
        end

        return f.combine_losses(losses)
    end

    # if the simulation throws an error, the result will be `nothing` => return infinite loss
    # the ABC algorithm will remove these before calculating the rejection threshold
    # for optimization algorithms, this is not guaranteed to work as intended (depends on the implementation)

    lfun(data::Any, sim::Nothing) = Inf

    return lfun
   
end

#### Functions for rejection ABC ####

function rand(prior::Prior)
    return [rand(p) for p in prior.dists]
end

"""
    posterior_sample(accepted::DataFrame; reserved_colnames::Vector{String} = RESERVED_COLNAMES)

Take posterior sample from a data frame of accepted values.
"""
function posterior_sample(accepted::DataFrame; reserved_colnames::Vector{String} = RESERVED_COLNAMES)::Vector{Float64}
    ω =  accepted.weight
    selectcols = filter(x -> !(x in reserved_colnames), names(accepted)) 
    sampled_values = accepted[sample(axes(accepted, 1), Weights(ω)),selectcols]
    return Vector{Float64}(sampled_values)
end

#posterior_sample(res::SMCResult; kwargs...) = posterior_sample(res.accepted; kwargs...)

# "old" posterior sampling method based on dataframes. maybe still need it.
function posterior_sample!(
    p::Any, 
    accepted::DataFrame; 
    reserved_colnames::Vector{String} = RESERVED_COLNAMES
    )

    ω =  accepted.weight
    selectcols = filter(x -> !(x in reserved_colnames), names(accepted))
    sampled_values = accepted[sample(axes(accepted, 1), Weights(ω)),selectcols]
    param_names = names(sampled_values)
    assign!.(p, param_names, sampled_values)
end

function posterior_sample(
    samples::Matrix{Float64}, 
    weights::Vector{Float64}
    )
    return samples[:,sample(1:size(f.samples)[2], Weights(weights))]
end

# dispatches to Matrix method
function posterior_sample(f::ModelFit)
    return posterior_sample(f.samples, f.weights)
end

"""
    bestfit(defparams::AbstractParams, accepted::AbstractDataFrame)
Get the best f from `accepted` (particle with minimum loss) and assign to a copy of `defparams`.
"""
function bestfit(accepted::AbstractDataFrame)
    return posterior_sample(accepted[accepted.loss.==minimum(accepted.loss),:])
end

function prior_predictive_check(
    f::ModelFit;
    compute_loss::Bool = true,
    loss = f.loss,
    n::Int64 = 100
    )::NamedTuple

    losses = Vector{Union{Float64,Vector{Float64}}}(undef, n)
    predictions = Vector{Any}(undef,n)
    samples = Vector{Vector{Float64}}(undef, n)

    @info "#### ---- Evaluating $n prior samples on $(Threads.nthreads()) threads ---- ####"

    @showprogress @threads for i in 1:n
        
        prior_sample = rand(f.prior)
        prediction = f.simulator(prior_sample)

        L = NaN

        if compute_loss
            L = loss(f.data, prediction)
        end

        predictions[i] = prediction
        losses[i] = L
        samples[i] = prior_sample

    end

    return (
        predictions = predictions,
        losses = losses,
        samples = samples
    )
end

function kernel(
    dists::AbstractVector, 
    thresholds::Vector{Float64}
    )

    w = Vector{Float64}(undef, length(thresholds))
    for i in eachindex(thresholds)
        dist = dists[i]
        threshold = thresholds[i]
        if dist > threshold
            w[i] = 0
        else
            w[i] = 1/threshold * (1 - (dist/threshold)^2)
        end
    end

    return w
end

compute_thresholds(q_dist::Float64, losses::Matrix{Float64}) = [quantile(d, q_dist) for d in eachrow(losses)]

"""
    compute_weights(losses::Matrix{Float64}, thresholds::Vector{Float64})::Vector{Float64}

Compute sampling weights for ABC rejection sampling  from `losses` and `thresholds`. <br>
`losses` is 2-dimensional and has a row for each response variable and a column for each evaluated sample. <br>
`thresholds` is 1-dimensional and has a cut-off value for the loss for each response variable.  
"""
function compute_weights(losses::Matrix{Float64}, thresholds::Vector{Float64})::Vector{Float64}

    # apply kernel function for each column in losses, 
    # i.e., for each response variable (endpoint, summary statistics, whatchamacallit...)

    weights = [kernel(d, thresholds) for d in eachcol(losses)] |> 
    x -> hcat(x...)

    # normalize weights for each response variable
    for i in 1:size(weights)[1]
        weights[i,:] ./ sum(weights[i,:])
    end

    # combine weights over response variables by multiplication
    weights_combined = [prod(w_i) for w_i in eachcol(weights)]

    return weights_combined
end

function apply_rejection(q_dist::Float64, losses::Matrix{Float64})

    thresholds = compute_thresholds(q_dist, losses)
    weights = compute_weights(losses, thresholds)

    return weights

end

function apply_rejection!(f::ModelFit, q_dist::Float64)::Nothing
    
    f.thresholds = compute_thresholds(q_dist, f.losses)
    f.weights = compute_weights(f.losses, f.thresholds)

    return nothing
end

"""
    run_ABC!(
        f::ModelFit; 
        n::Int = 1000,
        q_dist::Float64 = .1
        )::Nothing

Execute parameter inference with basic rejection ABC.
"""
function run_ABC!(
    f::ModelFit; 
    loss = f.loss, 
    n::Int = 1000,
    q_dist::Float64 = .1,
    append::Bool = false,
    savetag::Union{Nothing,String} = nothing
    )::Nothing

    @info "#### ---- Evaluating $n ABC samples on $(Threads.nthreads()) threads ---- ####"
        
    losses = Vector{Union{Float64,Vector{Float64}}}(undef, n)
    samples = Vector{Vector{Float64}}(undef, n)

    @showprogress @threads for i in 1:n
        prior_sample = rand(f.prior)
        
        prediction = f.simulator(prior_sample)
        L = loss(f.data, prediction)

        losses[i] = L
        samples[i] = prior_sample
    end

    @info "#### ---- Computing weights ---- ####"

    if !append
        f.samples = hcat(samples...)
        f.losses = hcat(losses...)
    else
        f.samples = hcat(f.samples, hcat(samples...))
        f.losses = hcat(f.losses, hcat(losses...))
    end

    valid_losses = [sum(isinf.(x) .|| isnan.(x))==0 for x in eachcol(f.losses)]

    @info "Retained $(sum(1 .- valid_losses)) valid samples in $n total samples"

    f.samples = f.samples[:,valid_losses]
    f.losses = f.losses[:,valid_losses]

    apply_rejection!(f, q_dist)

    #if !isnothing(savetag)
    #    if !isdir(datadir("sims", savetag))
    #        mkdir(datadir("sims", savetag))
    #    end
    #
    #    if !append
    #        CSV.write(datadir("sims", savetag, "samples.csv"))
    #    end
    #
    #end

    return nothing
end

norm(x) = x ./ sum(x)

#prod(pdf.(f.prior.dists, θ_i))

"""
    prior_prob(prior::Prior, theta::AbstractVector)

This function calculates the prior probability of a particle, using the unit-scaled prior distributions and particle. 
"""
function prior_prob(prior::Prior, theta::AbstractVector)

    probs = pdf.(
        prior.scaled_dists, 
        scale_param.(theta, prior.μs, prior.σs)
        )

    return prod(probs)

end



# Define a logging function that writes messages to a file
function setup_logging(log_file_name)
    # Open the log file for appending
    log_file = open(log_file_name, "a")

    # Define a file logging backend
    file_logger = FileLogger(log_file)

    # Create a logger with global level of Debug to capture all levels, especially Warning
    custom_logger = MinLevelLogger(file_logger, Logging.Debug)

    # Set custom logger as current logger
    global_logger(custom_logger)

    return log_file  # Return the log file so it can be closed later
end

"""
    generate_posterior_summary(
        particles::Matrix{Float64}, 
        losses::Matrix{Float64}, 
        weights::Vector{Float64}; 
        tex = false,
        paramlabels::Union{Nothing,Dict} = nothing,
        savetag::Union{Nothing,String} = nothing
        )::DataFrame

Generate summary of marginal posterior distributions. 

For `tex=true`, the result will aditionaly be converted to LaTeX and saved as `posterior_summary.tex` within the `savetag` directory. 
This option asssumes that DrWatson is in use and the datadir() function is defined.
If no `savetag` is provided, `tex=true` will be ignored.
"""
function generate_posterior_summary(
    labels::Vector{String},
    particles::Matrix{Float64}, 
    losses::Matrix{Float64}, 
    weights::Vector{Float64}; 
    tex = false,
    paramlabels::Union{Nothing,AbstractDict} = nothing,
    savetag::Union{Nothing,String} = nothing
    )::DataFrame

    if tex & isnothing(savetag)
        tex = false
        "No savetag provided, ignoring tex=true"
    end
        
    best_fit = particles[:,argmin(vec(losses))]
    medians = mapslices(x -> median(x, Weights(weights)), particles, dims = 2) |> vec
    q05 = mapslices(x -> quantile(x, Weights(weights), 0.05), particles, dims=2) |> vec
    q95 = mapslices(x -> quantile(x, Weights(weights), 0.95), particles, dims=2) |> vec


    posterior_summary = DataFrame(
        param = labels,
        best_fit = best_fit, 
        median = medians, 
        q05 = q05, 
        q95 = q95, 

    )

    if tex
        if !isnothing(paramlabels)
            parnames = [p in keys(paramlabels) ? paramlabels[p] : p for p in labels]
            tex_df = @transform(posterior_summary, :param = parnames)
            _df_to_tex(tex_df, datadir("sims", savetag, "posterior_summary.tex"), colnames = ["Parameter", "Best fit", "Median", L"$P_{05}$", L"$P_{95}$"])
        end
    end

    if !isnothing(savetag)
        CSV.write(datadir("sims", savetag, "posterior_summary.csv"), posterior_summary)
    end

    return posterior_summary
    
end



"""
    run_PMC!(
        f::ModelFit; 
        dist = f.loss, 
        n::Int = 1000,
        n_init::Union{Nothing,Int} = nothing,
        q_dist::Float64 = .1,
        t_max = 3,
        savetag::Union{Nothing,String} = nothing,
        continue_from::Union{Nothing,String} = nothing,
    )::NamedTuple

Model fitting with Approximate Bayesian Computation Population Monte Carlo (ABC-PMC) 
as described by Beaumont et al. (2009). 

args

- `f`: A `ModelFit` object, containing simulator, priors, data and loss function. 

kwargs

- `dist`: A distance funtion. Default is `f.loss`. 
- `n`: Number of evaluated samples per population
- `n_init`: Number of evaluated samples in the initial population. The initial population may contain more non-finite losses, so it can make sense to choose `n_init>n`.
- `q_dist`: Distance quantile to determine next acceptance threshold. A lower `q_dist` value leads to more agressive rejection and faster convergence to a solution, with the risk of identifying a local minimum. If all samples return a finite loss, the number of accepted particles is `n*q_eps`. If there are Infs or NaNs in the losses, the number of accepted particles will be lower. 
- `t_max`: Number of PMC updates to perform (not counting initial population). 
- `evals_per_sample`: Number of evaluations to perform for each parameter sample. `evals_per_sample>1` is relevant for stochastic models.
- `savetag`: Tag under which results are saved. Assumes that function `datadir()` is known (provided by `DrWatson`)
- `continue_from`: Path to a checkpoint file from which to continue the fitting. 
- `paramlabels`: Formatted parameter labels used to generate a summary of the posterior distribution as latex table. Labels have to be LaTeX-compatible.  
"""
function run_PMC!(
    f::ModelFit; 
    dist = f.loss, 
    n::Int = 1000,
    n_init::Union{Nothing,Int} = nothing,
    q_dist::Float64 = .1,
    t_max = 3,
    evals_per_sample::Int64 = 1,
    savetag::Union{Nothing,String} = nothing,
    continue_from::Union{Nothing,String} = nothing,
    paramlabels::Union{Nothing,AbstractDict} = nothing
    )::NamedTuple

    t = 0

    all_particles = Matrix{Float64}[]
    all_losses = Matrix{Float64}[]
    all_weights = Vector{Float64}[]
    all_vars = Vector{Float64}[]
    all_thresholds = Vector{Float64}[]

    lower = minimum.(f.prior.dists)
    upper = maximum.(f.prior.dists)

    if isnothing(n_init)
        n_init = n
    end

    if !isnothing(savetag)
        @info "Saving results to $(datadir("sims", savetag))"

        if !isdir(datadir("sims", savetag))
            mkdir(datadir("sims", savetag))
        end
    end
    
    while (t <= t_max) 
        t += 1
        if t == 1       
            # if there is no previous checkpoint given via the continue_from argument, 
            # run the initial population     
            if isnothing(continue_from)
                @info "#### ---- Evaluating $n_init initial samples on $(Threads.nthreads()) threads ---- ####"
                
                particles = Vector{Vector{Float64}}(undef, n_init)
                weights = fill(1/n_init, n_init)
                losses = Vector{Union{Float64,Vector{Float64}}}(undef, n_init)

                @showprogress @threads for i in 1:n_init
                    θ = rand.(f.prior.dists) # take a prior sample
                    L = 0. # initialize loss for this sample
                    for eval in 1:evals_per_sample
                        sim = f.simulator(θ) # run the simulation
                        L += dist(f.data, sim) # add up losses
                    end
                    particles[i] = θ
                    losses[i] = L/evals_per_sample # average the losses retrieved from repeated simulations
                end
                
                particles = hcat(particles...)
                losses = hcat(losses...)

                valid_dist_idxs = [sum(isinf.(x) .|| isnan.(x))==0 for x in eachcol(losses)]
                
                particles = particles[:,valid_dist_idxs]
                losses = losses[:,valid_dist_idxs]
                weights = weights[valid_dist_idxs]

                # apply (possibly multivariate) rejection 

                threshold = compute_thresholds(q_dist, losses)
                accepted_particle_idxs = [sum(d .> threshold)==0 for d in eachcol(losses)]
                
                particles = particles[:,accepted_particle_idxs]
                losses = losses[:,accepted_particle_idxs]
                weights = weights[accepted_particle_idxs] |> norm
                
                # take τ to be twice the empirical variance of Θ
                vars = [2*var(tht) for tht in eachrow(particles)]

                push!(all_particles, particles)
                push!(all_weights, weights)
                push!(all_losses, losses)
                push!(all_vars, vars)
                push!(all_thresholds, threshold)

                # save data to checkpoint
                if !(isnothing(savetag))
                    save(datadir("sims", savetag, "checkpoint.jld2"), Dict(
                        "particles" => all_particles, 
                        "weights" => all_weights, 
                        "losses" => all_losses,
                        "variances" => all_vars, 
                        "thresholds" => all_thresholds,
                        "settings" => Dict(
                            "n" => n, 
                            "n_init" => n_init, 
                            "t_max" => t_max, 
                            "q_dist" => q_dist
                        )
                    ))
                end
            # if we have a previous checkpoint to continue from, load data and continue
            else
                @info "Continuing model fit from checkpoint $(continue_from)"
                chk = load(continue_from)
                all_particles = chk["particles"]
                all_weights = chk["weights"]
                all_losses = chk["losses"]
                all_vars = chk["variances"]
                all_thresholds = chk["thresholds"]
            end
        else
            
            @info "#### ---- PMC step $(t-1)/$(t_max) ---- ####"

            old_particles = all_particles[end]
            old_weights = all_weights[end]
            old_vars = all_vars[end]
            old_vars = max.(1e-100, old_vars)
            
            particles = Vector{Vector{Float64}}(undef, n)
            weights = fill(1/n, n)
            losses = Vector{Union{Float64,Vector{Float64}}}(undef, n)

            @showprogress @threads for i in 1:n

                # take a weighted sample of previously accepted particles
                idx = sample(1:length(old_weights), Weights(old_weights))
                 
                # get associated Θ and weights
                θ_i_ast = old_particles[:,idx]
                ω = old_weights[idx]
                θ_i = similar(θ_i_ast)

                # perturb the particle
                for (k,(tht_k,var_k)) in enumerate(zip(θ_i_ast, old_vars))
                    θ_i[k] = rand(truncated(Normal(tht_k, var_k .+ 1e-100), lower[k], upper[k]))
                end

                # calculate the weight 
                weight_num = prior_prob(f.prior, θ_i) # numerator is the prior probability
                weight_denom = 1e-300 # initialize denominator

                # set ω_it ∝ π(θ)/∑(ω_jt-1 K(θ_i | θ_j, τ^2)} (cf. Beaumont et al. 2009)
                for j in eachindex(old_weights)
                    ω_j = old_weights[j]
                    θ_j = old_particles[:,j]
                    ϕ = prod(pdf.(Normal.(), (θ_j .- θ_i)./old_vars))
                    weight_denom += ω_j * ϕ
                end

                # using log-weights often stabilizes the posterior
                ω = log((weight_num/weight_denom) + 1)

                # run the simulations
                L = 0.
                for eval in 1:evals_per_sample
                    sim = f.simulator(θ_i)
                    L += dist(f.data, sim) # generate ρ(x,y)
                end

                particles[i] = θ_i
                weights[i] = ω
                losses[i] = L/evals_per_sample
            end

            particles = hcat(particles...)
            losses = hcat(losses...)

            valid_dist_idxs = [sum(isinf.(x) .|| isnan.(x))==0 for x in eachcol(losses)]

            particles = particles[:,valid_dist_idxs]
            losses = losses[:,valid_dist_idxs]
            weights = weights[valid_dist_idxs]

            threshold = compute_thresholds(q_dist, losses)
            accepted_particle_idxs = [sum(l .> threshold)==0 for l in eachcol(losses)]
            particles = particles[:,accepted_particle_idxs]
            weights = weights[accepted_particle_idxs] |> norm
            losses = losses[:,accepted_particle_idxs]

            # take τ^2_t+1 as twice the weighted empirical variance of the θ_its 
            vars = [2*var(tht, Weights(weights)) for tht in eachrow(particles)]

            push!(all_particles, particles)
            push!(all_weights, weights)
            push!(all_losses, losses)
            push!(all_vars, vars)
            
            # save data to checkpoint
            if !(isnothing(savetag))
                save(datadir("sims", savetag, "checkpoint.jld2"), Dict(
                    "particles" => all_particles, 
                    "weights" => all_weights, 
                    "losses" => all_losses,
                    "variances" => all_vars, 
                    "thresholds" => all_thresholds,
                    "prior" => f.prior,
                    "settings" => Dict(
                        "n" => n, 
                        "n_init" => n_init, 
                        "t_max" => t_max, 
                        "q_dist" => q_dist
                    )))
            end
        end
    end

    f.samples = all_particles[end]
    f.weights = all_weights[end]
    f.losses = all_losses[end]
    
    if !isnothing(savetag)
        @info "Saving results to $(datadir("sims", savetag))"
    
        accepted = DataFrame(f.samples', f.prior.labels)
        accepted[!,:weight] .= f.weights
        accepted[!,:loss] = vcat(f.losses...)
        
        settings = DataFrame(n = n, q_dist = q_dist, t_max = t_max)

        CSV.write(datadir("sims", savetag, "accepted.csv"), accepted)
        CSV.write(datadir("sims", savetag, "settings.csv"), settings)

        # saving posterior summary to csv + tex  
        posterior_summary = generate_posterior_summary(
            f.prior.labels,
            f.samples, 
            f.losses, 
            f.weights;
            tex = !isnothing(savetag),
            savetag = savetag,
            paramlabels = paramlabels
            )

    end

    return (
        particles = all_particles, 
        weights = all_weights, 
        dists = all_losses, 
        vars = all_vars,
        posterior_summary = posterior_summary
    )
end

function posterior_predictions(f::ModelFit, n::Int64 = 100)

    predictions = Vector{Any}(undef, n)
    samples = Vector{Any}(undef, n)

    @showprogress @threads for i in 1:n        
        sample = posterior_sample(f)
        predictions[i] = f.simulator(sample)
        samples[i] = sample
    end

    return (predictions=predictions, samples=samples)
end

#### %%%% Functions for frequentist optimization %%%% ####

"""
    define_objective_function(
        data::OrderedDict, 
        simulator::Function,
        lossfun::Function
    )::Function
  

Define an objective function based on `data`, `simulator`, `lossfun`. <br>
`simulator` and `lossfun` are assumed to be generated by `define_fitting_simulator` 
and `define_loss_function`, respectively.
"""
function define_objective_function(
    data::OrderedDict, 
    simulator::Function,
    lossfun::Function
    )::Function


    function objective_function(pvec; kwargs...)::Float64

        sim = simulator(pvec; kwargs...) 
        l = lossfun(data, sim)

        return l
    end

    return objective_function
end

function define_objective_function(f::ModelFit)

    return define_objective_function(
        f.data, 
        f.simulator,
        f.loss 
    )

end

function load_pmchist_from_checkpoint(path::String)
    chk = load(path)

    return (
        particles = chk["particles"],
        weights = chk["weights"], 
        dists = chk["losses"], 
        vars = chk["variances"]

    )
end

function load_pmchist_from_checkpoint!(f::ModelFit, path::String)::Nothing

    pmchist = load_pmchist_from_checkpoint(path)
    f.samples = pmchist.particles[end]
    f.weights = pmchist.weights[end]
    f.losses = pmchist.dists[end]

    return nothing

end

function assign_value_by_label!(p, label, value)::Nothing

    labels = EcotoxSystems.ComponentArrays.labels(p)
    idx = findfirst(x -> x == label, labels)
    
    p[idx] = value

    return nothing

end


function assign_values_from_file!(p, file; exceptions = OrderedDict())::Nothing

    posterior_summary = CSV.read(file, DataFrame)

    for (label,value) in zip(posterior_summary.param, posterior_summary.best_fit)
        if !(label in keys(exceptions))
            assign_value_by_label!(p, label, value)
        else
            exceptions[label](p, label, value)
        end
    end

    return nothing

end

function show_intguess(f)
    p0 = mode.(f.prior.dists)
    sim0 = [f.simulator(p0) for _ in 1:10]
    plt = f.plot_data(); plot_sims!(plt, sim0)
    return plt
end