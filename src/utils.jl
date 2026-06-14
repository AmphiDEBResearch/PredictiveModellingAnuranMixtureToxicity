using DataFrames 
using Latexify
using DataStructures
import EcotoxModelFitting: Hyperdist

function getjoinvars(f::ModelFit, i::Int64)

    if !(f.time_resolved[i])
        return f.grouping_vars[i]
    else
        return vcat(f.time_vars, f.grouping_vars[i])
    end

end

function plot_prior_opt(f, prior_check; plot_sims! = plot_sims!)
    prior_opt = get_prior_opt(prior_check)
    sims = [f.simulator(prior_opt) for _ in 1:10]
    plt = f.plot_data()
    plot_sims!(plt, sims)
    return plt
end

function get_prior_opt(prior_check)
     return prior_check.samples[prior_check.losses .== minimum(filter(isfinite, prior_check.losses))][1]
end


function save_sims(
    predictions::AbstractVector, 
    savedir::AbstractString,
    savetag::AbstractString, 
    prefix::AbstractString
    )::Nothing

    for key in keys(predictions[1])
        df = vcat([@transform(p[key], :num_sample = i) for (i,p) in enumerate(predictions)]...)
        CSV.write(joinpath(savedir, savetag, "$(prefix)_$(key).csv"), df)
    end

    return nothing
end


"""
Replaces the stressor index "1" in parameter labels with `idx`. 
"""
function set_stressor_idx!(
    posterior_summary::AbstractDataFrame, 
    idx::Int64
    )

    posterior_summary.param = [replace(label, "[1,"=>"[$idx,") for label in posterior_summary.param]
    
    return nothing
end

"""
    plot_metam_phase(sim::DataFrame)

Plot state varariables with indication of life stages
"""
function plot_metam_phase(sim::DataFrame)

    plt = @df sim plot(
        plot(:t, :S, color = :black, lw = 1.5, ylabel = "S"), 
        plot(:t, :E_mt, color = :black, lw = 1.5, ylabel = "E_mt"), 
        xlabel = "Time (d)", xlim = (0,100), leg = false
    )

    for (i,var) in enumerate([:S, :E_mt])
        @df sim plot!(subplot = i, :t, :larva .* maximum(sim[:,var]), fill = true, fillalpha = .2, lw = 0, color = :purple, linetype = :stepmid)
        @df sim plot!(subplot = i, :t, :metamorph .* maximum(sim[:,var]), fill = true, fillalpha = .2, lw = 0, color = :steelblue, linetype = :stepmid)
        @df sim plot!(subplot = i, :t, :juvenile .* maximum(sim[:,var]), fill = true, fillalpha = .2, lw = 0, color = :teal, linetype = :stepmid)
    end

    return plt
end

macro h(x)
    quote
        display("text/markdown", @doc $x)
    end    
end

function is_finite_row(
    row::DataFrameRow
    )::Bool

    return sum(.!(check_for_nonfinite.(Vector(row)))) == 0

end

check_for_nonfinite(x::Number)::Bool = isfinite(x)
check_for_nonfinite(x::Any)::Bool = true

"""
    clean(df::AbstractDataFrame)

Removes all rows with any non-finite and missing values from dataframe. 
"""
function clean(df::AbstractDataFrame)
    
    valid_idxs = [is_finite_row(row) for row in eachrow(df)]

    return dropmissing(df[valid_idxs,:])

end

"""
    match_order(a::Vector, b::Vector)

Match the order of elements in `a` to the order of elements in `b`, assuming that `a` is a subset of `b`.
"""
function match_order(a::Vector, b::Vector)
    
    index_map = Dict(elem => idx for (idx, elem) in enumerate(a))

    return [a[index_map[elem]] for elem in b]
end


"""
    fround(x; sigdigits=2)
    
Formatted rounding to significant digits (omitting decimal point when appropriate). 
Returns rounded number as string.

"""
function fround(x; sigdigits=2)
    xround = string(round(x, sigdigits = sigdigits))
    if xround[end-1:end]==".0"
        xround = string(xround[1:end-2])
    end
    return xround
end


function _df_to_tex(df::AbstractDataFrame, fname::AbstractString; colnames::Union{Nothing,Vector{AbstractString}} = nothing)::Nothing

    tex_table = @chain df begin
        !isnothing(colnames) ? rename(_, colnames) : _
        latexify(env = :table, booktabs = true, latex = false, fmt = FancyNumberFormatter(3))   
    end 
    
    open(fname, "w") do f
        write(f, tex_table)
    end

    @info "Writing latex table to $fname"

    return nothing
end


"""
Convert parameter object to table (`DataFrame`).
"""
function _as_table(p::EcotoxSystems.ComponentVector; printtable = true)

    df = DataFrame(
        param = EcotoxSystems.ComponentArrays.labels(p), 
        value = vcat(p...)
    )

    if printtable
        show(df, allrows = true)
    end

    return df
end


function reverse(od::OrderedDict)
  new_od = OrderedDict()
  for (k, v) in od
    new_od[v] = k
  end
  return new_od
end

nrmsd(a,b) = sqrt(sum((a .- b).^2) ./ length(b)) / iqr(b) 
mre(a,b) = mean((a .- b) ./ b)
mrae(a,b) = mean(abs.(a .- b) ./ b)

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