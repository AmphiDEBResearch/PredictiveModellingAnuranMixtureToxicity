# ========================================================== #
# traits.jl
# functions to infer traits from trajectories and parameters
# ========================================================== #

using DataFrames, DataFramesMeta
using AmphiDEB.ComponentArrays

"""
    calc_S_max_lrv(spc::ComponentVector)

Calculates maximum structural mass of larvae from parameters. 
"""
calc_S_max_lrv(spc::ComponentVector) = ((spc.kappa_emb * spc.dI_max_lrv * spc.eta_IA) / spc.k_M_emb)^3

"""
    calc_E_mt_max(spc::ComponentVector)

Calculates maximum reserve level of larvae from parameters.
"""
calc_E_mt_max(spc::ComponentVector) = (calc_S_max_lrv(spc)*spc.gamma)/(1-spc.gamma)


"""
    calc_Wdry_max_lrv

Calculates maximum dry mass of larvae from parameters.
"""
calc_Wdry_max_lrv(spc::ComponentVector) = calc_S_max_lrv(spc) + calc_E_mt_max(spc)


"""
    metamorphosis_energy_budgets(sim::AbstractDataFrame)

Infer the the energy budget during metamorphic climax from simulation output `sim`. <br>
Note that the relative amounts do not add up to 1 because we also have dissipation (e.g. through growth efficiency).
"""
function metamorphosis_energy_budgets(sim::AbstractDataFrame)
    
    tj1, tj2 = robustextrema(sim[sim.metamorph .> 0.5,:t])

    sim_tj2 = sim[sim.t .== tj2,:][1,:]
    sim_tj1 = sim[sim.t .== tj1,:][1,:]

    total_energy_spent = (sim_tj2.A - sim_tj1.A) + (sim_tj1.E_mt - sim_tj2.E_mt)

    dS = sim_tj2.S[1] - sim_tj1.S[1]
    dM = sim_tj2.M[1] - sim_tj1.M[1]
    dH = sim_tj2.H[1] - sim_tj1.H[1]
    dJ = sim_tj2.J[1] - sim_tj1.J[1]

    dS_rel = dS ./ total_energy_spent[1]
    dM_rel = dM ./ total_energy_spent[1]
    dH_rel = dH ./ total_energy_spent[1]
    dJ_rel = dJ ./ total_energy_spent[1]

    return (dS_rel = dS_rel, dM_rel = dM_rel, dH_rel = dH_rel, dJ_rel = dJ_rel)

end

function robustextrema(v::AbstractVector)

    if length(v)==0
        return Inf,Inf
    end

    return extrema(v)
end

"""
Infer timing of metamorphosis from simulation output `sim`. 
"""
metamorphosis_timing(sim::AbstractDataFrame)::Tuple{Real,Real} = robustextrema(sim[sim.metamorph .> 0.5,:t])

"""
Infer duration of metamorphosis from simulation output `sim`.
"""
metamorphosis_duration(sim::AbstractDataFrame) = diff([metamorphosis_timing(sim)...])[1]

"""
Infer duration of larval development from simulation output `sim`.
"""
larval_period_duration(sim::AbstractDataFrame) = diff([robustextrema(sim[sim.larva .> 0.5,:t])...])[1]


"""
    metamorphosis_weightchange(sim::AbstractDataFrame)::Real

Compute the relative weight change during metamorphic climax. 
"""
function metamorphosis_weightchange(sim::AbstractDataFrame)::Real
    tj1, tj2 = metamorphosis_timing(sim)

    sim_tj2 = sim[sim.t .== tj2,:]
    sim_tj1 = sim[sim.t .== tj1,:]

    W_tot2 = nrow(sim_tj2) > 0 ? sim_tj2.S[1] + sim_tj2.E_mt[1] : NaN
    W_tot1 = nrow(sim_tj1) > 0 ? sim_tj1.S[1] + sim_tj1.E_mt[1] : NaN

    return (W_tot2 - W_tot1)/W_tot1
end

"""
    metamorphosis_reserve_fraction(sim::AbstractDataFrame)::Real

Compute the fraction of total body mass which is reserves at the beginning of metamorphosis. 
"""
function metamorphosis_reserve_fraction(sim::AbstractDataFrame)::Real

    sim[sim.t .== metamorphosis_timing(sim)[1],:] |>
    x -> nrow(x)>0 ? x.E_mt[1]/(x.S[1] + x.E_mt[1]) : NaN
end

"""
    age_at_birth(sim::AbstractDataFrame; rule = "X_emb")

Compute age at birth based on the trajectory of the vitellus `X_emb` (`rule="X_emb"`), or as the first time-point in the series (`rule="init"`, for output that has been processed to omit embryonic period).
"""
function age_at_birth(sim::AbstractDataFrame; rule = "X_emb")

    # infers age at birth based on DEBkiss rule (when egg buffer is used up)
    if rule == "X_emb" 
        return sim[sim.X_emb .<= 0,:].t |>  EcotoxSystems.robustmin
    end

    # infers age at birth assuming that time-series starts at birth
    if rule == "init"
        return sim.t |> EcotoxSystems.robustminimum
    end

    error("Unkown rule to compute age at birth: $rule")

end

"""
    estimate_emb_dev_time(p::AmphiDEB.ComponentVector)


Estimate embryonic development time from parameters. Internally calls `age_at_birth()`.
"""
function estimate_emb_dev_time(p::AmphiDEB.ComponentVector; kwargs...)

    emb_dev_time = AmphiDEB.ODE_simulator(p; kwargs...) |> age_at_birth
    
    return emb_dev_time
end

"""
    drymass_at_birth(sim::AbstractDataFrame; rule = "X_emb")

Compute dry mass at birth. `rule` propagates to `age_at_birth`
"""
function drymass_at_birth(sim::AbstractDataFrame; rule = "X_emb")

    age_at_birth = age_at_birth(sim, rule = rule)

    if isfinite(age_at_birth)
        return sim[sim.t .== age_at_birth,:].W[1] 
    else 
        return NaN
    end
end

function drymass_at_G42(sim::AbstractDataFrame)

    t_G42, _ = metamorphosis_timing(sim)
    
    if isfinite(t_G42)
        df_G42 = sim[sim.t .== t_G42,:]
        return df_G42.S[1] + df_G42.E_mt[1]
    else
        return Inf
    end
end

# NOTE: this assumes that wet mass has been inferred before, based on known or estimated water content for all life stages
function wetmass_at_G42(sim::AbstractDataFrame)

    t_G42, _ = metamorphosis_timing(sim)
    
    if isfinite(t_G42)
        df_G42 = sim[sim.t .== t_G42,:]
        return df_G42.wetmass_mg[1]
    else
        return Inf
    end
end

robustmax(x::AbstractVector) = length(x)>0 ? maximum(x) : Inf
robustmin(x::AbstractVector) = length(x)>0 ? minimum(x) : Inf
robustmean(x::AbstractVector) = length(x)>0 ? mean(x) : Inf
robustmedian(x::AbstractVector) = length(x)>0 ? median(x) : Inf
robustiqr(x::AbstractVector) = length(x)>0 ? iqr(x) : Inf

function drymass_at_G46(sim::AbstractDataFrame)

    _, t_G46 = metamorphosis_timing(sim)

    if isfinite(t_G46)
        df_G46 = sim[sim.t .== t_G46,:]
        return df_G46.S[1] + df_G46.E_mt[1]
    else
        return Inf
    end
end


# NOTE: this assumes that wet mass has been inferred before, based on known or estimated water content for all life stages
function wetmass_at_G46(sim::AbstractDataFrame)

    _, t_G46 = metamorphosis_timing(sim)

    if isfinite(t_G46)
        df_G46 = sim[sim.t .== t_G46,:]
        return df_G46.wetmass_mg[1]
    else
        return Inf
    end
end

function larval_growth_rate(sim::AbstractDataFrame)
    let lrv = sim[sim.larva .> 0.5,:] 

        dt = larval_period_duration(sim) 
        W0, W1 = lrv.S[1]+lrv.E_mt[1], lrv.S[end]+lrv.E_mt[end]

        return (log(W1) - log(W0))/dt
    end
end


"""
    calc_S_max_juv(spc::ComponentVector)

Calculates maximum structural mass of larvae from parameters. 
"""
calc_S_max_juv(spc::EcotoxSystems.ComponentVector) = ((spc.kappa_juv * spc.dI_max_juv * spc.eta_IA) / spc.k_M_juv)^3

function calc_wetmass(r::DataFrameRow, watercontent_larvae, watercontent_juveniles)

    try 
        let lifestage = ["embryo", "larva", "metamorph", "juvenile", "adult"][argmax([r.embryo, r.larva, r.metamorph, r.juvenile, r.adult])]
            # assuming same water content for larvae and embryos 
            # (irrelevant for this study, but we need to enter some value for embryos)
            if (lifestage=="embryo") || (lifestage=="larva") #(r.larva > 0.5) || (r.embryo > 0.5) 
                return r.drymass_mg/(1 - watercontent_larvae)
            # water content for juveniles and adults
            elseif (lifestage=="juvenile") || (lifestage=="adult") #(r.juvenile > 0.5) || (r.adult > 0.5)
                return r.drymass_mg/(1 - watercontent_juveniles)
            # using linear interpolation based on reserve level to calculate water content of metamorphs 
            # the interpolation is irrelevant when we have fixed-stage measurements for G42 and G46 and none in between (as is the case in this study) 
            # what matters is that we have watercontent_larvae at the end of larval period 
            # and switch to watercontent_juveniles after metamorphosis
            elseif lifestage=="metamorph" #r.metamorph > 0.5
                weight_larva = r.E_mt / r.E_mt_max
                weight_juvenile = 1 - weight_larva
                watercontent_metamorph = weight_larva*watercontent_larvae + weight_juvenile*watercontent_juveniles
                return r.drymass_mg/(1 - watercontent_metamorph)
            end
            error("Did not match any life stages")
        end
    catch e
        println((embryo = r.embryo, larva = r.larva, metamorph = r.metamorph, juvenile = r.juvenile, adult = r.adult))
        global dfrow = r 
        error(e)
    end
end

"""
    clutchsize_at_specific_wetmass(sim::AbstractDataFrame)

Infer clutch size produced by female at given body size, expressed in wet mass. 

Currently assumes that the global constants `WETMASS_AT_REPRO_MEASURE_MG`, `REPRODUCTION_PERIOD` and `EGG_DRYMASS_MG` are defined. 

The function extracts the simulation output row for which simulated wetmass is closest to `WETMASS_AT_REPRO_MEASURE_MG`, 
then back-tracks the simulation output so that the time-difference between both output rows matches `REPRODUCTION PERIOD`. 
If the second time-point falls into a pre-adult life stage, it is replaced with time at maturity. 
The returned clutch size is then calculated from the change in reproduction buffer over the extracted period, 
divided by `EGG_DRYMASS_MG`.
"""
function clutchsize_at_specific_wetmass(sim::AbstractDataFrame)
    let idx1

        # extract time at maturity 
        t_p = EcotoxSystems.robustmin(@subset(sim, :adult .> 0.5).t)

        # if maturity has never been reached, we can stop here and return NaN
        if isinf(t_p)
            return NaN
        end

        # extract output row at time-point where simulated wetmass is closest to the specified wetmass
        
        idx2 = argmin(abs.(sim.wetmass_mg .- WETMASS_AT_REPRO_MEASURE_MG))
        df2 = sim[idx2,[:t,:R]]
        t2 = df2.t[1]

        #@assert nrow(df2)==1 "df2 has too many rows - missing grouping variable?"
        
        # extrat output row one year before that (or whatever the value of REPRODUCTION_PERIOD is)
        
        t1 = max(t2 - REPRODUCTION_PERIOD, t_p)
        idx1 = argmin(abs.(sim.t .- t1))
        df1 = sim[idx1,[:t,:R]]

        #@assert nrow(df1)==1 "df1 has too many rows - missing grouping variable?"

        # convert to clutch size corresponding to the amount of reproduction buffer accumulated in the last year

        return (df2.R[1]-df1.R[1])/EGG_DRYMASS_MG
    end
end

"""
    post_metamorphic_weightchange(sim::AbstractDataFrame; dt = 10)

Infer post-metamorphic relative weight change from simulation output, calculated from weight measurements 
at two time points `t1` and `t2` as `(W2-W1)/dt`.

## Arguments

- `t1_post_metam` is the time difference between completion of metamorphosis (G46) and the first measurement.
- `dt` is the time difference between both measurements.
"""
function post_metamorphic_weightchange(sim::AbstractDataFrame; dt::Real = 10, t1_post_metam::Real = 0)::Float64

    t1 = metamorphosis_timing(sim)[2] + t1_post_metam
    t2 = t1 + dt

    W1 = sim[isapprox.(t1, sim.t),:] |> x -> nrow(x) >0 ? x.drymass_mg[1] : NaN
    W2 = sim[isapprox.(t2, sim.t),:] |> x -> nrow(x) >0 ? x.drymass_mg[1] : NaN

    return (W2-W1)/dt
end


