using DrWatson
@quickactivate "."
#using Pkg; Pkg.instantiate()

using DataFrames, DataFramesMeta, CSV 
using StatsPlots, Plots.Measures
default(leg = false)
theme(:default)

using EcotoxSystems, AmphiDEB
include(scriptsdir("utils.jl"));

using Revise

const SAVETAG_LARVALFIT = "input/Discoglossus_larvae" # directory from which larval/metamorph parameters are loaded
const SAVETAG_JUVENILEFIT = "input/Discoglossus_juveniles" # directory from which juvenile/adult parameters are loaded
const SAVETAG = "Discoglossus_Flupyradifurone"

using Revise

includet(scriptsdir("Discoglossus_galganoi_Flupyradifurone", "fit.jl")) 

f = setup_modelfit("A")

prior_check = prior_predictive_check(f);

plt = f.plot_data()
plot_predictions!(plt, prior_check.predictions)
plt

i = 1
f = setup_modelfit(PMOAS[i]);

@time @suppress global pmchist, posterior_check = fit_model!(
    f; 
    savetag = "$(SAVETAG)_$(PMOAS[i])",
    n_init = 20_000, 
    n = 10_000, 
    t_max = 10, 
    q_dist = .01,
); 

i = 7
f = setup_modelfit(PMOAS[i]);

@time @suppress global pmchist, posterior_check = fit_model!(
    f; 
    savetag = "$(SAVETAG)_$(PMOAS[i])",
    n_init = 20_000, 
    n = 10_000, 
    t_max = 10, 
    q_dist = .01,
); 
