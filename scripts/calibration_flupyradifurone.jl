include("boilerplate.jl")

using DataFrames, DataFramesMeta, CSV 
using StatsPlots, Plots.Measures
default(leg = false)
theme(:default)

using EcotoxSystems, AmphiDEB

using Revise

const SAVETAG_LARVALFIT = "input/Discoglossus_larvae" # directory from which larval/metamorph parameters are loaded
const SAVETAG_JUVENILEFIT = "input/Discoglossus_juveniles" # directory from which juvenile/adult parameters are loaded
const SAVETAG = "Discoglossus_Flupyradifurone"

using Revise

include(scriptsdir("utils.jl"));
includet(scriptsdir("Discoglossus_galganoi_Flupyradifurone", "fit.jl")) 

# ======================================== #
# Perior checks
# ======================================== #

f = setup_modelfit("A")

prior_check = prior_predictive_check(f);

plt = f.plot_data()
plot_sims!(plt, prior_check.predictions)
plt

# ======================================== #
# Model fits
# ======================================== #

# ---- PMoA G

i = 1
f = setup_modelfit(PMOAS[i]);

pmcsettings =  (
    n_init = 20_000, 
    n = 10_000, 
    t_max = 10, 
    q_dist = .01
)

@time @suppress pmchist, posterior_check = fit_model!(
    f; 
    pmcsettings = pmcsettings,
    savetag = "$(SAVETAG)_$(PMOAS[i])",
); 

# ---- PMoA κ

i = 7
f = setup_modelfit(PMOAS[i]);

@time @suppress global pmchist, posterior_check = fit_model!(
    f; 
    savetag = "$(SAVETAG)_$(PMOAS[i])",
    pmcsettings
); 
