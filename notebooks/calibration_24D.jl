using Pkg; Pkg.activate(".")
include("boilerplate.jl")

using DataFrames, DataFramesMeta, CSV 
using StatsPlots, Plots.Measures
default(leg = false)
theme(:default)

using EcotoxSystems, AmphiDEB
include(srcdir("utils.jl"));

using Revise
# directory from which larval/metamorph parameters are loaded
const SAVETAG_LARVALFIT = "input/Discoglossus_larvae" 
 # directory from which juvenile/adult parameters are loaded
const SAVETAG_JUVENILEFIT = "input/Discoglossus_juveniles"
const SAVETAG = "Discoglossus_24D"
using Revise

include(notebooksdir("Discoglossus_galganoi_24D", "fit.jl")) 

# ======================================== #
# Prior checks per PMoA
# ======================================== #

f = setup_modelfit("G", sigma_factor = 2) 

plt = f.plot_data()
prior_check = prior_predictive_check(f, n = 100)
plot_sims!(plt, prior_check.predictions)
plt

f = setup_modelfit("M", sigma_factor = 2) 

plt = f.plot_data()
prior_check = prior_predictive_check(f, n = 100)
plot_sims!(plt, prior_check.predictions)
plt

f = setup_modelfit("A", sigma_factor = 2) 

plt = f.plot_data()
prior_check = prior_predictive_check(f, n = 100)
plot_sims!(plt, prior_check.predictions)

plt

# ======================================== #
# Model fits per PMoA
# ======================================== #

i = 1 
f = setup_modelfit(PMOAS[i]);

 pmcsettings = (
    n_init = 2_000,#20_000, 
    n = 1_000,#10_000, 
    t_max = 5,#10, 
    q_dist = .1,#.01,
    )


pmchist = run_PMC!(
    f; 
    pmcsettings...
)

@time @suppress global pmchist, posterior_check = fit_model!(
    f; 
    pmcsettings = pmcsettings,
    savetag = "$(SAVETAG)_$(PMOAS[i])",
); 

i = 2 
f = setup_modelfit(PMOAS[i]; loss_functions = loss_mse);

@time @suppress global pmchist, posterior_check = fit_model!(
    f; 
    savetag = "$(SAVETAG)_$(PMOAS[i])",
    n_init = 20_000, 
    n = 10_000, 
    t_max = 10, 
    q_dist = .01,
);

i = 3
f = setup_modelfit(PMOAS[i]; loss_functions = loss_mse);

@time @suppress global pmchist, posterior_check = fit_model!(
    f; 
    savetag = "$(SAVETAG)_$(PMOAS[i])",
    n_init = 20_000, 
    n = 10_000, 
    t_max = 10, 
    q_dist = .01,
);
