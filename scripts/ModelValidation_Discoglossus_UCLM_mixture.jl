using DrWatson
@quickactivate "."

using DataFrames, DataFramesMeta, CSV 
using StatsPlots, Plots.Measures
default(leg = false)
theme(:default)

using EcotoxSystems, AmphiDEB
include(srcdir("utils.jl"));

using Revise

const SAVETAG_LARVALFIT = joinpath("input", "Discoglossus_larvae")
const SAVETAG_JUVENILEFIT = joinpath("input", "Discoglossus_juveniles")
const SAVETAG_24DFIT = joinpath("Discoglossus_24D_2024-06-23_numtadpoles", "Discoglossus_24D_M")
const SAVETAG_FLPFIT = joinpath("Discoglossus_Flupyradifurone_2024-06-23_numtadpoles", "Discoglossus_Flupyradifurone_G") 

using Revise

includet(scriptsdir("ModelValidation_Discoglossus_UCLM_mixture.jl")) # code to set up the mixture simulations

using Pkg
Pkg.status("AmphiDEB") # show version of AmphiDEB used

Pkg.status("EcotoxSystems") # show version of EcotoxSystems used

Pkg.status("EcotoxModelFitting")

p = define_defaultparams_UCLM_mix()
sims = [simulator_UCLM_mixture(p) for _ in 1:100];
save_sims(sims, datadir("sims", "ModelValidation_Discoglossus_UCLM_mixture"), "")

plt = plot_data_UCLM_mix_growth()
plot_sims_UCLM_mix_growth!(plt, sims)

savefig(plot(plt, dpi = 400), plotsdir("ModelValidation_Discoglossus_UCLM_mixture_growth.png"))
display(plt)

plt = plot_data_UCLM_mix_fracttadpoles()
plot_sims_UCLM_mix_fracttadpoles!(plt, sims)

savefig(plot(plt, dpi = 400), plotsdir("ModelValidation_Discoglossus_UCLM_mixture_fracttadpoles.png"))
display(plt)

# extracting simulated metamorph data from simulation vector
sim_metamorphs = extract_simkey(sims, :metamorphs)
leftjoin!(sim_metamorphs, TREATMENT_IDS_MIX, on = :treatment_id);

# reading data
paths = OrderedDict(
        :aquatic => [datadir("exp_raw", "UCLM", "discoglossus_galganoi", "disco_flp_24D_aquatic.csv"), 1], # number indicates row where data header is located (omitting metadata)
        :metamorphs => [datadir("exp_raw", "UCLM", "discoglossus_galganoi", "disco_flp_24D_metamorphs.csv"), 1]
)
data = OrderedDict()

for (key, info) in pairs(paths)
    path, header = info
    data[key] = CSV.read(path, DataFrame, header=header)
end

data[:metamorphs].wetmass_G42_mg = float.(data[:metamorphs].wetmass_G42_mg)
data[:metamorphs].wetmass_G46_mg = float.(data[:metamorphs].wetmass_G46_mg)

leftjoin!(
    data[:metamorphs], TREATMENT_IDS_MIX, on=[:D_ppm, :F_ppm]
);

# plotting data + predictions
plt = plot_data_UCLM_mix_metamorphs()
plot_sims_UCLM_mix_metamorphs!(plt, sims)
savefig(plot(plt, dpi = 400), plotsdir("ModelValidation_Discoglossus_UCLM_mixture_metamorphs.png"))
display(plt)
