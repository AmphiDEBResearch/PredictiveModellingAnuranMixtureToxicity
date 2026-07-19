include("boilerplate.jl")

using Revise

const SAVETAG_LARVALFIT = "input/Discoglossus_larvae" # directory from which larval/metamorph parameters are loaded
const SAVETAG_JUVENILEFIT = "input/Discoglossus_juveniles" # directory from which juvenile/adult parameters are loaded
const SAVETAG_TKTD = "Discoglossus_Flupyradifurone"
const SAVETAG = "Discoglossus_Flupyradifurone"

using Revise

includet(scriptsdir("Discoglossus_galganoi_Flupyradifurone", "fit.jl"))
includet(scriptsdir("Discoglossus_galganoi_Flupyradifurone", "cross_validation.jl")) 

savedir = joinpath("Discoglossus_Flupyradifurone_2025-06-23_numtadpoles")

acc_G = CSV.read(datadir("sims", savedir, "$(SAVETAG_TKTD)_G", "accepted.csv"), DataFrame) |> x->@transform(x, :pmoa = "G")
acc_kap = CSV.read(datadir("sims", savedir, "$(SAVETAG_TKTD)_KAP", "accepted.csv"), DataFrame) |> x->@transform(x, :pmoa = "KAP")



using ProgressMeter

@showprogress let sims = [], f = setup_modelfit("G")
    for _ in 1:100
        p = posterior_sample(
            acc_G, 
            reserved_colnames = vcat(EcotoxModelFitting.RESERVED_COLNAMES, "pmoa")
            )

        sim_p = f.simulator(p_opt_G)
        push!(sims, sim_p)
    end
end


combined = vcat(
    acc_G[:,[:pmoa, :loss]], 
    acc_kap[:,[:pmoa, :loss]]
)



@df @subset(combined, :pmoa .== "G") histogram(:loss, normalize = :pdf)
@df @subset(combined, :pmoa .== "KAP") histogram!(:loss, normalize = :pdf)



@df histogram(acc_G, normalize = :pdf, weight = :weight)

acc_G
