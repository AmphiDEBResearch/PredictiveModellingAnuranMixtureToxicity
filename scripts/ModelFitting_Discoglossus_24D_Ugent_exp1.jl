# ModelFitting_Discoglossus_24D_Ugent_exp1.jl
# estimating 2,4-D parameters from effects on juvenile mass after larval pre-expousre (UGent exp 1 data)
# this is DISCONTINUED because effects are veryt small

using OrdinaryDiffEq
import AmphiDEB: ComponentVector

includet(srcdir("ModelFitting.jl"))
includet(srcdir("utils.jl"))
includet(srcdir("loss.jl"))

function load_data_exp1(;
    paths::OrderedDict = OrderedDict(
        :juveniles => [datadir("exp_raw", "UGent", "exp1", "juveniles.csv"), 1]
    ))

    data = OrderedDict()

    for (key,info) in pairs(paths)
        path, header = info
        data[key] = CSV.read(path, DataFrame, header = header) |> 
        x -> @subset(x, :treatment_bd .== "uninfected") 
    end

    # add column for control-normalized weight change
    #data[:juveniles][!,:y_post_metam_weightchange] = data[:juveniles].post_metam_weightchange ./ mean(skipmissing(@subset(data[:juveniles], :pretreatment_24D .== 0).post_metam_weightchange))

    #treatment_ids = Dict(zip(
    #    sort(unique(data[:juveniles].pretreatment_24D)),
    #    collect(eachindex(sort(unique(data[:juveniles].pretreatment_24D))))
    #))
    #
    #data[:juveniles][!,:treatment_id] = [treatment_ids[x] for x in data[:juveniles].pretreatment_24D]

    return data

end
