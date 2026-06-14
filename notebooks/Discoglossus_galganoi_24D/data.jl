
function load_data(;
    paths::OrderedDict = OrderedDict(
        :aquatic => [datadir("exp_raw", "UCLM", "discoglossus_galganoi", "disco_flp_24D_aquatic.csv"), 1], # number indicates row where data header is located (omitting metadata)
        :metamorphs => [datadir("exp_raw", "UCLM", "discoglossus_galganoi", "disco_flp_24D_metamorphs.csv"), 1],
        ##:adults => [datadir("exp_raw", "Discoglossus_03_adults.csv"), 1]
    ))

    data = OrderedDict()

    for (key,info) in pairs(paths)
        path, header = info
        data[key] = CSV.read(path, DataFrame, header = header)
    end

    # process aquatic data

    dropmissing!(data[:aquatic])

    data[:aquatic] = @subset(
        data[:aquatic], 
        :F_ppm .== 0, # omit Flupyradifurone treatments
        :D_ppm .< maximum(:D_ppm) # omit the highest treatment due to 100% mortality
        ) 

    data[:aquatic].num_tadpoles = float.(data[:aquatic].num_tadpoles)
    data[:aquatic].fract_tadpoles = data[:aquatic].num_tadpoles ./ data[:aquatic].survival

    rename!(data[:aquatic], :D_ppm => :C_W_1) 
    data[:aquatic] = EcotoxSystems.relative_response(
        data[:aquatic], 
        [:wetmass_mg, :num_tadpoles], 
        :C_W_1; 
        groupby_vars = [:t_exp, :aquarium]
    )

    # process metamorph data

    data[:metamorphs].wetmass_G42_mg = float.(data[:metamorphs].wetmass_G42_mg)
    data[:metamorphs].wetmass_G46_mg = float.(data[:metamorphs].wetmass_G46_mg)

    data[:metamorphs] = @subset(
        data[:metamorphs], 
        :F_ppm .== 0, # omit Flupyradifurone treatments
        :D_ppm .< maximum(:D_ppm) # omit the highest treatment due to 100% mortality
        ) 
    data[:metamorphs][!,:treatment_id] = [TREATMENT_IDS[t] for t in data[:metamorphs].D_ppm]
    rename!(data[:metamorphs], :D_ppm => :C_W_1)
    dropmissing!(data[:metamorphs])

    return data
end

function plot_data(;kwargs...)

    plt_aqua = @df f.data[:aquatic] plot(
        groupedlineplot(
            :t_exp, :wetmass_mg, :C_W_1, 
            layout = (1,length(unique(:C_W_1))), 
            leg = [:bottomright false false false false], 
            label = "Observed mean (P₅-P₉₅)",
            title = hcat(["$x mg/L" for x in unique(:C_W_1)]...), 
            xlim = (-5, 30), 
            ylim = (100,350), 
            fillalpha = .2, lw = 2, 
            color = :black, marker = true, 
            ylabel = ["Wet mass (mg)" "" "" "" ""],
        ), 
        groupedlineplot(
            :t_exp, :fract_tadpoles, :C_W_1, 
            layout = (1,length(unique(:C_W_1))), 
            leg = false, 
            fillalpha = .2, lw = 2, color = :black, 
            marker = true, 
            xlabel = "Time (d)",  
            ylabel = ["Fraction of \n tadpoles" "" "" "" ""],
            xlim = (-5, 30), 
        ), 
        layout = (2,1), size = (1250,500), leftmargin = 7.5mm, bottommargin = 5mm
    )

    c = 0
    num_concs = length(unique(f.data[:aquatic].C_W_1))

    for (i,C_W) in enumerate(unique(f.data[:aquatic].C_W_1))
        c += 1
        df = @subset(f.data[:aquatic], :C_W_1 .== 0)
        @df df lineplot!(plt_aqua, :t_exp, :wetmass_mg, (.5, .5), subplot = c, color = :gray, lw = 2, linestyle = :dash, label = "")
        @df df lineplot!(plt_aqua, :t_exp, :fract_tadpoles, (.5, .5), subplot = c+num_concs, color = :gray, lw = 2, linestyle = :dash)
    end

    plt = plot(
        plt_aqua;
        kwargs...
        )

    return plt
end

function plot_metamorphs(;kwargs...)

    @df f.data[:metamorphs] plot(
        boxplot(string.(:C_W_1), :t_exp_G42),
        boxplot(string.(:C_W_1), :t_exp_G42),
        boxplot(string.(:C_W_1), :wetmass_G42_mg),
        boxplot(string.(:C_W_1), :wetmass_G46_mg)
    )

    @df f.data[:metamorphs] plot(
        boxplot(string.(:C_W_1), :t_exp_G42, ylim = (0,20)),
        boxplot(string.(:C_W_1), :t_exp_G42, ylim = (0,20)),
        boxplot(string.(:C_W_1), :wetmass_G42_mg, ylim = (50,400)),
        boxplot(string.(:C_W_1), :wetmass_G46_mg, ylim = (50,400)), 
        layout = (1,4), size = (1200,350),
        leg = [true false false false], label = "Observed", 
        xlabel = "2,4D (mg/L)", leftmargin = 7.5mm, bottommargin = 5mm, 
        fillcolor = :gray, markercolor = :black, fillalpha = .5,
        ylabel = ["Time since \n start of experiment (d)" "Time since \n start of experiment (d)" "Wet mass (mg)" "Wet mass (mg)"], 
        title = ["Timing of Gosner 42" "Timing of Gosner 46" "Mass at Gosner 42" "Mass at Gonser 46"], titlefontsize = 10, labelfontsize = 10, 
    )


end

function plot_sims!(plt, predictions::AbstractVector; label = "Simulation", color = :steelblue)::Nothing

    df_aqua = sort(vcat([df[:aquatic] for df in predictions]...), :t_exp)

    c = 0
    num_concs = length(unique(df_aqua.C_W_1))

    for (i,C_W) in enumerate(unique(df_aqua.C_W_1))
        c += 1
        df = @subset(df_aqua, :C_W_1 .== C_W)
        @df df lineplot!(plt, :t_exp, :wetmass_mg, subplot = c, color = color, lw = 2, fillalpha = .2, label = label)
        @df df lineplot!(plt, :t_exp, :fract_tadpoles, subplot = c+num_concs, color = color, lw = 2, fillalpha = .2)
    end

    return nothing
end