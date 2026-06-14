
function define_defaultparams()::AmphiDEB.ComponentVector

    p = AmphiDEB.ComponentVector(
        glb = AmphiDEB.defaultparams.glb, 
        pth = AmphiDEB.defaultparams.pth,
        spc = EcotoxSystems.ComponentVector(
            AmphiDEB.defaultparams.spc; 
            # auxiliary parameters
            log_k_D_G = log(1e-10),
            log_k_D_M = log(1e-10),
            log_k_D_A = log(1e-10),
            log_k_D_KAP = log(1e-10),
            log_e_G = log(1e10),
            log_e_M = log(1e10),
            log_e_A = log(1e10),
            log_e_KAP = log(1e10),
            watercontent_larvae = 0.93, 
            watercontent_juveniles = 0.85,
            time_since_birth = 15.,
            emb_dev_time = 2. 
        ))

    # setting global parameters

    p.glb.t_max = 100. # setting simulation time conservatively, for cases where metamorphosis is delayed a lot
    p.glb.pathogen_inoculation_time = Inf # no pathogen inoculation
    p.glb.dX_in = [1e10, 1e10] # ad libitum feeding conditions

    p.spc.Z = truncated(Normal(1, 0.1), 0, Inf)
    # propagation of zoom factor to H_j1 is turned off => we want variability in the transition to metamorphs
    p.spc.propagate_zoom.H_j1 = 0.

    # adding point estimates as defaults

    posterior_summary_larvalfit = CSV.read(datadir("sims", SAVETAG_LARVALFIT, "posterior_summary.csv"), DataFrame)
    posterior_summary_juvenilefit = CSV.read(datadir("sims", SAVETAG_JUVENILEFIT, "posterior_summary.csv"), DataFrame)

    @info "Overwriting default values of $(posterior_summary_larvalfit.param)"
    for (label,value) in zip(posterior_summary_larvalfit.param, posterior_summary_larvalfit.best_fit)
        if label == "spc.Z"
            p.spc.Z = truncated(Normal(1, value), 0, Inf)
        else
            assign_value_by_label!(p, label, value)
        end
    end

    @info "Overwriting default values of $(posterior_summary_juvenilefit.param)"
    for (label,value) in zip(posterior_summary_juvenilefit.param, posterior_summary_juvenilefit.best_fit)
        if label == "spc.Z"
            p.spc.Z = truncated(Normal(1, value), 0, Inf)
        else
            assign_value_by_label!(p, label, value)
        end
    end

    p.spc.k_M_juv = p.spc.k_M_emb
    p.spc.X_emb_int = 1. # ≈ initial dry mass of an egg (mg)

    p.spc.KD .= 0.
    p.spc.B .= 2.
    p.spc.E .= 1e10

    p.spc.emb_dev_time = estimate_emb_dev_time(p)

    return p
end


# misc

paramlabels = OrderedDict(
    "spc.log_k_D_G" => L"ln k_{D,G}",
    "spc.KD[1,1]" => L"k_{D,G}",
    "spc.log_e_G" => L"ln e_G",
    "spc.B[1,1]" => L"b_G",
    
    "spc.log_k_D_M" => L"ln k_{D,M}",
    "spc.KD[1,2]" => L"k_{D,M}",
    "spc.log_e_M" => L"ln e_M",
    "spc.B[1,2]" => L"b_M",
    
    "spc.log_k_D_A" => L"ln k_{D,A}",
    "spc.KD[1,3]" => L"k_{D,A}",
    "spc.log_e_A" => L"ln(e_A)",
    "spc.B[1,3]" => L"b_A",

    "spc.log_k_D_R" => L"ln k_{D,R}",
    "spc.KD[1,4]" => L"\k_{D,R}",
    "spc.log_e_R" => L"ln(e_R)",
    "spc.B[1,4]" => L"b_R",

    "spc.log_k_D_Hneg" => L"ln k_{D,H^-}",
    "spc.KD[1,5]" => L"\k_{D,H^-}",
    "spc.log_e_KAP" => L"ln(e_{H^-})",
    "spc.B[1,5]" => L"b_{H^-}",

    "spc.log_k_D_Hpos" => L"ln k_{D,H^+}",
    "spc.KD[1,6]" => L"\k_{D,H^+}",
    "spc.log_e_KAP" => L"ln(e_{H^+})",
    "spc.B[1,6]" => L"b_{H^+}",

    "spc.log_k_D_KAP" => L"ln k_{D,\kappa}",
    "spc.KD[1,7]" => L"\k_{D,\kappa}",
    "spc.log_e_KAP" => L"ln(e_{\kappa})",
    "spc.B[1,7]" => L"b_{\kappa}",

    "spc.E[1,1]" => L"e_G",
    "spc.E[1,2]" => L"e_M",
    "spc.E[1,3]" => L"e_A",
    "spc.E[1,5]" => L"e_{H^{-}}",
    "spc.E[1,6]" => L"e_{H^{+}}",
    "spc.E[1,7]" => L"e_{\kappa^{-}}",

    "spc.Z_mean_UGent" => L"\overline{Z}_{corr}", 
    "spc.eta_AS_juv" => L"\eta_{AS}^{juv}"
)
