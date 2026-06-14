
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


paramlabels = OrderedDict(
    "spc.eta_AS_emb" => L"\eta_{AS}^{emb}",
    "spc.eta_AS_juv" => L"\eta_{AS}^{juv}",
    "spc.eta_AR" => L"\eta_{AR}",
    "spc.dI_max_lrv" => L"\{ \dot{I} \}_{max}^{lrv}",
    "spc.dI_max_juv" => L"\{ \dot{I} \}_{max}^{juv}",
    "spc.k_M_emb" => L"k_M^{emb}",
    "spc.Z" => L"\sigma_{Z,M}",
    "spc.H_j1_prime" => L"H^{j'}",
    "spc.H_p_prime" => L"H^{p'}",
    "spc.gamma" => L"\gamma",
    "spc.kappa_emb" => L"\kappa",
    "spc.watercontent_larvae" => L"Larval\ water\ content",
    "spc.watercontent_juveniles" => L"Juvenile\ water\ content"
)