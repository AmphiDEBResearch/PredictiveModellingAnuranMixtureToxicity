
function define_defaultparams_UCLM_mix(;
    posterior_summary_larvalfit=datadir("sims", SAVETAG_LARVALFIT, "posterior_summary.csv"),
    posterior_summary_juvenilefit=datadir("sims", SAVETAG_JUVENILEFIT, "posterior_summary.csv"),
    posterior_summary_24dfit=datadir("sims", SAVETAG_24DFIT, "posterior_summary.csv"),
    posterior_summary_flpfit=datadir("sims", SAVETAG_FLPFIT, "posterior_summary.csv")
    )

    p = ComponentVector(
        glb = ComponentVector(
            t_max = 56.0,
            N0 = 1.0,
            dX_in = [20.0, 20.0],
            k_V = [0.0, 0.0],
            V_patch = [1.0, 1.0],
            T = 293.15,
            C_W = [0.0, 0.0],
            pathogen_inoculation_dose = 0.0,
            pathogen_inoculation_time = 30.0,
            medium_renewals = [0.0]
            ),
        pth = AmphiDEB.defaultparams.pth,
        spc = ComponentVector(
            Z = Dirac(1.0),
            propagate_zoom = (
                dI_max_emb = 0.3333333333333333,
                dI_max_lrv = 0.3333333333333333,
                dI_max_juv = 0.3333333333333333,
                X_emb_int = 1.0,
                H_j1 = 1.0,
                H_p = 1.0,
                K_X_lrv = 0.3333333333333333,
                K_X_juv = 0.3333333333333333
                ),
            X_emb_int = 1,
            K_X_lrv = 1.0,
            K_X_juv = 1.0,
            dI_max_emb = 1,
            dI_max_lrv = 1,
            dI_max_juv = 1,
            kappa_emb = 0.8,
            kappa_juv = 0.8,
            gamma = 0.5,
            eta_IA = 0.54,
            eta_AS_emb = 0.4,
            eta_AS_juv = 0.4,
            eta_AR = 0.95,
            eta_SA = 0.8,
            k_M_emb = 0.11,
            k_M_juv = 0.11,
            delta_k_M_mt = 1.0,
            k_J_emb = 0.027,
            k_J_juv = 0.027,
            H_j1 = 1,
            H_p = 55.0,
            delta_E = 1.0,
            T_A = 8000.0,
            T_ref = 293.15,
            b_T = 40.0,
            fb_G = 0.0,
            h_b = 0.0,
            
            KD = [0. 0. 0. 0. 0. 0. 0.; 0. 0. 0. 0. 0. 0. 0.], # k_D - value per PMoA (G,M,A,R,H,kap) and stressor (1 row = 1 stressor)
            B = [2. 2. 2. 2. 2. 2. 2.; 2. 2. 2. 2. 2. 2. 2.], # slope parameters
            E = [1e10 1e10 1e10 1e10 1e10 1e10 1e10; 1e10 1e10 1e10 1e10 1e10 1e10 1e10], # sensitivity parameters (thresholds)
            KD_h = [0.; 0.], # k_D - value for GUTS-Sd module (1 row = 1 stressor)
            E_h = [1e10; 1e10], # sensitivity parameter (threshold) for GUTS-SD module
            B_h = [1.; 1.], # slope parameter for GUTS-SD module 
            C_h = [1.; 1.], # proportionality constant to convert relative response to hazard rate 
            
            S_rel_crit = 0.66,
            h_S = 0.6,
            a_max = Truncated(Normal(5475.0, 547.5), 0.0, Inf),
            tau_R = 365.0,
            Chi = LogNormal(1.0, 1.0),
            E_P = [Inf, Inf, Inf, Inf],
            B_P = [2.0, 2.0, 2.0, 2.0],
                
            # auxiliary parameters
            watercontent_larvae = 0.93,
            watercontent_juveniles = 0.85,
            time_since_birth = 15.,
            emb_dev_time = 2.
            )
    )

    # setting global parameters

    p.glb.t_max = 100. # setting simulation time conservatively, for cases where metamorphosis is delayed a lot
    p.glb.pathogen_inoculation_time = Inf # no pathogen inoculation
    p.glb.dX_in = [1e10, 1e10] # ad libitum feeding conditions

    p.spc.Z = truncated(Normal(1, 0.1), 0, Inf)
    p.spc.propagate_zoom.H_j1 = 0. # propagation of zoom factor to H_j1 is turned off => we want variability in the transition to metamorphs

    # adding point estimates from calibrations as defaults

    p.spc.KD .= 0.
    p.spc.B .= 2.
    p.spc.E .= 1e10

    assign_values_from_file!(
        p,
        posterior_summary_larvalfit;
        exceptions = OrderedDict(
            "spc.Z" => (p, label, value) -> p.spc.Z = truncated(Normal(1, value), 0, Inf)
        )
    )

    assign_values_from_file!(
        p, 
        posterior_summary_juvenilefit; 
        exceptions = OrderedDict()
        )
    
    assign_values_from_file!(
        p, 
        posterior_summary_24dfit; 
        exceptions = OrderedDict()
        )
    
    
    postsum_flp = CSV.read(posterior_summary_flpfit, DataFrame)
    set_stressor_idx!(postsum_flp, 2) # flp = stressor 2
    CSV.write("postsum_flp.csv", postsum_flp) # dummy file for assigning parameters with correct stressor idx
    assign_values_from_file!(p, "postsum_flp.csv"; exceptions = OrderedDict())
    rm("postsum_flp.csv") # remove dummy file

    p.spc.k_M_juv = p.spc.k_M_emb
    p.spc.X_emb_int = 1. # ≈ initial dry mass of an egg (mg)

    return p
end