MODULE DMFT_ED
  USE ED_INPUT_VARS  , only: &
       ed_read_input , &
       Norb          , &
       Nspin         , &
       Nloop         , &
       Nph           , &
       Uloc          , &
       Ust           , &
       Jh            , &
       Jx            , &
       Jp            , &
       xmu           , &
       beta          , &
       g_ph          , &
       w0_ph         , &
       eps           , &
       wini          , &
       wfin          , &
       xmin          , &
       xmax          , &
       Nsuccess      , &
       dmft_error    , &
       sb_field      , &
       cg_Scheme     , &
       nread         , &
       Lmats         , &
       Lreal         , &
       Lpos          , &
       Hfile         , &
       HLOCfile      , &
       LOGfile       , &
       ed_mode       , &
       bath_type 

  USE ED_BATH, only:                                                    &
       ed_set_Hreplica                 => set_Hreplica                 , &
       ed_get_bath_dimension           => get_bath_dimension           , &
       ed_get_bath_component_dimension => get_bath_component_dimension , &
       ed_get_bath_component           => get_bath_component           , &
       ed_set_bath_component           => set_bath_component           , &
       ed_copy_bath_component          => copy_bath_component          , &
       ed_spin_symmetrize_bath         => spin_symmetrize_bath         , &
       ed_orb_symmetrize_bath          => orb_symmetrize_bath          , &
       ed_orb_equality_bath            => orb_equality_bath            , &
       ed_ph_symmetrize_bath           => ph_symmetrize_bath           , &
       ed_ph_trans_bath                => ph_trans_bath                , &
       ed_break_symmetry_bath          => break_symmetry_bath          , &
       ed_enforce_normal_bath          => enforce_normal_bath


  USE ED_AUX_FUNX, only:                        &
       ed_set_suffix                                                     , &
       ed_reset_suffix                                                   , &
       ed_so2os_reshape                 => so2os_reshape                 , &
       ed_os2so_reshape                 => os2so_reshape                 , &
       ed_search_variable                                                , &
       ed_search_chemical_potential     => search_chemical_potential     , &
       ed_atomic_SOC                    => atomic_SOC                    , &
       ed_atomic_SOC_rotation           => atomic_SOC_rotation           , &
       ed_orbital_Lz_rotation_Norb      => orbital_Lz_rotation_Norb      , &
       ed_orbital_Lz_rotation_NorbNspin => orbital_Lz_rotation_NorbNspin , &
       ed_atomic_j                      => atomic_j


  USE ED_IO, only: &
       ed_get_gimp            , &
       ed_get_sigma           , &
       ed_get_sigma_matsubara , &
       ed_get_sigma_realaxis  , &
       ed_get_self_matsubara  , &
       ed_get_self_realaxis   , &
       ed_get_gimp_matsubara  , &
       ed_get_gimp_realaxis   , &       
       ed_get_fimp_matsubara  , &
       ed_get_fimp_realaxis   , &
       ed_get_g0imp_matsubara , &
       ed_get_f0imp_matsubara , &
       ed_get_g0imp_realaxis  , &
       ed_get_f0imp_realaxis  , &
       ed_get_delta_function  , &
       ed_get_fdelta_function , &
       ed_get_dens            , &
       ed_get_mag             , &
       ed_get_docc            , &
       ed_get_eimp            , &
       ed_get_epot            , &
       ed_get_eint            , &
       ed_get_ehartree        , &
       ed_get_eknot           , &
       ed_get_doubles         , &
       ed_get_density_matrix  , &
       ed_get_quantum_SOC_operators_single, ed_get_quantum_SOC_operators_lattice


  USE ED_MAIN, only:    &
       ed_init_solver , &
       ed_rebuild_gf  , &
       ed_solve

  USE ED_BATH_FIT,  only: ed_chi2_fitgf


END MODULE DMFT_ED

