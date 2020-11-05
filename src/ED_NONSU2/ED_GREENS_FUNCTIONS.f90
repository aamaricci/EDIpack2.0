MODULE ED_GREENS_FUNCTIONS
  USE ED_GF_SHARED
  USE ED_GF_NORMAL
  USE ED_GF_SUPERC
  USE ED_GF_NONSU2
  USE ED_GF_CHISPIN
  USE ED_GF_CHIDENS
  USE ED_GF_CHIPAIR
  !
  implicit none
  private 

  public :: buildGf_impurity
  public :: buildChi_impurity

contains


  !+------------------------------------------------------------------+
  ! GF CALCULATIONS
  !+------------------------------------------------------------------+
  subroutine buildGF_impurity()
    !
    call allocate_grids()
    !
    impGmats=zero
    impGreal=zero
    impFmats=zero
    impFreal=zero
    !
    impSmats = zero
    impSreal = zero
    impSAmats = zero
    impSAreal = zero
    !
    impG0mats=zero
    impG0real=zero
    impF0mats=zero
    impF0real=zero
    !
    !
    write(LOGfile,"(A)")"Get impurity Greens functions:"
    select case(ed_mode)
    case default
       call build_gf_normal()
       call build_sigma_normal()
    case ("superc")
       call build_gf_superc()
       call build_sigma_superc()
    case ("nonsu2")
       call build_gf_nonsu2()
       call build_sigma_nonsu2()
    end select
    !
    if(MPIMASTER)then
       if(ed_print_Sigma)call ed_print_impSigma()
       if(ed_print_G)call ed_print_impG()
       if(ed_print_G0)call ed_print_impG0()
    endif
    !
    call deallocate_grids()
    !
  end subroutine buildGF_impurity





  !+------------------------------------------------------------------+
  ! SUSCEPTIBILITY CALCULATIONS
  !+------------------------------------------------------------------+
  subroutine buildChi_impurity()
    !
    call allocate_grids
    !
    !
    !BUILD SPIN SUSCEPTIBILITY    
    spinChi_tau=zero
    spinChi_w=zero
    spinChi_iv=zero
    if(chispin_flag)call build_chi_spin()
    !
    !BUILD CHARGE SUSCEPTIBILITY
    densChi_tau=zero
    densChi_w=zero
    densChi_iv=zero
    densChi_mix_tau=zero
    densChi_mix_w=zero
    densChi_mix_iv=zero
    densChi_tot_tau=zero
    densChi_tot_w=zero
    densChi_tot_iv=zero
    if(chidens_flag)call build_chi_dens()
    !
    !BUILD PAIR SUSCEPTIBILITY
    pairChi_tau=zero
    pairChi_w=zero
    pairChi_iv=zero
    if(chipair_flag)call build_chi_pair()
    !
    !PRINTING:
    if(MPIMASTER.and.any([chispin_flag,chidens_flag,chipair_flag]))call ed_print_impChi()
    !
    !
    call deallocate_grids
  end subroutine buildChi_impurity




end MODULE ED_GREENS_FUNCTIONS
