MODULE ED_IO
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_SETUP
  USE ED_BATH
  USE ED_BATH_FUNCTIONS
  !
  USE SF_LINALG
  USE SF_ARRAYS,  only: linspace,arange
  USE SF_IOTOOLS, only: str,reg,free_unit,splot,sread
  USE SF_MISC,    only: assert_shape
  implicit none
  private

  !Retrieve self-energy through routines:
  interface ed_get_sigma_matsubara
     module procedure ed_get_sigma_matsubara_main
     module procedure ed_get_sigma_matsubara_site
     module procedure ed_get_sigma_matsubara_lattice
  end interface ed_get_sigma_matsubara

  interface ed_get_self_matsubara
     module procedure ed_get_self_matsubara_main
     module procedure ed_get_self_matsubara_site
     module procedure ed_get_self_matsubara_lattice
  end interface ed_get_self_matsubara

  interface ed_get_sigma_realaxis
     module procedure ed_get_sigma_realaxis_main
     module procedure ed_get_sigma_realaxis_site
     module procedure ed_get_sigma_realaxis_lattice
  end interface ed_get_sigma_realaxis

  interface ed_get_self_realaxis
     module procedure ed_get_self_realaxis_main
     module procedure ed_get_self_realaxis_site
     module procedure ed_get_self_realaxis_lattice
  end interface ed_get_self_realaxis


  !Retrieve imp GF through routines.
  interface ed_get_gimp_matsubara
     module procedure ed_get_gimp_matsubara_main
     module procedure ed_get_gimp_matsubara_site
     module procedure ed_get_gimp_matsubara_lattice
  end interface ed_get_gimp_matsubara

  interface ed_get_fimp_matsubara
     module procedure ed_get_fimp_matsubara_main
     module procedure ed_get_fimp_matsubara_site
     module procedure ed_get_fimp_matsubara_lattice
  end interface ed_get_fimp_matsubara

  interface ed_get_gimp_realaxis
     module procedure ed_get_gimp_realaxis_main
     module procedure ed_get_gimp_realaxis_site
     module procedure ed_get_gimp_realaxis_lattice
  end interface ed_get_gimp_realaxis

  interface ed_get_fimp_realaxis
     module procedure ed_get_fimp_realaxis_main
     module procedure ed_get_fimp_realaxis_site
     module procedure ed_get_fimp_realaxis_lattice
  end interface ed_get_fimp_realaxis


  !Retrieve imp GF_0 (G0_and) through routines.
  interface ed_get_g0imp_matsubara
     module procedure ed_get_g0imp_matsubara_main
  end interface ed_get_g0imp_matsubara

  interface ed_get_f0imp_matsubara
     module procedure ed_get_f0imp_matsubara_main
  end interface ed_get_f0imp_matsubara

  interface ed_get_g0imp_realaxis
     module procedure ed_get_g0imp_realaxis_main
  end interface ed_get_g0imp_realaxis

  interface ed_get_f0imp_realaxis
     module procedure ed_get_f0imp_realaxis_main
  end interface ed_get_f0imp_realaxis



  interface ed_get_gimp
     module procedure rebuild_gimp_single
     module procedure rebuild_gimp_ineq
  end interface ed_get_gimp

  interface ed_get_sigma
     module procedure rebuild_sigma_single
     module procedure rebuild_sigma_ineq
  end interface ed_get_sigma



  !Retrieve static common observables  
  interface ed_get_dens
     module procedure ed_get_dens_1
     module procedure ed_get_dens_2
     module procedure ed_get_dens_lattice_1
     module procedure ed_get_dens_lattice_2
  end interface ed_get_dens

  interface ed_get_mag
     module procedure ed_get_mag_1
     module procedure ed_get_mag_2
     module procedure ed_get_mag_lattice_1
     module procedure ed_get_mag_lattice_2
  end interface ed_get_mag

  interface ed_get_docc
     module procedure ed_get_docc_1
     module procedure ed_get_docc_2
     module procedure ed_get_docc_lattice_1
     module procedure ed_get_docc_lattice_2
  end interface ed_get_docc

  interface ed_get_phisc
     module procedure ed_get_phisc_1
     module procedure ed_get_phisc_2
     module procedure ed_get_phisc_lattice_1
     module procedure ed_get_phisc_lattice_2
  end interface ed_get_phisc

  interface ed_get_eimp
     module procedure :: ed_get_eimp_
     module procedure :: ed_get_eimp_lattice
  end interface ed_get_eimp

  interface ed_get_epot
     module procedure :: ed_get_epot_
     module procedure :: ed_get_epot_lattice
  end interface ed_get_epot

  interface ed_get_eint
     module procedure :: ed_get_eint_
     module procedure :: ed_get_eint_lattice
  end interface ed_get_eint

  interface ed_get_ehartree
     module procedure :: ed_get_ehartree_
     module procedure :: ed_get_ehartree_lattice
  end interface ed_get_ehartree

  interface ed_get_eknot
     module procedure :: ed_get_eknot_
     module procedure :: ed_get_eknot_lattice
  end interface ed_get_eknot

  interface ed_get_doubles
     module procedure :: ed_get_doubles_
     module procedure :: ed_get_doubles_lattice
  end interface ed_get_doubles

  interface ed_get_dust
     module procedure :: ed_get_dust_
     module procedure :: ed_get_dust_lattice
  end interface ed_get_dust

  interface ed_get_dund
     module procedure :: ed_get_dund_
     module procedure :: ed_get_dund_lattice
  end interface ed_get_dund

  interface ed_get_dse
     module procedure :: ed_get_dse_
     module procedure :: ed_get_dse_lattice
  end interface ed_get_dse

  interface ed_get_dph
     module procedure :: ed_get_dph_
     module procedure :: ed_get_dph_lattice
  end interface ed_get_dph

  interface ed_get_density_matrix
     module procedure :: ed_get_density_matrix_single
     module procedure :: ed_get_density_matrix_lattice
  end interface ed_get_density_matrix

  interface ed_read_impSigma
     module procedure :: ed_read_impSigma_single
     module procedure :: ed_read_impSigma_lattice
  end interface ed_read_impSigma


  public :: ed_get_sigma_matsubara
  public :: ed_get_self_matsubara
  public :: ed_get_sigma_realaxis
  public :: ed_get_self_realaxis

  public :: ed_get_gimp_matsubara
  public :: ed_get_fimp_matsubara
  public :: ed_get_gimp_realaxis
  public :: ed_get_fimp_realaxis

  public :: ed_get_g0imp_matsubara
  public :: ed_get_f0imp_matsubara
  public :: ed_get_g0imp_realaxis
  public :: ed_get_f0imp_realaxis

  public :: ed_get_g0and_function
  public :: ed_get_f0and_function
  public :: ed_get_delta_function
  public :: ed_get_fdelta_function

  public :: ed_get_gimp
  public :: ed_get_sigma

  public :: ed_get_dens
  public :: ed_get_mag
  public :: ed_get_docc
  public :: ed_get_phisc

  public :: ed_get_eimp
  public :: ed_get_epot
  public :: ed_get_eint 
  public :: ed_get_ehartree
  public :: ed_get_eknot

  public :: ed_get_doubles
  public :: ed_get_dust
  public :: ed_get_dund
  public :: ed_get_dse
  public :: ed_get_dph

  public :: ed_get_density_matrix
  public :: ed_get_quantum_SOC_operators_single
  public :: ed_get_quantum_SOC_operators_lattice

  public :: ed_get_neigen_total

  public :: ed_read_impSigma
  public :: ed_read_impGmatrix

  !****************************************************************************************!
  !****************************************************************************************!

  public :: ed_print_impSigma
  public :: ed_print_impG
  public :: ed_print_impG0
  public :: ed_print_impD
  public :: ed_print_impChi
  public :: ed_print_impGmatrix


  !****************************************************************************************!
  !****************************************************************************************!


  character(len=64)                :: suffix





contains


  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the impurity self-energy 
  !+-----------------------------------------------------------------------------+!
  include "get_sigma.f90"


  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the impurity green's functions 
  !+-----------------------------------------------------------------------------+!
  include "get_gimp.f90"

  !+--------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the impurity green's functions 
  !+--------------------------------------------------------------------------+!
  include "get_g0imp.f90"

  !+--------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve Anderson non-interacting green's functions 
  !+--------------------------------------------------------------------------+!
  include "get_gand_bath.f90"

  !+--------------------------------------------------------------------------+!
  ! PURPOSE: Re-build the Impurity green's functions and self-energy at
  !          arbitrary complex zeta
  !+--------------------------------------------------------------------------+!
  include "rebuild_gimp.f90"
  include "rebuild_sigma.f90"

  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the local observables
  !+-----------------------------------------------------------------------------+!
  include "get_dens.f90"
  include "get_mag.f90"
  include "get_docc.f90"
  include "get_phisc.f90"
  include "get_eimp.f90"
  include "get_doubles.f90"
  include "get_imp_dm.f90"
  include "get_imp_SOC_op.f90"
  include "get_lanc_info.f90"




  !+------------------------------------------------------------------+
  !PURPOSE  : Print impurity Functions case:
  ! - impSigma
  ! - impG
  ! - impG0
  ! NORMAL - SUPERConducting - NONSU2
  !+------------------------------------------------------------------+
  include "print_impSigma.f90"
  subroutine ed_print_impSigma
    call allocate_grids
    select case(ed_mode)
    case ("normal");call print_impSigma_normal
    case ("superc");call print_impSigma_superc
    case ("nonsu2");call print_impSigma_nonsu2
    case default;stop "ed_print_impSigma error: ed_mode not in the list"
    end select
    call deallocate_grids
  end subroutine ed_print_impSigma


  include "print_impG.f90"
  subroutine ed_print_impG
    call allocate_grids
    select case(ed_mode)
    case ("normal");call print_impG_normal
    case ("superc");call print_impG_superc
    case ("nonsu2");call print_impG_nonsu2
    case default;stop "ed_print_impG error: ed_mode not in the list"
    end select
    call deallocate_grids
  end subroutine ed_print_impG


  include "print_impG0.f90"
  subroutine ed_print_impG0
    call allocate_grids
    select case(ed_mode)
    case ("normal");call print_impG0_normal
    case ("superc");call print_impG0_superc
    case ("nonsu2");call print_impG0_nonsu2
    case default;stop "ed_print_impG0 error: ed_mode not in the list"
    end select
    call deallocate_grids
  end subroutine ed_print_impG0

  subroutine ed_print_impD
    call allocate_grids()
    !Print the impurity functions:
    call splot("impDph_iw.ed"   ,vm,impDmats_ph(:))
    call splot("impDph_realw.ed",vr,impDreal_ph(:))
    call deallocate_grids()
  end subroutine ed_print_impD


  include "print_impChi.f90"
  subroutine ed_print_impChi
    call allocate_grids
    if(chispin_flag)call print_chi_spin
    if(chidens_flag)call print_chi_dens
    if(chipair_flag)call print_chi_pair
    if(chiexct_flag)call print_chi_exct
    call deallocate_grids
  end subroutine ed_print_impChi


  subroutine ed_print_impGmatrix(file)
    character(len=*),optional :: file
    character(len=256)        :: file_
    if(.not.allocated(impGmatrix))stop "ED_PRINT_IMPGFMATRIX ERROR: impGmatrix not allocated!"
    file_="gfmatrix";if(present(file))file_=str(file)
    call write_GFmatrix(impGmatrix,str(file_)//str(ed_file_suffix)//".restart")
  end subroutine ed_print_impGmatrix





  ! PURPOSE: Read self-energy function(s) - also for inequivalent sites.
  !+-----------------------------------------------------------------------------+!

  include "read_impSigma.f90"
  subroutine ed_read_impSigma_single
    !
    if(allocated(impSmats))deallocate(impSmats)
    if(allocated(impSreal))deallocate(impSreal)
    if(allocated(impSAmats))deallocate(impSAmats)
    if(allocated(impSAreal))deallocate(impSAreal)
    allocate(impSmats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impSreal(Nspin,Nspin,Norb,Norb,Lreal))
    allocate(impSAmats(Nspin,Nspin,Norb,Norb,Lmats)) !THIS SHOULD NOT DEPEND ON SPIN: NSPIN=>1
    allocate(impSAreal(Nspin,Nspin,Norb,Norb,Lreal)) !THIS SHOULD NOT DEPEND ON SPIN: NSPIN=>1
    impSmats=zero
    impSreal=zero
    impSAmats=zero
    impSAreal=zero
    !
    select case(ed_mode)
    case ("normal");call read_impSigma_normal
    case ("superc");call read_impSigma_superc
    case ("nonsu2");call read_impSigma_nonsu2
    case default;stop "ed_read_impSigma error: ed_mode not in the list"
    end select
  end subroutine ed_read_impSigma_single

  subroutine ed_read_impSigma_lattice(Nineq)
    integer :: Nineq
    integer :: ilat
    !
    if(allocated(Smats_ineq))deallocate(Smats_ineq)
    if(allocated(Sreal_ineq))deallocate(Sreal_ineq)
    if(allocated(SAmats_ineq))deallocate(SAmats_ineq)
    if(allocated(SAreal_ineq))deallocate(SAreal_ineq)
    allocate(Smats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Sreal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(SAmats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(SAreal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
    Smats_ineq  = zero 
    Sreal_ineq  = zero 
    SAmats_ineq = zero 
    SAreal_ineq = zero
    !
    do ilat=1,Nineq
       ed_file_suffix=reg(ineq_site_suffix)//str(ilat,site_indx_padding)
       call ed_read_impSigma_single
       Smats_ineq(ilat,:,:,:,:,:)  = impSmats
       Sreal_ineq(ilat,:,:,:,:,:)  = impSreal
       SAmats_ineq(ilat,:,:,:,:,:) = impSAmats
       SAreal_ineq(ilat,:,:,:,:,:) = impSAreal
    enddo
    ed_file_suffix=""
  end subroutine ed_read_impSigma_lattice


  !+-------------------------------------------------------------------+
  !PURPOSE  : Read cluster GF from file
  !+-------------------------------------------------------------------+
  subroutine ed_read_impGmatrix(file)
    character(len=*),optional :: file
    character(len=256)        :: file_
    !
    if(allocated(impGmatrix))call deallocate_GFmatrix(impGmatrix)
    if(allocated(impGmatrix))deallocate(impGmatrix)
    allocate(impGmatrix(Nspin,Nspin,Norb,Norb))
    file_="gfmatrix";if(present(file))file_=str(file)
    call read_GFmatrix(impGmatrix,str(file_)//str(ed_file_suffix)//".restart")
  end subroutine ed_read_impGmatrix


END MODULE ED_IO







