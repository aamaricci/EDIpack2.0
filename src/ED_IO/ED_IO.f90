MODULE ED_IO
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_SETUP
  USE ED_BATH
  !
  USE SF_LINALG
  USE SF_SPIN
  USE SF_ARRAYS,  only: linspace,arange
  USE SF_IOTOOLS, only: str,reg,free_unit,splot,sread
  USE SF_MISC,    only: assert_shape
  implicit none
  private


  interface ed_get_sigma
!| This subrotine gets from the EDIpack2 library the value of the self-energy calculated 
! on the Matsubara or real-frequency axis, with number of frequencies :f:var:`lmats` or :f:var:`lreal` .
!| The self-energy is an array having the following possible dimensions:
!
!  * [:f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nspin`:math:`\cdot`:f:var:`norb`, :f:var:`lmats` / :f:var:`lreal`]
!  * [:f:var:`nlat` :math:`\cdot` :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nlat` :math:`\cdot` :f:var:`nspin` 
!    :math:`\cdot` :f:var:`norb`, :f:var:`lmats` / :f:var:`lreal`] 
!  * [:f:var:`nlat`, :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`lmats` / :f:var:`lreal`]
!  * [:f:var:`nspin`, :f:var:`nspin`, :f:var:`norb`, :f:var:`norb`, :f:var:`lmats` / :f:var:`lreal`]
!  * [:f:var:`nlat`, :f:var:`nspin`, :f:var:`nspin`, :f:var:`norb`, :f:var:`norb`, :f:var:`lmats` / :f:var:`lreal`]
!
     module procedure :: ed_get_sigma_site_n3
     module procedure :: ed_get_sigma_site_n5
     module procedure :: ed_get_sigma_lattice_n3
     module procedure :: ed_get_sigma_lattice_n4
     module procedure :: ed_get_sigma_lattice_n6
  end interface ed_get_sigma


  interface ed_get_gimp
!This subroutine gets from the EDIpack2 library the value of the impurity Green's function calculated 
!on the Matsubara or real-frequency axis, with number of frequencies :f:var:`lmats` or :f:var:`lreal` .
!
!The impurity Green's function is an array having the following possible dimensions:
!
!  * [:f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nspin`:math:`\cdot`:f:var:`norb`, :f:var:`lmats` / :f:var:`lreal`]
!  * [:f:var:`nlat` :math:`\cdot` :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nlat` :math:`\cdot` :f:var:`nspin` 
!    :math:`\cdot` :f:var:`norb`, :f:var:`lmats` / :f:var:`lreal`] 
!  * [:f:var:`nlat`, :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`lmats` / :f:var:`lreal`]
!  * [:f:var:`nspin`, :f:var:`nspin`, :f:var:`norb`, :f:var:`norb`, :f:var:`lmats` / :f:var:`lreal`]
!  * [:f:var:`nlat`, :f:var:`nspin`, :f:var:`nspin`, :f:var:`norb`, :f:var:`norb`, :f:var:`lmats` / :f:var:`lreal`]
!
     module procedure :: ed_get_gimp_site_n3
     module procedure :: ed_get_gimp_site_n5
     module procedure :: ed_get_gimp_lattice_n3
     module procedure :: ed_get_gimp_lattice_n4
     module procedure :: ed_get_gimp_lattice_n6
  end interface ed_get_gimp



  interface ed_get_g0imp
!| This subroutine gets from the EDIpack2 library the value of the impurity non-interacting Green's function calculated 
! on the Matsubara or real-frequency axis, with number of frequencies :f:var:`lmats` or :f:var:`lreal` .
!| It autonomously decides whether the system is single-impurity or real-space DMFT based on the :f:var:`bath` shape
!
!The impurity non-interacting Green's function is an array having the following possible dimensions:
! 
!  * [:f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nspin`:math:`\cdot`:f:var:`norb`, :f:var:`lmats` / :f:var:`lreal`]  
!  * [:f:var:`nlat` :math:`\cdot` :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nlat` :math:`\cdot` :f:var:`nspin` 
!    :math:`\cdot` :f:var:`norb`, :f:var:`lmats` / :f:var:`lreal`] 
!  * [:f:var:`nlat`, :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`lmats` / :f:var:`lreal`]
!  * [:f:var:`nspin`, :f:var:`nspin`, :f:var:`norb`, :f:var:`norb`, :f:var:`lmats` / :f:var:`lreal`]
!  * [:f:var:`nlat`, :f:var:`nspin`, :f:var:`nspin`, :f:var:`norb`, :f:var:`norb`, :f:var:`lmats` / :f:var:`lreal`]
!
!The bath is an array having the following possible dimensions:
!
!  * [:f:var:`nb`] for single-impurity DMFT
!  * [:f:var:`nlat`, :f:var:`nb`] for real-space DMFT, with :f:var:`nlat` the number of inequivalent impurity sites
!
!Where :f:var:`nb` is the length of the :f:var:`bath` array.
 !
     module procedure :: ed_get_g0imp_site_n3
     module procedure :: ed_get_g0imp_site_n5
     module procedure :: ed_get_g0imp_lattice_n3
     module procedure :: ed_get_g0imp_lattice_n4
     module procedure :: ed_get_g0imp_lattice_n6
  end interface ed_get_g0imp


  interface ed_build_gimp
!| This subroutine returns to the user the impurity Green's function matrix calculated at any provided frequency
! in the complex plane, by obtaining it from the stored poles and weights.
!
!The impurity Green's function is an array having the following possible dimensions:
! 
!  * [:f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nspin`:math:`\cdot`:f:var:`norb`, :code:`size(zeta)`]  
!  * [:f:var:`nlat` :math:`\cdot` :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nlat` :math:`\cdot` :f:var:`nspin` 
!    :math:`\cdot` :f:var:`norb`, :code:`size(zeta)`]  
!  * [:f:var:`nlat`, :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :code:`size(zeta)`] 
!  * [:f:var:`nspin`, :f:var:`nspin`, :f:var:`norb`, :f:var:`norb`, :code:`size(zeta)`]
!  * [:f:var:`nlat`, :f:var:`nspin`, :f:var:`nspin`, :f:var:`norb`, :f:var:`norb`, :code:`size(zeta)`]
!
     module procedure :: rebuild_gimp_single_n3
     module procedure :: rebuild_gimp_single_n5
     module procedure :: rebuild_gimp_ineq_n3
     module procedure :: rebuild_gimp_ineq_n4
     module procedure :: rebuild_gimp_ineq_n6
  end interface ed_build_gimp

  interface ed_build_sigma
!| This subroutine returns to the user the self-energy matrix calculated at any provided frequency
! in the complex plane, by obtaining it from the stored poles and weights
!
!The self-energy is an array having the following possible dimensions:
! 
!  * [:f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nspin`:math:`\cdot`:f:var:`norb`, :code:`size(zeta)`]  
!  * [:f:var:`nlat` :math:`\cdot` :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nlat` :math:`\cdot` :f:var:`nspin` 
!    :math:`\cdot` :f:var:`norb`, :code:`size(zeta)`]  
!  * [:f:var:`nlat`, :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :code:`size(zeta)`] 
!  * [:f:var:`nspin`, :f:var:`nspin`, :f:var:`norb`, :f:var:`norb`, :code:`size(zeta)`]
!  * [:f:var:`nlat`, :f:var:`nspin`, :f:var:`nspin`, :f:var:`norb`, :f:var:`norb`, :code:`size(zeta)`]
!
     module procedure :: rebuild_sigma_single_n3
     module procedure :: rebuild_sigma_single_n5
     module procedure :: rebuild_sigma_ineq_n3
     module procedure :: rebuild_sigma_ineq_n4
     module procedure :: rebuild_sigma_ineq_n6
  end interface ed_build_sigma


  !Build Gand/Delta from a user bath
  interface ed_get_g0and
     module procedure :: ed_get_g0and_n3
     module procedure :: ed_get_g0and_n5
  end interface ed_get_g0and

  interface ed_get_delta
     module procedure :: ed_get_delta_n3
     module procedure :: ed_get_delta_n5
  end interface ed_get_delta


  !Observables
  interface ed_get_dens
 !This subroutine gets from the EDIpack2 library the value of the charge density and passes it to the user.
 !
 !The :f:var:`self` variable can have the following dimensions:
 ! 
 !  * scalar: if :f:var:`iorb` is provided for single-impurity DMFT, density for that orbital
 !  * [:f:var:`norb`]: if no optional variable is provided for single-impurity DMFT, density for all orbitals
 !  * [:f:var:`nlat`]: if :f:var:`iorb` (default = 1) is provided for real-space DMFT with :f:var:`nlat` impurities, density for that orbital for all impurity sites
 !  * [:f:var:`nlat`, :f:var:`norb`]: if :f:var:`nlat` is provided for real-space DMFT, density for all impurity sites and orbitals
 !
     module procedure :: ed_get_dens_n0
     module procedure :: ed_get_dens_n1
     module procedure :: ed_get_dens_n2
  end interface ed_get_dens

  interface ed_get_mag
 !This subroutine gets from the EDIpack2 library the value of the magnetization and passes it to the user.
 !
 !The :f:var:`self` variable can have the following dimensions:
 ! 
 !  * scalar: if :f:var:`component` and :f:var:`iorb` are provided for single-impurity DMFT, given magnetization component for that orbital
 !  * [:f:var:`norb`]: for single-impurity DMFT, one magnetization component for all orbitals
 !  * [:f:var:`nlat`]: for real-space DMFT with :f:var:`nlat` impurities, magnetization for that orbital for all impurity sites
 !  * [:f:var:`nlat`, :f:var:`norb`]: if :f:var:`nlat` is provided for real-space DMFT, one magnetization component for all orbitals and impurity sites
 !  * [:f:var:`nlat`, :code:`3`, :f:var:`norb`]: if :f:var:`nlat` is provided for real-space DMFT, all magnetization components for all orbitals and sites
 !
     module procedure :: ed_get_mag_n0
     module procedure :: ed_get_mag_n1
     module procedure :: ed_get_mag_n2
     module procedure :: ed_get_mag_n3
  end interface ed_get_mag

  interface ed_get_docc
 !This subroutine gets from the EDIpack2 library the value of the double occupation and passes it to the user.
 !
 !The :f:var:`self` variable can have the following dimensions:
 ! 
 !  * scalar: if :f:var:`iorb` is provided for single-impurity DMFT, dobule-occupation for that orbital
 !  * [:f:var:`norb`]: if no optional variable is provided for single-impurity DMFT, double-occupation for all orbitals
 !  * [:f:var:`nlat`]: if :f:var:`iorb` (default = 1) is provided for real-space DMFT with :f:var:`nlat` impurities, double-occupation for that orbital for all impurity sites
 !  * [:f:var:`nlat`, :f:var:`norb`]: if :f:var:`nlat` is provided for real-space DMFT, double-occupation for all impurity sites and orbitals
 !
     module procedure :: ed_get_docc_n0
     module procedure :: ed_get_docc_n1
     module procedure :: ed_get_docc_n2
  end interface ed_get_docc

  interface ed_get_phi
     module procedure :: ed_get_phisc_n0
     module procedure :: ed_get_phisc_n1
     module procedure :: ed_get_phisc_n2
  end interface ed_get_phi


  !Get Energies
  interface ed_get_eimp
!This subroutine gets from the EDIpack2 library and passes to the user the array [ :f:var:`ed_epot` , :f:var:`ed_eint` , :f:var:`ed_ehartree` , :f:var:`ed_eknot` ].
!These are the expectation values various contribution to the internal energy
!
!  * :f:var:`ed_epot` = energy contribution from the interaction terms, **including** the Hartree term
!  * :f:var:`ed_eint` = energy contribution from the interaction terms, **excluding** the Hartree term
!  * :f:var:`ed_ehartree` = :math:`-\frac{U}{2} \sum_{i} \langle n_{i\uparrow} + n_{i\downarrow} \rangle 
!    -\frac{2U^{'}-J_{H}}{2} \sum_{i < j} \langle n_{i\uparrow}+n_{i\downarrow} + n_{i\downarrow}+n_{j\downarrow} \rangle
!    +\frac{U}{4} + \frac{2U^{'}-J_{H}}{2}` for :math:`i,j` orbitals
!  * :f:var:`ed_eknot` = kinetic term from the **local** 1-body Hamiltonian
!
!The returned array can have the following dimensions:
!
!  * [:code:`4`]: for single-site DMFT
!  * [:f:var:`nlat`, :code:`4`]: for real-space DMFT with :f:var:`nlat` impurities
!
     module procedure :: ed_get_eimp_n1
     module procedure :: ed_get_eimp_n2
  end interface ed_get_eimp

  interface ed_get_epot
!This subroutine gets from the EDIpack2 library and passes to the user the value of 
!:f:var:`ed_epot`, the energy contribution from the interaction terms, **including** the Hartree term.
!The returned array can have the following dimensions:
!
!  * scalar: for single-site DMFT
!  * [:f:var:`nlat`]: for real-space DMFT with :f:var:`nlat` impurities
!
     module procedure :: ed_get_epot_n0
     module procedure :: ed_get_epot_n1
  end interface ed_get_epot

  interface ed_get_eint
!This subroutine gets from the EDIpack2 library and passes to the user the value of 
!:f:var:`ed_int`, the energy contribution from the interaction terms, **excluding** the Hartree term.
!The returned array can have the following dimensions:
!
!  * scalar: for single-site DMFT
!  * [:f:var:`nlat`]: for real-space DMFT with :f:var:`nlat` impurities
!
     module procedure :: ed_get_eint_n0
     module procedure :: ed_get_eint_n1
  end interface ed_get_eint

  interface ed_get_ehartree
!This subroutine gets from the EDIpack2 library and passes to the user the value of the Hartree potential 
!:f:var:`ed_ehartree`. The returned array can have the following dimensions:
!
!  * scalar: for single-site DMFT
!  * [:f:var:`nlat`]: for real-space DMFT with :f:var:`nlat` impurities
!
     module procedure :: ed_get_ehartree_n0
     module procedure :: ed_get_ehartree_n1
  end interface ed_get_ehartree

  interface ed_get_eknot
!This subroutine gets from the EDIpack2 library and passes to the user the value
!:f:var:`ed_eknot`, the kinetic term from the **local** 1-body Hamiltonian
!The returned array can have the following dimensions:
!
!  * scalar: for single-site DMFT
!  * [:f:var:`nlat`]: for real-space DMFT with :f:var:`nlat` impurities
!
     module procedure :: ed_get_eknot_n0
     module procedure :: ed_get_eknot_n1
  end interface ed_get_eknot


  interface ed_get_doubles
!This subroutine gets from the EDIpack2 library and passes to the user the array [ :f:var:`ed_dust` , :f:var:`ed_dund` , :f:var:`ed_dse` , :f:var:`ed_dph` ].
!These are the expectation values of the two-body operators associated with the density-density inter-orbital interaction (with opposite and parallel spins), 
!spin-exchange and pair-hopping.
!
!  * :f:var:`ed_dust` = :math:`\sum_{i < j} n_{i\uparrow}n_{j\downarrow} + n_{i\downarrow}n_{j\uparrow}` for :math:`i,j` orbitals
!  * :f:var:`ed_dund` = :math:`\sum_{i < j} n_{i\uparrow}n_{j\uparrow}  + n_{i\downarrow}n_{j\downarrow}` for :math:`i,j` orbitals
!  * :f:var:`ed_dse` = :math:`\sum_{i < j} c^{\dagger}_{i\uparrow}c^{\dagger}_{j\uparrow}c_{i\downarrow}c_{j\uparrow}` for :math:`i,j` orbitals
!  * :f:var:`ed_dph` = :math:`\sum_{i < j} c^{\dagger}_{i\uparrow}c^{\dagger}_{i\downarrow}c_{j\downarrow}c_{j\uparrow}` for :math:`i,j` orbitals
!
!The returned array can have the following dimensions:
!
!  * [:code:`4`]: for single-site DMFT
!  * [:f:var:`nlat`, :code:`4`]: for real-space DMFT with :f:var:`nlat` impurities
!
     module procedure :: ed_get_doubles_n1
     module procedure :: ed_get_doubles_n2
  end interface ed_get_doubles

  interface ed_get_dust
!This subroutine gets from the EDIpack2 library and passes to the user the value of 
!:f:var:`ed_dust` = :math:`\sum_{i < j} n_{i\uparrow}n_{j\downarrow} + n_{i\downarrow}n_{j\uparrow}` for :math:`i,j` orbitals
!The returned array can have the following dimensions:
!
!  * scalar: for single-site DMFT
!  * [:f:var:`nlat`]: for real-space DMFT with :f:var:`nlat` impurities
!
     module procedure :: ed_get_dust_n0
     module procedure :: ed_get_dust_n1
  end interface ed_get_dust

  interface ed_get_dund
!This subroutine gets from the EDIpack2 library and passes to the user the value of 
!:f:var:`ed_dund` = :math:`\sum_{i < j} n_{i\uparrow}n_{j\uparrow}  + n_{i\downarrow}n_{j\downarrow}` for :math:`i,j` orbitals
!The returned array can have the following dimensions:
!
!  * scalar: for single-site DMFT
!  * [:f:var:`nlat`]: for real-space DMFT with :f:var:`nlat` impurities
!
     module procedure :: ed_get_dund_n0
     module procedure :: ed_get_dund_n1
  end interface ed_get_dund

  interface ed_get_dse
!This subroutine gets from the EDIpack2 library and passes to the user the value of 
!:f:var:`ed_dse` = :math:`\sum_{i < j} c^{\dagger}_{i\uparrow}c^{\dagger}_{j\uparrow}c_{i\downarrow}c_{j\uparrow}` for :math:`i,j` orbitals
!The returned array can have the following dimensions:
!
!  * scalar: for single-site DMFT
!  * [:f:var:`nlat`]: for real-space DMFT with :f:var:`nlat` impurities
!
     module procedure :: ed_get_dse_n0
     module procedure :: ed_get_dse_n1
  end interface ed_get_dse

  interface ed_get_dph
!This subroutine gets from the EDIpack2 library and passes to the user the value of 
!:f:var:`ed_dph` = :math:`\sum_{i < j} c^{\dagger}_{i\uparrow}c^{\dagger}_{i\downarrow}c_{j\downarrow}c_{j\uparrow}` for :math:`i,j` orbitals
!The returned array can have the following dimensions:
!
!  * scalar: for single-site DMFT
!  * [:f:var:`nlat`]: for real-space DMFT with :f:var:`nlat` impurities
!
     module procedure :: ed_get_dph_n0
     module procedure :: ed_get_dph_n1
  end interface ed_get_dph



  

  interface ed_get_density_matrix
     module procedure :: ed_get_density_matrix_single
     module procedure :: ed_get_density_matrix_lattice
  end interface ed_get_density_matrix

  interface ed_read_impSigma
     module procedure :: ed_read_impSigma_single
     module procedure :: ed_read_impSigma_lattice
  end interface ed_read_impSigma


  public :: ed_get_sigma
  public :: ed_get_gimp
  public :: ed_get_g0imp
  public :: ed_get_g0and
  public :: ed_get_delta

  public :: ed_build_gimp
  public :: ed_build_sigma

  public :: ed_get_dens
  public :: ed_get_mag
  public :: ed_get_docc
  public :: ed_get_phi
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
  public :: ed_get_neigen_total
  !  
  public :: ed_get_density_matrix
  public :: ed_get_quantum_SOC_operators_single
  public :: ed_get_quantum_SOC_operators_lattice


  !****************************************************************************************!
  !****************************************************************************************!

  public :: ed_print_impSigma
  public :: ed_print_impG
  public :: ed_print_impG0
  public :: ed_print_impD
  public :: ed_print_impChi
  public :: ed_print_impGmatrix


  public :: ed_read_impSigma
  public :: ed_read_impGmatrix


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
  include "rebuild_impG.f90"
  include "rebuild_impSigma.f90"


  !+--------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve impurity density matrices and SOC operators
  !+--------------------------------------------------------------------------+!
  include "get_imp_dm.f90"
  include "get_imp_SOC_op.f90"

  !+--------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the local observables
  !+--------------------------------------------------------------------------+!
  include "get_dens.f90"
  include "get_mag.f90"
  include "get_docc.f90"
  include "get_phi.f90"
  include "get_energy.f90"
  include "get_doubles.f90"
  include "get_neigen.f90"






  !+------------------------------------------------------------------+
  !PURPOSE  : Print impurity Functions case INTERNAL USE
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




  !+--------------------------------------------------------------------------+!
  ! PURPOSE: Read self-energy function(s) - also for inequivalent sites.
  !+--------------------------------------------------------------------------+!
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







