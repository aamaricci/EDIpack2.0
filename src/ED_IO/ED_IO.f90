MODULE ED_IO
  !Contains a set of routines that retrieve quantities such as Green's functions, self-energies (see :f:mod:`ed_greens_functions` ) and observables (from :f:mod:`ed_observables` ) and pass them to the user, as well ass routines to read and store Green's function and self-energies.
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
     !| This subroutine returns to the user the normal non-interacting Green's function :math:`G_0(x)` and
     ! the anomalous non-interacting Green's function :math:`F_0(x)` on a given set of frequencies. It does so
     ! by calling :f:func:`g0and_bath_function` and :f:func:`g0and_bath_function`.
     !
     !The non-interacting Green's function is an array having the following possible dimensions:
     ! 
     !  * [:f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nspin`:math:`\cdot`:f:var:`norb`, :code:`size(x)`]  
     !  * [:f:var:`nspin`, :f:var:`nspin`, :f:var:`norb`, :f:var:`norb`, :code:`size(x)`]
     !
     module procedure :: ed_get_g0and_n3
     module procedure :: ed_get_g0and_n5
  end interface ed_get_g0and

  interface ed_get_delta
     !| This subroutine returns to the user the normal hybridization function :math:`\Delta(x)` and
     ! the anomalous hybridization function :math:`\Theta(x)` on a given set of frequencies. It does so
     ! by calling :f:func:`delta_bath_function` and :f:func:`fdelta_bath_function`.
     !
     !The hybridization function is an array having the following possible dimensions:
     ! 
     !  * [:f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nspin`:math:`\cdot`:f:var:`norb`, :code:`size(x)`]  
     !  * [:f:var:`nspin`, :f:var:`nspin`, :f:var:`norb`, :f:var:`norb`, :code:`size(x)`]
     !
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
     !This subroutine gets from the EDIpack2 library the value of the superconducting order parameter :math:`\phi` ( :f:var:`ed_mode` = :code:`superc` ) and passes it to the user.
     !
     !The :f:var:`self` variable can have the following dimensions:
     ! 
     !  * scalar: if :f:var:`iorb` is provided for single-impurity DMFT, :math:`\phi` for that orbital
     !  * [:f:var:`norb`]: if no optional variable is provided for single-impurity DMFT, :math:`\phi` for all orbitals
     !  * [:f:var:`nlat`]: if :f:var:`iorb` (default = 1) is provided for real-space DMFT with :f:var:`nlat` impurities, :math:`\phi` for that orbital for all impurity sites
     !  * [:f:var:`nlat`, :f:var:`norb`]: if :f:var:`nlat` is provided for real-space DMFT, :math:`\phi` for all impurity sites and orbitals
     !
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
     !This subroutine returns to the user the impurity density matrix.
     !The density matrix is an array having the following possible dimensions:
     ! 
     !  * [:f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nspin`:math:`\cdot`:f:var:`norb`]  for single-impurity DMFT
     !  * [:code:`nlat`, :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nspin`:math:`\cdot`:f:var:`norb`] for real-space DMFT
     !
     module procedure :: ed_get_density_matrix_single
     module procedure :: ed_get_density_matrix_lattice
  end interface ed_get_density_matrix

  interface ed_read_impSigma
     !This subroutine reads the impurity Sigmas from files in the execution folder and stores them in the global variables 
     ! 
     !  * :f:var:`impsmats` normal self-energy, Matsubara axis
     !  * :f:var:`impsreal` normal self-energy, real frequency axis
     !  * :f:var:`impsamats` anomalous self-energy, Matsubara axis
     !  * :f:var:`impsareal` anomalous self-energy, real frequency axis
     !  * :f:var:`smats_ineq` normal self-energy, Matsubara axis, real-space DMFT
     !  * :f:var:`sreal_ineq` normal self-energy, real frequency axis, real-space DMFT
     !  * :f:var:`samats_ineq` anomalous self-energy, Matsubara axis, real-space DMFT
     !  * :f:var:`sareal_ineq` anomalous self-energy, real frequency axis, real-space DMFT
     !
     !The files have to be formatted to be compatible with the EDIpack2 library, that is :math:`[\omega,\mathrm{Im}\Sigma,\mathrm{Re}\Sigma]` .
     !One file per self-energy component, with the name
     !
     !  * :code:`"impSigma_l"//str(iorb)[str(jorb)]//_s"//str(ispin)[str(jspin)]"_iw"//reg(ed_file_suffix)//".ed"` normal self-energy, Matsubara axis
     !  * :code:`"impSigma_l"//str(iorb)[str(jorb)]//_s"//str(ispin)[str(jspin)]"_realw"//reg(ed_file_suffix)//".ed"` normal self-energy, real frequency axis
     !  * :code:`"impSelf_l"//str(iorb)[str(jorb)]//_s"//str(ispin)"_iw"//reg(ed_file_suffix)//".ed"` anomalous self-energy, Matsubara axis
     !  * :code:`"impSelf_l"//str(iorb)[str(jorb)]//_s"//str(ispin)"_realw"//reg(ed_file_suffix)//".ed"` anomalous self-energy, real frequency axis
     !
     !The variable :f:var:`ed_file_suffix` is :code:`"_ineq_Nineq"` padded with 4 zeros in the case of inequivalent sites, as per documentation.
     !
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
  public :: ed_get_quantum_SOC_operators

  !****************************************************************************************!
  !****************************************************************************************!

  public :: ed_print_impSigma
  public :: ed_print_impG
  public :: ed_print_impG0
  public :: ed_print_impDelta
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






  !+------------------------------------------------------------------+
  !PURPOSE  : Print impurity Functions case INTERNAL USE
  ! - impSigma
  ! - impG
  ! - impG0
  ! NORMAL - SUPERConducting - NONSU2
  !+------------------------------------------------------------------+
  include "print_impSigma.f90"
  subroutine ed_print_impSigma
    !This subroutine print the impurity self-energy on plain text files in the execution folder.
    !The files are formatted like :math:`[\omega,\mathrm{Im}\Sigma,\mathrm{Re}\Sigma]` .
    !One file per self-energy component, with the name
    !
    !  * :code:`"impSigma_l"//str(iorb)[str(jorb)]//_s"//str(ispin)[str(jspin)]"_iw"//reg(ed_file_suffix)//".ed"` normal self-energy, Matsubara axis
    !  * :code:`"impSigma_l"//str(iorb)[str(jorb)]//_s"//str(ispin)[str(jspin)]"_realw"//reg(ed_file_suffix)//".ed"` normal self-energy, real frequency axis
    !  * :code:`"impSelf_l"//str(iorb)[str(jorb)]//_s"//str(ispin)"_iw"//reg(ed_file_suffix)//".ed"` anomalous self-energy, Matsubara axis
    !  * :code:`"impSelf_l"//str(iorb)[str(jorb)]//_s"//str(ispin)"_realw"//reg(ed_file_suffix)//".ed"` anomalous self-energy, real frequency axis
    !
    !The variable :f:var:`ed_file_suffix` is :code:`"_ineq_Nineq"` padded with 4 zeros in the case of inequivalent sites, as per documentation.
    !
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
    !This subroutine print the impurity Green's function on plain text files in the execution folder.
    !The files are formatted like :math:`[\omega,\mathrm{Im}G,\mathrm{Re}G]` .
    !One file per Green'sfunction component, with the name
    !
    !  * :code:`"impG_l"//str(iorb)[str(jorb)]//_s"//str(ispin)[str(jspin)]"_iw"//reg(ed_file_suffix)//".ed"` normal G, Matsubara axis
    !  * :code:`"impG_l"//str(iorb)[str(jorb)]//_s"//str(ispin)[str(jspin)]"_realw"//reg(ed_file_suffix)//".ed"` normal G, real frequency axis
    !  * :code:`"impF_l"//str(iorb)[str(jorb)]//_s"//str(ispin)"_iw"//reg(ed_file_suffix)//".ed"` anomalous G, Matsubara axis
    !  * :code:`"impF_l"//str(iorb)[str(jorb)]//_s"//str(ispin)"_realw"//reg(ed_file_suffix)//".ed"` anomalous G, real frequency axis
    !
    !The variable :f:var:`ed_file_suffix` is :code:`"_ineq_Nineq"` padded with 4 zeros in the case of inequivalent sites, as per documentation.
    !
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
    !This subroutine print the non-interacting impurity Green's function on plain text files in the execution folder.
    !The files are formatted like :math:`[\omega,\mathrm{Im}G_{0},\mathrm{Re}G_{0}]` .
    !One file per Green's function component, with the name
    !
    !  * :code:`"impG0_l"//str(iorb)[str(jorb)]//_s"//str(ispin)[str(jspin)]"_iw"//reg(ed_file_suffix)//".ed"` normal G, Matsubara axis
    !  * :code:`"impG0_l"//str(iorb)[str(jorb)]//_s"//str(ispin)[str(jspin)]"_realw"//reg(ed_file_suffix)//".ed"` normal G, real frequency axis
    !  * :code:`"impF0_l"//str(iorb)[str(jorb)]//_s"//str(ispin)"_iw"//reg(ed_file_suffix)//".ed"` anomalous G, Matsubara axis
    !  * :code:`"impF0_l"//str(iorb)[str(jorb)]//_s"//str(ispin)"_realw"//reg(ed_file_suffix)//".ed"` anomalous G, real frequency axis
    !
    !The variable :f:var:`ed_file_suffix` is :code:`"_ineq_Nineq"` padded with 4 zeros in the case of inequivalent sites, as per documentation.
    !
    call allocate_grids
    select case(ed_mode)
    case ("normal");call print_impG0_normal
    case ("superc");call print_impG0_superc
    case ("nonsu2");call print_impG0_nonsu2
    case default;stop "ed_print_impG0 error: ed_mode not in the list"
    end select
    call deallocate_grids
  end subroutine ed_print_impG0



  include "print_impDelta.f90"
  subroutine ed_print_impDelta
    !This subroutine print the Hybridzation function on plain text files in the execution folder.
    !The files are formatted like :math:`[\omega,\mathrm{Im}\Delta_{0},\mathrm{Re}\Delta_{0}]` .
    !One file per Green's function component, with the name
    !
    !  * :code:`"impDelta_l"//str(iorb)[str(jorb)]//_s"//str(ispin)[str(jspin)]"_iw"//reg(ed_file_suffix)//".ed"` normal, Matsubara axis
    !  * :code:`"impDelta_l"//str(iorb)[str(jorb)]//_s"//str(ispin)[str(jspin)]"_realw"//reg(ed_file_suffix)//".ed"` normal, real frequency axis
    !  * :code:`"impFDelta_l"//str(iorb)[str(jorb)]//_s"//str(ispin)"_iw"//reg(ed_file_suffix)//".ed"` anomalous, Matsubara axis
    !  * :code:`"impFDelta_l"//str(iorb)[str(jorb)]//_s"//str(ispin)"_realw"//reg(ed_file_suffix)//".ed"` anomalous, real frequency axis
    !
    !The variable :f:var:`ed_file_suffix` is :code:`"_ineq_Nineq"` padded with 4 zeros in the case of inequivalent sites, as per documentation.
    !
    call allocate_grids
    select case(ed_mode)
    case ("normal");call print_impDelta_normal
    case ("superc");call print_impDelta_superc
    case ("nonsu2");call print_impDelta_nonsu2
    case default;stop "ed_print_impG0 error: ed_mode not in the list"
    end select
    call deallocate_grids
  end subroutine ed_print_impDelta



  subroutine ed_print_impD
    !This subroutine print the impurity phonon self-energy on the files
    !  * :code:`"impDph_iw.ed"`  matsubara axis
    !  * :code:`impDph_realw.ed"` real frequency axis
    !
    call allocate_grids()
    !Print the impurity functions:
    call splot("impDph_iw.ed"   ,vm,impDmats_ph(:))
    call splot("impDph_realw.ed",vr,impDreal_ph(:))
    call deallocate_grids()
  end subroutine ed_print_impD


  include "print_impChi.f90"
  subroutine ed_print_impChi
    !This subroutine prints the susceptibilities.
    !The files are formatted like :math:`[\omega,\mathrm{Im}\\chi,\mathrm{Re}\\chi]` .
    !Which susceptibilities are printed depends on the values of :f:var:`chispin_flag` (spin), :f:var:`chidens_flag` (charge), :f:var:`chipair_flag` (pair), :f:var:`chiexct_flag` (exciton).
    !One file per component. The name of the files are
    !
    !  * :code:`"[spin/dens/pair/exct]Chi_[singlet/tripletXY,tripletZ]_l"//str(iorb)[str(jorb)]//_s"//str(ispin)[str(jspin)]"_tau"//reg(ed_file_suffix)//".ed"` imaginary time
    !  * :code:`"[spin/dens/pair/exct]Chi_[singlet/tripletXY,tripletZ]_l"//str(iorb)[str(jorb)]//_s"//str(ispin)[str(jspin)]"_iw"//reg(ed_file_suffix)//".ed"` Matsubara axis
    !  * :code:`"[spin/dens/pair/exct]Chi_[singlet/tripletXY,tripletZ]_l"//str(iorb)[str(jorb)]//_s"//str(ispin)[str(jspin)]"_realw"//reg(ed_file_suffix)//".ed"` real frequency axis axis
    !
    !The variable :f:var:`ed_file_suffix` is :code:`"_ineq_Nineq"` padded with 4 zeros in the case of inequivalent sites, as per documentation.
    !
    call allocate_grids
    if(chispin_flag)call print_chi_spin
    if(chidens_flag)call print_chi_dens
    if(chipair_flag)call print_chi_pair
    if(chiexct_flag)call print_chi_exct
    call deallocate_grids
  end subroutine ed_print_impChi


  subroutine ed_print_impGmatrix(file)
    !This subroutine prints weights and poles of the impurity Green's function by calling :f:func:`write_GFmatrix`. These are stored
    !one a file named :code:`"file"//str(ed_file_suffix)//.restart"` taking into account the value of the global variable :f:var:`ed_file_suffix` ,
    !which is :code:`"_ineq_Nineq"` padded with 4 zeros in the case of inequivalent sites, as per documentation
    !
    character(len=*),optional :: file !filename prefix (default :code:`gfmatrix`)
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
    integer :: Nineq !number of inequivalent impurity sites for real-space DMFT
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


  !+-------------------------------------------------------------------+
  !PURPOSE  : get SOC operators
  !+-------------------------------------------------------------------+

  subroutine ed_get_quantum_soc_operators()
    !This subroutine gets and prints the values of the components :math:`\overrightarrow{L}`, :math:`\overrightarrow{S}`, :math:`\overrightarrow{J}`
    !in the chosen basis depending on :f:var:`jz_basis`, and prints them on the files :code:`"L_imp_"//reg(str(ndx))//".dat"` , 
    !:code:`"S_imp_"//reg(str(ndx))//".dat"` and :code:`"J_imp_"//reg(str(ndx))//".dat"` , where :code:`ndx` is the inequivalent
    !impurity site for real-space DMFT (if that is the case). The ordering of the results in the output files is described by comments
    !in the files themselves
    !
    if (allocated(imp_density_matrix_ineq)) then
       call ed_get_quantum_soc_operators_lattice()
    else
       call ed_get_quantum_SOC_operators_single()
    endif
  end subroutine ed_get_quantum_soc_operators

  !+----------------------------------------------------------------------------+
  !PURPOSE  : get number of spectrum eigenstates for inequivalent impurity sites
  !+----------------------------------------------------------------------------+ 


  subroutine ed_get_neigen_total(nlii,Nlat) 
    !In the case of inequivalent impurity sites, this function returns the number of eigenstates per impurity
    !site in the ED spectrum.
    integer                      :: Nlat !number of inequivalent impurity sites for real-space DMFT
    integer,dimension(Nlat)      :: nlii !array containing the number of eigenstates per inequivalent impurity site
    nlii=0d0
    if(allocated(neigen_total_ineq))then
       if(Nlat>size(neigen_total_ineq)) stop "ed_get_neigen_total error: required N_sites > evaluated N_sites"
       nlii=neigen_total_ineq
    endif
  end subroutine ed_get_neigen_total


END MODULE ED_IO







