MODULE ED_GREENS_FUNCTIONS
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: str,free_unit,reg,free_units,txtfy
  USE SF_LINALG,  only: inv,eigh,eye
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_IO                     !< this contains the routine to print GF,Sigma and G0
  USE ED_EIGENSPACE
  USE ED_BATH
  USE ED_BATH_FUNCTIONS
  USE ED_SETUP
  USE ED_HAMILTONIAN
  USE ED_AUX_FUNX
  !
  USE ED_GF_NORMAL
  USE ED_GF_SUPERC
  USE ED_GF_NONSU2
  !
  implicit none
  private 

  public :: buildGf_impurity
  public :: rebuildGf_impurity

contains


  !+------------------------------------------------------------------+
  ! GF CALCULATIONS
  !+------------------------------------------------------------------+
  subroutine buildGF_impurity()
    !
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG build_GF: build GFs"
#endif
    call allocate_grids
    !
    call deallocate_GFmatrix(impGmatrix)
    impGmats=zero
    impGreal=zero
    !
    impSmats = zero
    impSreal = zero
    !
    impG0mats=zero
    impG0real=zero
    !
    impDmats_ph=zero
    impDreal_ph=zero
    !
    write(LOGfile,"(A)")"Get impurity Greens functions:"
    select case(ed_mode)
    case default  ;call build_gf_normal()
    case("superc");call build_gf_superc()
    case("nonsu2");call build_gf_nonsu2()
    end select
    !
    select case(ed_mode)
    case default  ;call build_sigma_normal()
    case("superc");call build_sigma_superc()
    case("nonsu2");call build_sigma_nonsu2()
    end select
    !
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG build_GF: writing results"
    write(Logfile,"(A)")""
#endif
    if(MPIMASTER)then
       if(ed_print_Sigma)        call ed_print_impSigma()
       if(ed_print_G)            call ed_print_impG()
       if(ed_print_G.AND.DimPh>1)call ed_print_impD()
       if(ed_print_G0)           call ed_print_impG0()
    endif
    if(MPIMASTER)                call ed_print_impGmatrix()
    !
    call deallocate_GFmatrix(impGmatrix)
    call deallocate_grids
    !
  end subroutine buildGF_impurity








  subroutine rebuildGF_impurity()
    !
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG rebuild_GF: re-building GF"
#endif
    !
    call ed_read_impGmatrix()
    !
    if(allocated(impSmats))deallocate(impSmats)
    if(allocated(impSAmats))deallocate(impSAmats)
    allocate(impSmats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impSAmats(Nspin,Nspin,Norb,Norb,Lmats))
    impSmats=zero
    impSAmats=zero
    !
    if(allocated(impGmats))deallocate(impGmats)
    if(allocated(impFmats))deallocate(impFmats)
    allocate(impGmats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impFmats(Nspin,Nspin,Norb,Norb,Lmats))
    impGmats=zero
    impFmats=zero
    !
    if(allocated(impG0mats))deallocate(impG0mats)
    if(allocated(impF0mats))deallocate(impF0mats)
    allocate(impG0mats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impF0mats(Nspin,Nspin,Norb,Norb,Lmats))
    impG0mats=zero
    impF0mats=zero
    if(allocated(impSreal))deallocate(impSreal)
    if(allocated(impSAreal))deallocate(impSAreal)
    allocate(impSreal(Nspin,Nspin,Norb,Norb,Lreal))
    allocate(impSAreal(Nspin,Nspin,Norb,Norb,Lreal))
    impSreal=zero
    impSAreal=zero
    !
    if(allocated(impGreal))deallocate(impGreal)
    if(allocated(impFreal))deallocate(impFreal)
    allocate(impGreal(Nspin,Nspin,Norb,Norb,Lreal))
    allocate(impFreal(Nspin,Nspin,Norb,Norb,Lreal))
    impGreal=zero
    impFreal=zero
    !
    if(allocated(impG0real))deallocate(impG0real)
    if(allocated(impF0real))deallocate(impF0real)
    allocate(impG0real(Nspin,Nspin,Norb,Norb,Lreal))
    allocate(impF0real(Nspin,Nspin,Norb,Norb,Lreal))
    impG0real=zero
    impF0real=zero
    !
    call allocate_grids
    !
    write(LOGfile,"(A)")"Get impurity Greens functions:"
    select case(ed_mode)
    case default  ;call rebuild_gf_normal()
    case("superc");call rebuild_gf_superc()
    case("nonsu2");call rebuild_gf_nonsu2()
    end select
    !
    select case(ed_mode)
    case default  ;call build_sigma_normal()
    case("superc");call build_sigma_superc()
    case("nonsu2");call build_sigma_nonsu2()
    end select
    !
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG rebuild_GF: writing results"
    write(Logfile,"(A)")""
#endif
    if(MPIMASTER)then
       if(ed_print_Sigma)  call ed_print_impSigma()
       if(ed_print_G)      call ed_print_impG()
       if(ed_print_G0)     call ed_print_impG0()
    endif
    !
    call deallocate_GFmatrix(impGmatrix)
    call deallocate_grids
    !
  end subroutine rebuildGF_impurity









end MODULE ED_GREENS_FUNCTIONS
