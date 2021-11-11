module ED_MAIN
  USE SF_IOTOOLS, only: str,reg
  USE SF_TIMER,only: start_timer,stop_timer
  USE SF_MISC,only: assert_shape
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE, only: state_list,es_delete_espace
  USE ED_AUX_FUNX
  USE ED_SETUP
  USE ED_BATH
  USE ED_HAMILTONIAN
  USE ED_GREENS_FUNCTIONS
  USE ED_CHI_FUNCTIONS
  USE ED_OBSERVABLES
  USE ED_DIAG

  implicit none
  private

  !
  !>INIT ED SOLVER
  !
  interface ed_init_solver
     module procedure :: ed_init_solver_single
     module procedure :: ed_init_solver_lattice
#ifdef _MPI
     module procedure :: ed_init_solver_single_mpi
     module procedure :: ed_init_solver_lattice_mpi
#endif
  end interface ed_init_solver
  !>
  public :: ed_init_solver


  !
  !> ED SOLVER
  !
  interface ed_solve
     module procedure :: ed_solve_single
     module procedure :: ed_solve_lattice
#ifdef _MPI
     module procedure :: ed_solve_single_mpi
     module procedure :: ed_solve_lattice_mpi
#endif
  end interface ed_solve
  !>
  public :: ed_solve



  !
  !> ED REBUILD GF
  !
  interface ed_rebuild_gf
     module procedure :: ed_rebuild_gf_single
     module procedure :: ed_rebuild_gf_lattice
#ifdef _MPI
     module procedure :: ed_rebuild_gf_single_mpi
     module procedure :: ed_rebuild_gf_lattice_mpi
#endif
  end interface ed_rebuild_gf
  !>
  public :: ed_rebuild_gf



contains

  

  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: allocate and initialize one or multiple baths -+!
  !+-----------------------------------------------------------------------------+!
  !                              SINGLE SITE                                      !
  !+-----------------------------------------------------------------------------+!
  subroutine ed_init_solver_single(bath)
    real(8),dimension(:),intent(inout) :: bath
    logical                            :: check
    logical,save                       :: isetup=.true.
    integer                            :: i
    logical                            :: MPI_MASTER=.true.
    integer                            :: MPI_ERR
    !
    write(LOGfile,"(A)")"INIT SOLVER FOR "//trim(ed_file_suffix)
    !
    !Init ED Structure & memory
    if(isetup)call init_ed_structure()
    !
    check = check_bath_dimension(bath)
    if(.not.check)stop "init_ed_solver_single error: wrong bath dimensions"
    !
    bath = 0d0
    !
    call allocate_dmft_bath(dmft_bath)
    call init_dmft_bath(dmft_bath)
    call get_dmft_bath(dmft_bath,bath)
    !
    if(isetup)then
       call setup_global
    endif
    call deallocate_dmft_bath(dmft_bath)
    isetup=.false.
    !
  end subroutine ed_init_solver_single

#ifdef _MPI
  subroutine ed_init_solver_single_mpi(MpiComm,bath)
    integer                            :: MpiComm
    real(8),dimension(:),intent(inout) :: bath
    logical                            :: check 
    logical,save                       :: isetup=.true.
    integer                            :: i
    !
    !SET THE LOCAL MPI COMMUNICATOR :
    call ed_set_MpiComm(MpiComm)
    !
    write(LOGfile,"(A)")"INIT SOLVER FOR "//trim(ed_file_suffix)
    !
    !Init ED Structure & memory
    if(isetup)call init_ed_structure() 
    !
    check = check_bath_dimension(bath)
    if(.not.check)stop "init_ed_solver_single error: wrong bath dimensions"
    !
    bath = 0d0
    !
    call allocate_dmft_bath(dmft_bath)
    call init_dmft_bath(dmft_bath)
    call get_dmft_bath(dmft_bath,bath)
    !
    if(isetup)then
       call setup_global
    endif
    call deallocate_dmft_bath(dmft_bath)
    isetup=.false.
    !
    call ed_del_MpiComm()
    !
  end subroutine ed_init_solver_single_mpi
#endif


  !+-----------------------------------------------------------------------------+!
  !                           INEQUVALENT SITES                                   !
  !+-----------------------------------------------------------------------------+!
  subroutine ed_init_solver_lattice(bath)
    real(8),dimension(:,:) :: bath        ![Nlat,Nb]
    integer                :: ilat,Nineq
    logical                :: check_dim
    integer                :: MPI_ERR
    !
    if(allocated(dens_ineq))deallocate(dens_ineq)
    if(allocated(docc_ineq))deallocate(docc_ineq)
    if(allocated(mag_ineq))deallocate(mag_ineq)
    if(allocated(phisc_ineq))deallocate(phisc_ineq)
    if(allocated(e_ineq))deallocate(e_ineq)
    if(allocated(dd_ineq))deallocate(dd_ineq)
    if(allocated(Smats_ineq))deallocate(Smats_ineq)
    if(allocated(Sreal_ineq))deallocate(Sreal_ineq)
    if(allocated(SAmats_ineq))deallocate(SAmats_ineq)
    if(allocated(SAreal_ineq))deallocate(SAreal_ineq)
    if(allocated(Gmats_ineq))deallocate(Gmats_ineq)
    if(allocated(Greal_ineq))deallocate(Greal_ineq)
    if(allocated(Fmats_ineq))deallocate(Fmats_ineq)
    if(allocated(Freal_ineq))deallocate(Freal_ineq)
    if(allocated(Dmats_ph_ineq))deallocate(Dmats_ph_ineq)
    if(allocated(Dreal_ph_ineq))deallocate(Dreal_ph_ineq)
    if(allocated(neigen_sector_ineq))deallocate(neigen_sector_ineq)
    if(allocated(neigen_total_ineq))deallocate(neigen_total_ineq)
    !
    Nineq = size(bath,1)
    if(bath_type=='replica' .AND. .not.allocated(Hreplica_lambda_ineq))&
         stop "ERROR ed_init_solver: replica parameters lambda not defined for all sites"
    !
    allocate(dens_ineq(Nineq,Norb))
    allocate(docc_ineq(Nineq,Norb))
    allocate(mag_ineq(Nineq,3,Norb))
    allocate(phisc_ineq(Nineq,Norb))
    allocate(e_ineq(Nineq,4))
    allocate(dd_ineq(Nineq,4))
    !
    allocate(Smats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Sreal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(SAmats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(SAreal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
    !
    allocate(Gmats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Greal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(Fmats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Freal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
    !
    allocate(Dmats_ph_ineq(Nineq,Lmats))
    allocate(Dreal_ph_ineq(Nineq,Lreal))    
    !
    allocate(imp_density_matrix_ineq(Nineq,Nspin,Nspin,Norb,Norb))
    !
    do ilat=1,Nineq
       call ed_set_suffix(ilat)
       if(bath_type=='replica')call Hreplica_site(ilat)
       !set the ilat-th lambda vector basis for the replica bath
       call ed_init_solver_single(bath(ilat,:))
    enddo
    call ed_reset_suffix
    !
    allocate(neigen_sector_ineq(Nineq,Nsectors))
    allocate(neigen_total_ineq(Nineq))
    do ilat=1,Nineq
       neigen_sector_ineq(ilat,:) = neigen_sector(:)
       neigen_total_ineq(ilat)    = lanc_nstates_total
    end do
    !
  end subroutine ed_init_solver_lattice

#ifdef _MPI
  subroutine ed_init_solver_lattice_mpi(MpiComm,bath)
    integer                :: MpiComm
    real(8),dimension(:,:) :: bath ![Nlat][:]
    integer                :: ilat,Nineq,Nsect
    logical                :: check_dim
    character(len=5)       :: tmp_suffix
    integer                :: MPI_ERR
    !
    if(allocated(dens_ineq))deallocate(dens_ineq)
    if(allocated(docc_ineq))deallocate(docc_ineq)
    if(allocated(mag_ineq))deallocate(mag_ineq)
    if(allocated(phisc_ineq))deallocate(phisc_ineq)
    if(allocated(e_ineq))deallocate(e_ineq)
    if(allocated(dd_ineq))deallocate(dd_ineq)
    if(allocated(Smats_ineq))deallocate(Smats_ineq)
    if(allocated(Sreal_ineq))deallocate(Sreal_ineq)
    if(allocated(SAmats_ineq))deallocate(SAmats_ineq)
    if(allocated(SAreal_ineq))deallocate(SAreal_ineq)
    if(allocated(Gmats_ineq))deallocate(Gmats_ineq)
    if(allocated(Greal_ineq))deallocate(Greal_ineq)
    if(allocated(Fmats_ineq))deallocate(Fmats_ineq)
    if(allocated(Freal_ineq))deallocate(Freal_ineq)
    if(allocated(Dmats_ph_ineq))deallocate(Dmats_ph_ineq)
    if(allocated(Dreal_ph_ineq))deallocate(Dreal_ph_ineq)
    if(allocated(imp_density_matrix_ineq))deallocate(imp_density_matrix_ineq)
    if(allocated(neigen_sector_ineq))deallocate(neigen_sector_ineq)
    if(allocated(neigen_total_ineq))deallocate(neigen_total_ineq)
    !
    Nineq = size(bath,1)
    if(bath_type=='replica'.AND..not.allocated(Hreplica_lambda_ineq))&
         stop "ERROR ed_init_solver: replica parameters lambda not defined for all sites"
    !
    allocate(dens_ineq(Nineq,Norb))
    allocate(docc_ineq(Nineq,Norb))
    allocate(mag_ineq(Nineq,3,Norb))
    allocate(phisc_ineq(Nineq,Norb))
    allocate(e_ineq(Nineq,4))
    allocate(dd_ineq(Nineq,4))
    !
    allocate(Smats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Sreal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(SAmats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(SAreal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
    !
    allocate(Gmats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Greal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(Fmats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Freal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
    !
    allocate(Dmats_ph_ineq(Nineq,Lmats))
    allocate(Dreal_ph_ineq(Nineq,Lreal))
    !
    allocate(imp_density_matrix_ineq(Nineq,Nspin,Nspin,Norb,Norb))
    !
    do ilat=1,Nineq
       call ed_set_suffix(ilat)
       if(bath_type=='replica')call Hreplica_site(ilat)
       call ed_init_solver_single(bath(ilat,:))
    enddo
    call MPI_Barrier(MpiComm,MPI_ERR)
    call ed_reset_suffix
    !
    !This follows because Nsectors is defined after ED is initialized
    allocate(neigen_sector_ineq(Nineq,Nsectors))
    allocate(neigen_total_ineq(Nineq))
    do ilat=1,Nineq       
       neigen_sector_ineq(ilat,:) = neigen_sector(:)
       neigen_total_ineq(ilat)    = lanc_nstates_total
    end do
    !
  end subroutine ed_init_solver_lattice_mpi
#endif







  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################








  !+-----------------------------------------------------------------------------+!
  !PURPOSE: solve the impurity problems for a single or many independent
  ! lattice site using ED. 
  !+-----------------------------------------------------------------------------+!
  !+-----------------------------------------------------------------------------+!
  !                              SINGLE SITE                                      !
  !+-----------------------------------------------------------------------------+!
  subroutine ed_solve_single(bath,Hloc,sflag)
    real(8),dimension(:),intent(in) :: bath
    complex(8),intent(in)           :: Hloc(Nspin,Nspin,Norb,Norb)
    logical,optional                :: sflag
    logical                         :: check,iflag
    !
    iflag=.true. ;if(present(sflag))iflag=sflag
    !
    if(MpiMaster)call save_input_file(str(ed_input_file))
    !
    call set_Himpurity(Hloc)
    !
    check = check_bath_dimension(bath)
    if(.not.check)stop "ED_SOLVE_SINGLE Error: wrong bath dimensions"
    !
    call allocate_dmft_bath(dmft_bath)
    call set_dmft_bath(bath,dmft_bath)
    call write_dmft_bath(dmft_bath,LOGfile)
    if(MpiMaster)call save_dmft_bath(dmft_bath,used=.true.)
    !
    !
    !SOLVE THE QUANTUM IMPURITY PROBLEM:
    call diagonalize_impurity()
    if(iflag)then
       call buildgf_impurity()
       call buildchi_impurity()
    endif
    call observables_impurity()
    call local_energy_impurity()
    !
    call deallocate_dmft_bath(dmft_bath)   
    call es_delete_espace(state_list)
    !
    nullify(spHtimesV_p)
    nullify(spHtimesV_cc)
    write(Logfile,"(A)")""
  end subroutine ed_solve_single



#ifdef _MPI
  !+-----------------------------------------------------------------------------+!
  !                              SINGLE SITE                                      !
  !+-----------------------------------------------------------------------------+!
  subroutine ed_solve_single_mpi(MpiComm,bath,Hloc,sflag)
    integer                         :: MpiComm
    real(8),dimension(:),intent(in) :: bath
    complex(8),intent(in)           :: Hloc(Nspin,Nspin,Norb,Norb)
    logical,optional                :: sflag
    logical                         :: check,iflag
    !
    iflag=.true. ;if(present(sflag))iflag=sflag
    !
    !SET THE LOCAL MPI COMMUNICATOR :
    call ed_set_MpiComm(MpiComm)
    !
    if(MpiMaster)call save_input_file(str(ed_input_file))
    !
    call set_Himpurity(Hloc)
    !
    check   = check_bath_dimension(bath)
    if(.not.check)stop "ED_SOLVE_SINGLE Error: wrong bath dimensions"
    !
    call allocate_dmft_bath(dmft_bath)
    call set_dmft_bath(bath,dmft_bath)
    call write_dmft_bath(dmft_bath,LOGfile)
    if(MpiMaster)call save_dmft_bath(dmft_bath,used=.true.)
    !
    !
    !SOLVE THE QUANTUM IMPURITY PROBLEM:
    call diagonalize_impurity()
    if(iflag)then
       call buildgf_impurity()
       call buildchi_impurity()
    endif
    call observables_impurity()
    call local_energy_impurity()
    !
    call deallocate_dmft_bath(dmft_bath)
    call es_delete_espace(state_list)
    !
    !DELETE THE LOCAL MPI COMMUNICATOR:
    call ed_del_MpiComm()
    !
    nullify(spHtimesV_cc)
    nullify(spHtimesV_p)
    write(Logfile,"(A)")""
  end subroutine ed_solve_single_mpi
#endif













  !+-----------------------------------------------------------------------------+!
  !                          INEQUIVALENT SITES                                   !
  !+-----------------------------------------------------------------------------+!
  subroutine ed_solve_lattice(bath,Hloc,Uloc_ii,Ust_ii,Jh_ii,Jp_ii,Jx_ii)
    !inputs
    real(8)          :: bath(:,:) ![Nlat][Nb]
    complex(8)       :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    real(8),optional :: Uloc_ii(size(bath,1),Norb)
    real(8),optional :: Ust_ii(size(bath,1))
    real(8),optional :: Jh_ii(size(bath,1))
    real(8),optional :: Jp_ii(size(bath,1))
    real(8),optional :: Jx_ii(size(bath,1))
    ! 
    integer          :: i,j,ilat,iorb,jorb,ispin,jspin
    integer          :: Nineq
    logical          :: check_dim
    character(len=5) :: tmp_suffix
    !
    logical          :: MPI_MASTER=.true.
    !
    ! Check dimensions !
    Nineq=size(bath,1)
    !
    if(size(neigen_sector_ineq,1)<Nineq)stop "ed_solve_lattice error: size(neigen_sector_ineq,1)<Nineq"
    if(size(neigen_total_ineq)<Nineq)stop "ed_solve_lattice error: size(neigen_total_ineq,1)<Nineq"
    !
    !Check the dimensions of the bath are ok:
    do ilat=1,Nineq
       check_dim = check_bath_dimension(bath(ilat,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    !
    Smats_ineq    = zero ; Sreal_ineq    = zero ; SAmats_ineq   = zero ; SAreal_ineq   = zero 
    Gmats_ineq    = zero ; Greal_ineq    = zero ; Fmats_ineq    = zero ; Freal_ineq    = zero
    Dmats_ph_ineq = zero ; Dreal_ph_ineq = zero 
    dens_ineq     = 0d0  ; docc_ineq     = 0d0  
    mag_ineq      = 0d0  ; phisc_ineq    = 0d0  
    e_ineq        = 0d0  ; dd_ineq       = 0d0 
    imp_density_matrix_ineq = zero
    !
    call start_timer
    !
    do ilat = 1, Nineq
       write(LOGfile,*)" SOLVING INEQ SITE: "//str(ilat,Npad=4)
       !
       ed_file_suffix=reg(ineq_site_suffix)//str(ilat,site_indx_padding)
       !
       !If required set the local value of U per each site
       if(present(Uloc_ii))Uloc(1:Norb) = Uloc_ii(ilat,1:Norb)
       if(present(Ust_ii)) Ust = Ust_ii(ilat)
       if(present(Jh_ii))  Jh  = Jh_ii(ilat)
       if(present(Jp_ii))  Jp  = Jp_ii(ilat)
       if(present(Jx_ii))  Jx  = Jx_ii(ilat)
       !
       !Solve the impurity problem for the ilat-th site
       neigen_sector(:)   = neigen_sector_ineq(ilat,:)
       lanc_nstates_total = neigen_total_ineq(ilat)
       !
       call ed_solve_single(bath(ilat,:),Hloc(ilat,:,:,:,:)) !Hloc(ilat,...)-->impHloc
       !
       neigen_sector_ineq(ilat,:)  = neigen_sector
       neigen_total_ineq(ilat)     = lanc_nstates_total
       Smats_ineq(ilat,:,:,:,:,:)  = impSmats
       Sreal_ineq(ilat,:,:,:,:,:)  = impSreal
       SAmats_ineq(ilat,:,:,:,:,:) = impSAmats
       SAreal_ineq(ilat,:,:,:,:,:) = impSAreal
       Gmats_ineq(ilat,:,:,:,:,:)  = impGmats
       Greal_ineq(ilat,:,:,:,:,:)  = impGreal
       Fmats_ineq(ilat,:,:,:,:,:)  = impFmats
       Freal_ineq(ilat,:,:,:,:,:)  = impFreal
       Dmats_ph_ineq(ilat,:)       = impDmats_ph
       Dreal_ph_ineq(ilat,:)       = impDreal_ph
       dens_ineq(ilat,1:Norb)      = ed_dens(1:Norb)
       docc_ineq(ilat,1:Norb)      = ed_docc(1:Norb)
       mag_ineq(ilat,:,1:Norb)     = ed_mag(:,1:Norb)
       phisc_ineq(ilat,1:Norb)     = ed_phisc(1:Norb)
       e_ineq(ilat,:)              = [ed_Epot,ed_Eint,ed_Ehartree,ed_Eknot]
       dd_ineq(ilat,:)             = [ed_Dust,ed_Dund,ed_Dse,ed_Dph]
       imp_density_matrix_ineq(ilat,:,:,:,:) = imp_density_matrix(:,:,:,:)
       !
    enddo
    !
    call stop_timer(unit=LOGfile)
    !
    call ed_reset_suffix
    !
  end subroutine ed_solve_lattice


  !FALL BACK: DO A VERSION THAT DOES THE SITES IN PARALLEL USING SERIAL ED CODE
#ifdef _MPI
  subroutine ed_solve_lattice_mpi(MpiComm,bath,Hloc,mpi_lanc,Uloc_ii,Ust_ii,Jh_ii,Jp_ii,Jx_ii)
    integer          :: MpiComm
    !inputs
    real(8)          :: bath(:,:) ![Nlat][Nb]
    complex(8)       :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    logical,optional :: mpi_lanc
    real(8),optional :: Uloc_ii(size(bath,1),Norb)
    real(8),optional :: Ust_ii(size(bath,1))
    real(8),optional :: Jh_ii(size(bath,1))
    real(8),optional :: Jp_ii(size(bath,1))
    real(8),optional :: Jx_ii(size(bath,1))
    !MPI  auxiliary vars
    complex(8)       :: Smats_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)       :: Sreal_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)       :: SAmats_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)       :: SAreal_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)       :: Gmats_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)       :: Greal_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)       :: Fmats_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)       :: Freal_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)       :: Dmats_tmp(size(bath,1),Lmats)
    complex(8)       :: Dreal_tmp(size(bath,1),Lreal)
    real(8)          :: dens_tmp(size(bath,1),Norb)
    real(8)          :: docc_tmp(size(bath,1),Norb)
    real(8)          :: mag_tmp(size(bath,1),3,Norb)
    real(8)          :: phisc_tmp(size(bath,1),Norb)
    real(8)          :: e_tmp(size(bath,1),4)
    real(8)          :: dd_tmp(size(bath,1),4)
    !    
    complex(8)       :: imp_density_matrix_tmp(size(bath,1),Nspin,Nspin,Norb,Norb)
    !
    integer          :: neigen_sectortmp(size(bath,1),Nsectors)
    integer          :: neigen_totaltmp(size(bath,1))
    ! 
    integer          :: i,j,ilat,iorb,jorb,ispin,jspin
    integer          :: Nineq
    logical          :: check_dim,mpi_lanc_
    character(len=5) :: tmp_suffix
    !
    integer          :: MPI_ID=0
    integer          :: MPI_SIZE=1
    logical          :: MPI_MASTER=.true.
    !
    integer          :: mpi_err 
    !
    MPI_ID     = get_Rank_MPI(MpiComm)
    MPI_SIZE   = get_Size_MPI(MpiComm)
    MPI_MASTER = get_Master_MPI(MpiComm)
    !
    mpi_lanc_=.false.;if(present(mpi_lanc))mpi_lanc_=mpi_lanc
    !
    ! Check dimensions !
    Nineq=size(bath,1)
    !
    if(size(neigen_sector_ineq,1)<Nineq)stop "ed_solve_lattice error: size(neigen_sectorii,1)<Nineq"
    if(size(neigen_total_ineq)<Nineq)stop "ed_solve_lattice error: size(neigen_totalii,1)<Nineq"
    !
    !Check the dimensions of the bath are ok.
    !This can always be done in parallel no issues with mpi_lanc
    do ilat=1+MPI_ID,Nineq,MPI_SIZE
       check_dim = check_bath_dimension(bath(ilat,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    !
    Smats_ineq    = zero ; Sreal_ineq    = zero ; SAmats_ineq   = zero ; SAreal_ineq   = zero 
    Gmats_ineq    = zero ; Greal_ineq    = zero ; Fmats_ineq    = zero ; Freal_ineq    = zero
    Dmats_ph_ineq = zero ; Dreal_ph_ineq = zero 
    dens_ineq     = 0d0  ; docc_ineq     = 0d0
    mag_ineq      = 0d0  ; phisc_ineq    = 0d0  
    e_ineq        = 0d0  ; dd_ineq       = 0d0 
    imp_density_matrix_ineq = zero
    !
    Smats_tmp  = zero ; Sreal_tmp  = zero ; SAmats_tmp = zero ; SAreal_tmp = zero
    Gmats_tmp  = zero ; Greal_tmp  = zero ; Fmats_tmp  = zero ; Freal_tmp  = zero
    Dmats_tmp  = zero ; Dreal_tmp  = zero
    dens_tmp   = 0d0  ; docc_tmp   = 0d0
    mag_tmp    = 0d0  ; phisc_tmp  = 0d0
    e_tmp      = 0d0  ; dd_tmp     = 0d0
    neigen_sectortmp = 0
    neigen_totaltmp  = 0
    imp_density_matrix_tmp = zero
    !
    select case(mpi_lanc_)
    case default              !mpi_lanc=False => solve sites with MPI
       if(MPI_MASTER)call start_timer
       do ilat = 1 + MPI_ID, Nineq, MPI_SIZE
          write(LOGfile,*)str(MPI_ID)//" SOLVES INEQ SITE: "//str(ilat,Npad=4)
          !
          ed_file_suffix=reg(ineq_site_suffix)//str(ilat,site_indx_padding)
          !
          !If required set the local value of U per each site
          if(present(Uloc_ii))Uloc(1:Norb) = Uloc_ii(ilat,1:Norb)
          if(present(Ust_ii)) Ust = Ust_ii(ilat)
          if(present(Jh_ii))  Jh  = Jh_ii(ilat)
          if(present(Jp_ii))  Jp  = Jp_ii(ilat)
          if(present(Jx_ii))  Jx  = Jx_ii(ilat)
          !
          !Solve the impurity problem for the ilat-th site, this are set at init time
          neigen_sector(:)   = neigen_sector_ineq(ilat,:)
          lanc_nstates_total = neigen_total_ineq(ilat)
          !
          call ed_solve_single(bath(ilat,:),Hloc(ilat,:,:,:,:))
          !
          neigen_sectortmp(ilat,:)   = neigen_sector(:)
          neigen_totaltmp(ilat)      = lanc_nstates_total
          Smats_tmp(ilat,:,:,:,:,:)  = impSmats
          Sreal_tmp(ilat,:,:,:,:,:)  = impSreal
          SAmats_tmp(ilat,:,:,:,:,:) = impSAmats
          SAreal_tmp(ilat,:,:,:,:,:) = impSAreal
          Gmats_tmp(ilat,:,:,:,:,:)  = impGmats
          Greal_tmp(ilat,:,:,:,:,:)  = impGreal
          Fmats_tmp(ilat,:,:,:,:,:)  = impFmats
          Freal_tmp(ilat,:,:,:,:,:)  = impFreal
          Dmats_tmp(ilat,:)          = impDmats_ph
          Dreal_tmp(ilat,:)          = impDreal_ph
          dens_tmp(ilat,1:Norb)      = ed_dens(1:Norb)
          docc_tmp(ilat,1:Norb)      = ed_docc(1:Norb)
          mag_tmp(ilat,:,1:Norb)     = ed_mag(:,1:Norb)
          phisc_tmp(ilat,1:Norb)     = ed_phisc(1:Norb)
          e_tmp(ilat,:)              = [ed_Epot,ed_Eint,ed_Ehartree,ed_Eknot]
          dd_tmp(ilat,:)             = [ed_Dust,ed_Dund,ed_Dse,ed_Dph]
          imp_density_matrix_tmp(ilat,:,:,:,:) = imp_density_matrix(:,:,:,:)
       enddo
       call MPI_Barrier(MpiComm,MPI_ERR)
       if(MPI_MASTER)call stop_timer(unit=LOGfile)
       call ed_reset_suffix
       !
       call AllReduce_MPI(MpiComm,Smats_tmp,Smats_ineq)
       call AllReduce_MPI(MpiComm,Sreal_tmp,Sreal_ineq)
       call AllReduce_MPI(MpiComm,SAmats_tmp,SAmats_ineq)
       call AllReduce_MPI(MpiComm,SAreal_tmp,SAreal_ineq)
       call AllReduce_MPI(MpiComm,Gmats_tmp,Gmats_ineq)
       call AllReduce_MPI(MpiComm,Greal_tmp,Greal_ineq)
       call AllReduce_MPI(MpiComm,Fmats_tmp,Fmats_ineq)
       call AllReduce_MPI(MpiComm,Freal_tmp,Freal_ineq)
       call AllReduce_MPI(MpiComm,Dmats_tmp,Dmats_ph_ineq)
       call AllReduce_MPI(MpiComm,Dreal_tmp,Dreal_ph_ineq)
       call AllReduce_MPI(MpiComm,dens_tmp,dens_ineq)
       call AllReduce_MPI(MpiComm,docc_tmp,docc_ineq)
       call AllReduce_MPI(MpiComm,mag_tmp,mag_ineq)
       call AllReduce_MPI(MpiComm,phisc_tmp,phisc_ineq)
       call AllReduce_MPI(MpiComm,e_tmp,e_ineq)
       call AllReduce_MPI(MpiComm,dd_tmp,dd_ineq)
       call AllReduce_MPI(MpiComm,imp_density_matrix_tmp,imp_density_matrix_ineq)
       neigen_sector_ineq=0
       neigen_total_ineq=0
       call AllReduce_MPI(MpiComm,neigen_sectortmp,neigen_sector_ineq)
       call AllReduce_MPI(MpiComm,neigen_totaltmp,neigen_total_ineq)
       !
       !       
    case(.true.)                !solve sites serial, Lanczos with MPI
       if(MPI_MASTER)call start_timer
       do ilat = 1, Nineq
          write(LOGfile,*)" SOLVES INEQ SITE: "//str(ilat,Npad=4)
          ed_file_suffix=reg(ineq_site_suffix)//str(ilat,site_indx_padding)
          !
          !If required set the local value of U per each site
          if(present(Uloc_ii))Uloc(1:Norb) = Uloc_ii(ilat,1:Norb)
          if(present(Ust_ii)) Ust = Ust_ii(ilat)
          if(present(Jh_ii))  Jh  = Jh_ii(ilat)
          if(present(Jp_ii))  Jp  = Jp_ii(ilat)
          if(present(Jx_ii))  Jx  = Jx_ii(ilat)
          !
          !Solve the impurity problem for the ilat-th site
          neigen_sector(:)   = neigen_sector_ineq(ilat,:)
          lanc_nstates_total = neigen_total_ineq(ilat)
          !
          call ed_solve_single_mpi(MpiComm,bath(ilat,:),Hloc(ilat,:,:,:,:))
          !
          neigen_sector_ineq(ilat,:)  = neigen_sector(:)
          neigen_total_ineq(ilat)     = lanc_nstates_total
          Smats_ineq(ilat,:,:,:,:,:)  = impSmats
          Sreal_ineq(ilat,:,:,:,:,:)  = impSreal
          SAmats_ineq(ilat,:,:,:,:,:) = impSAmats
          SAreal_ineq(ilat,:,:,:,:,:) = impSAreal
          Gmats_ineq(ilat,:,:,:,:,:)  = impGmats
          Greal_ineq(ilat,:,:,:,:,:)  = impGreal
          Fmats_ineq(ilat,:,:,:,:,:)  = impFmats
          Freal_ineq(ilat,:,:,:,:,:)  = impFreal
          Dmats_ph_ineq(ilat,:)       = impDmats_ph
          Dreal_ph_ineq(ilat,:)       = impDreal_ph
          dens_ineq(ilat,1:Norb)      = ed_dens(1:Norb)
          docc_ineq(ilat,1:Norb)      = ed_docc(1:Norb)
          mag_ineq(ilat,:,1:Norb)     = ed_mag(:,1:Norb)
          phisc_ineq(ilat,1:Norb)     = ed_phisc(1:Norb)
          e_ineq(ilat,:)              = [ed_Epot,ed_Eint,ed_Ehartree,ed_Eknot]
          dd_ineq(ilat,:)             = [ed_Dust,ed_Dund,ed_Dse,ed_Dph]
          imp_density_matrix_ineq(ilat,:,:,:,:) = imp_density_matrix(:,:,:,:)
       enddo
       if(MPI_MASTER)call stop_timer(unit=LOGfile)
       call ed_reset_suffix
    end select
    !
  end subroutine ed_solve_lattice_mpi
#endif






  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################




  


  !+-----------------------------------------------------------------------------+!
  !PURPOSE: rebuild impurity GF for a single or many sites
  !+-----------------------------------------------------------------------------+!
  !+-----------------------------------------------------------------------------+!
  !                              SINGLE SITE                                      !
  !+-----------------------------------------------------------------------------+!
  subroutine ed_rebuild_gf_single(Hloc)
    complex(8),intent(in)           :: Hloc(Nspin,Nspin,Norb,Norb)
    logical                         :: check
    !
    !
    if(MpiMaster)call save_input_file(str(ed_input_file))
    !
    call set_Himpurity(Hloc)
    !
    call allocate_dmft_bath(dmft_bath)
    call init_dmft_bath(dmft_bath,used=.true.)
    if(MpiMaster)call write_dmft_bath(dmft_bath,LOGfile)
    !
    !
    call rebuildGF_impurity()
    !
    !
    call deallocate_dmft_bath(dmft_bath)
    !
    nullify(spHtimesV_p)
    nullify(spHtimesV_cc)
  end subroutine ed_rebuild_gf_single



#ifdef _MPI
  subroutine ed_rebuild_gf_single_mpi(MpiComm,Hloc)
    integer                         :: MpiComm
    complex(8),intent(in)           :: Hloc(Nspin,Nspin,Norb,Norb)
    logical                         :: check
    !
    !
    !SET THE LOCAL MPI COMMUNICATOR :
    call ed_set_MpiComm(MpiComm)
    !
    if(MpiMaster)call save_input_file(str(ed_input_file))
    !
    call set_Himpurity(Hloc)
    !
    call allocate_dmft_bath(dmft_bath)
    call init_dmft_bath(dmft_bath,used=.true.)
    if(MpiMaster)call write_dmft_bath(dmft_bath,LOGfile)
    !
    !
    call rebuildGF_impurity()
    !
    !
    call deallocate_dmft_bath(dmft_bath)
    !
    !DELETE THE LOCAL MPI COMMUNICATOR:
    call ed_del_MpiComm()
    !
    nullify(spHtimesV_cc)
    nullify(spHtimesV_p)
  end subroutine ed_rebuild_gf_single_mpi
#endif




  !+-----------------------------------------------------------------------------+!
  !                          INEQUIVALENT SITES                                   !
  !+-----------------------------------------------------------------------------+!
  subroutine ed_rebuild_gf_lattice(Hloc)
    complex(8)              :: Hloc(:,:,:,:,:) ![Nineq,Nspin,Nspin,Norb,Norb]
    ! 
    integer                 :: i,j,ilat,iorb,jorb,ispin,jspin
    integer                 :: Nineq
    logical                 :: check_dim
    character(len=5)        :: tmp_suffix
    !
    logical                 :: MPI_MASTER=.true.
    !
    !
    ! Check dimensions !
    Nineq=size(Hloc,1)
    call assert_shape(Hloc,[Nineq,Nspin,Nspin,Norb,Norb],"ed_rebuild_gf_lattice","hloc")
    !
    !Check the dimensions of the bath are ok:
    if(allocated(Smats_ineq))deallocate(Smats_ineq)
    if(allocated(SAmats_ineq))deallocate(SAmats_ineq)
    if(allocated(Gmats_ineq))deallocate(Gmats_ineq)
    if(allocated(Fmats_ineq))deallocate(Fmats_ineq)
    allocate(Smats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(SAmats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Gmats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Fmats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
    Smats_ineq    = zero ; SAmats_ineq   = zero 
    Gmats_ineq    = zero ; Fmats_ineq    = zero
    !
    if(allocated(Sreal_ineq))deallocate(Sreal_ineq)
    if(allocated(SAreal_ineq))deallocate(SAreal_ineq)
    if(allocated(Greal_ineq))deallocate(Greal_ineq)
    if(allocated(Freal_ineq))deallocate(Freal_ineq)
    allocate(Sreal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(SAreal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(Greal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(Freal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
    Sreal_ineq    = zero ; SAreal_ineq   = zero 
    Greal_ineq    = zero ; Freal_ineq    = zero
    !
    call start_timer
    do ilat = 1, Nineq
       write(LOGfile,*)" SOLVING INEQ SITE: "//str(ilat,Npad=4)
       !
       call ed_set_suffix(ilat)
       call ed_rebuild_gf_single(Hloc(ilat,:,:,:,:))
       !
       Smats_ineq(ilat,:,:,:,:,:)  = impSmats
       SAmats_ineq(ilat,:,:,:,:,:) = impSAmats
       Gmats_ineq(ilat,:,:,:,:,:)  = impGmats
       Fmats_ineq(ilat,:,:,:,:,:)  = impFmats
       Sreal_ineq(ilat,:,:,:,:,:)  = impSreal
       SAreal_ineq(ilat,:,:,:,:,:) = impSAreal
       Greal_ineq(ilat,:,:,:,:,:)  = impGreal
       Freal_ineq(ilat,:,:,:,:,:)  = impFreal
       !
    enddo
    call stop_timer(unit=LOGfile)
    !
    call ed_reset_suffix
    !
  end subroutine ed_rebuild_gf_lattice


  !FALL BACK: DO A VERSION THAT DOES THE SITES IN PARALLEL USING SERIAL ED CODE
#ifdef _MPI
  subroutine ed_rebuild_gf_lattice_mpi(MpiComm,Hloc)
    integer                                       :: MpiComm
    !inputs
    complex(8)                                    :: Hloc(:,:,:,:,:) ![Nineq,Nspin,Nspin,Norb,Norb]
    !MPI  auxiliary vars
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Smats_tmp,Sreal_tmp
    complex(8),dimension(:,:,:,:,:,:),allocatable :: SAmats_tmp,SAreal_tmp
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Gmats_tmp,Greal_tmp
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Fmats_tmp,Freal_tmp
    ! 
    integer                                       :: i,j,ilat,iorb,jorb,ispin,jspin
    integer                                       :: Nineq
    logical                                       :: check_dim
    character(len=5)                              :: tmp_suffix
    !
    integer                                       :: MPI_ID=0
    integer                                       :: MPI_SIZE=1
    logical                                       :: MPI_MASTER=.true.
    !
    integer                                       :: mpi_err 
    !
    MPI_ID     = get_Rank_MPI(MpiComm)
    MPI_SIZE   = get_Size_MPI(MpiComm)
    MPI_MASTER = get_Master_MPI(MpiComm)
    !
    !
    ! Check dimensions !
    Nineq=size(hloc,1)
    call assert_shape(Hloc,[Nineq,Nspin,Nspin,Norb,Norb],"ed_rebuild_gf_lattice_mpi","hloc")
    !
    if(allocated(Smats_ineq))deallocate(Smats_ineq)
    if(allocated(SAmats_ineq))deallocate(SAmats_ineq)
    if(allocated(Gmats_ineq))deallocate(Gmats_ineq)
    if(allocated(Fmats_ineq))deallocate(Fmats_ineq)
    allocate(Smats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(SAmats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Gmats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Fmats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
    Smats_ineq    = zero ; SAmats_ineq   = zero 
    Gmats_ineq    = zero ; Fmats_ineq    = zero
    if(allocated(Sreal_ineq))deallocate(Sreal_ineq)
    if(allocated(SAreal_ineq))deallocate(SAreal_ineq)
    if(allocated(Greal_ineq))deallocate(Greal_ineq)
    if(allocated(Freal_ineq))deallocate(Freal_ineq)
    allocate(Sreal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(SAreal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(Greal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(Freal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
    Sreal_ineq    = zero ; SAreal_ineq   = zero 
    Greal_ineq    = zero ; Freal_ineq    = zero
    !
    allocate(Smats_tmp(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Sreal_tmp(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(SAmats_tmp(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(SAreal_tmp(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(Gmats_tmp(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Greal_tmp(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(Fmats_tmp(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Freal_tmp(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
    Smats_tmp  = zero ; Sreal_tmp  = zero ; SAmats_tmp = zero ; SAreal_tmp = zero
    Gmats_tmp  = zero ; Greal_tmp  = zero ; Fmats_tmp  = zero ; Freal_tmp  = zero
    !
    if(MPI_MASTER)call start_timer
    do ilat = 1 + MPI_ID, Nineq, MPI_SIZE
       write(LOGfile,*)str(MPI_ID)//" SOLVES INEQ SITE: "//str(ilat,Npad=4)
       !
       call ed_set_suffix(ilat)
       call ed_rebuild_gf_single(Hloc(ilat,:,:,:,:))
       !
       Smats_tmp(ilat,:,:,:,:,:)  = impSmats
       SAmats_tmp(ilat,:,:,:,:,:) = impSAmats
       Gmats_tmp(ilat,:,:,:,:,:)  = impGmats
       Fmats_tmp(ilat,:,:,:,:,:)  = impFmats
       Sreal_tmp(ilat,:,:,:,:,:)  = impSreal
       SAreal_tmp(ilat,:,:,:,:,:) = impSAreal
       Greal_tmp(ilat,:,:,:,:,:)  = impGreal
       Freal_tmp(ilat,:,:,:,:,:)  = impFreal
    enddo
    call MPI_Barrier(MpiComm,MPI_ERR)
    if(MPI_MASTER)call stop_timer(unit=LOGfile)
    call ed_reset_suffix
    !
    call AllReduce_MPI(MpiComm,Smats_tmp,Smats_ineq)
    call AllReduce_MPI(MpiComm,SAmats_tmp,SAmats_ineq)
    call AllReduce_MPI(MpiComm,Gmats_tmp,Gmats_ineq)
    call AllReduce_MPI(MpiComm,Fmats_tmp,Fmats_ineq)
    call AllReduce_MPI(MpiComm,Sreal_tmp,Sreal_ineq)
    call AllReduce_MPI(MpiComm,SAreal_tmp,SAreal_ineq)
    call AllReduce_MPI(MpiComm,Greal_tmp,Greal_ineq)
    call AllReduce_MPI(MpiComm,Freal_tmp,Freal_ineq)
    !       
  end subroutine ed_rebuild_gf_lattice_mpi
#endif




end module ED_MAIN




