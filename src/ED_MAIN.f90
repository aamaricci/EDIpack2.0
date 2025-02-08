module ED_MAIN
  !Contains routine that initialize, run and finalize the Impurity model solver
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
  USE ED_RDM
  USE ED_DIAG
  USE ED_IO
  implicit none
  private

  !>INIT ED SOLVER
  interface ed_init_solver
     !
     !Initialize the Exact Diagonalization solver of `EDIpack2.0`. This procedure reserves and allocates all the  
     !memory required by the solver, performs all the consistency check and initializes the bath instance guessing or reading from a file.      
     !It requires as an input a double precision array of rank-1 [ :f:var:`nb` ] for the single-impurity case or
     !or rank-2 [ :f:var:`nb` , :f:var:`nlat` ] for the Real space DMFT case. :f:var:`nlat` is the number of inequivalent impurity sites,
     !while :f:var:`nb` depends on the bath size and geometry and can be obtained from :f:func:`get_bath_dimension` .
     !
     module procedure :: ed_init_solver_single
     module procedure :: ed_init_solver_lattice
  end interface ed_init_solver


  !> ED SOLVER
  interface ed_solve
     !
     !Launch the Exact Diagonalizaton solver for the single-site and multiple-site (R-DMFT) quantum impurity problem.
     !It requires as an input a double precision array of rank-1 [ :f:var:`nb` ] for the single-impurity case or
     !or rank-2 [ :f:var:`nb` , :f:var:`nlat` ] for the Real space DMFT case. :f:var:`nlat` is the number of inequivalent impurity sites,
     !while :f:var:`nb` depends on the bath size and geometry and can be obtained from :f:func:`get_bath_dimension` .
     !
     ! The solution is achieved in this sequence:
     !
     !  #. setup the MPI environment, if any 
     !  #. Set the internal bath instance :f:var:`dmft_bath` copying from the user provided input :f:var:`bath`
     !  #. Get the low energy spectrum: call :f:func:`diagonalize_impurity`
     !  #. Get the impurity Green's functions: call :f:func:`buildgf_impurity` (if :f:var:`sflag` = :code:`.true.` )
     !  #. Get the impurity susceptibilities, if any: call :f:func:`buildchi_impurity` (if :f:var:`sflag` = :code:`.true.` )
     !  #. Get the impurity observables: call :f:func:`observables_impurity`
     !  #. Get the impurity local energy: call :f:func:`local_energy_impurity`
     !  #. Get the impurity reduced density matric: call :f:func:`rdm_impurity`
     !  #. Delete MPI environment and deallocate used structures :f:var:`state_list` and :f:var:`dmft_bath`
     !
     module procedure :: ed_solve_single
     module procedure :: ed_solve_lattice
  end interface ed_solve


  
  !> FINALIZE SOLVER AND CLEANUP ENVIRONMENT
  interface ed_finalize_solver
     ! 
     ! Finalize the Exact Diagonalization solver, clean up all the allocated memory. 
     !
     module procedure :: ed_finalize_solver_single
     module procedure :: ed_finalize_solver_lattice
  end interface ed_finalize_solver



  public :: ed_init_solver
  public :: ed_solve
  public :: ed_finalize_solver



  logical,save :: isetup=.true. !Allow :f:func:`init_ed_structure` and :f:func:`setup_global` to be called. Set to :f:code:`.false.` by :f:func:`ed_init_solver`, reset by :f:func:`ed_finalize_solver` . Default :code:`.true.`



contains



  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: allocate and initialize one or multiple baths -+!
  subroutine ed_init_solver_single(bath)
    real(8),dimension(:),intent(inout) :: bath !user bath input array
    logical                            :: check
    integer                            :: i
    !
    !SET THE MPI FRAMEWORK:
#ifdef _MPI
    if(check_MPI())call ed_set_MpiComm()
#endif
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
    !DELETE THE MPI FRAMEWORK:
#ifdef _MPI
    if(check_MPI())call ed_del_MpiComm()
#endif
    !
  end subroutine ed_init_solver_single

  subroutine ed_init_solver_lattice(bath)
    real(8),dimension(:,:),intent(inout) :: bath !user bath input array
    integer                              :: ilat,Nineq
    !
    if(allocated(dens_ineq))deallocate(dens_ineq)
    if(allocated(docc_ineq))deallocate(docc_ineq)
    if(allocated(mag_ineq))deallocate(mag_ineq)
    if(allocated(phisc_ineq))deallocate(phisc_ineq)
    if(allocated(e_ineq))deallocate(e_ineq)
    if(allocated(single_particle_density_matrix_ineq))deallocate(single_particle_density_matrix_ineq)
    if(allocated(impurity_density_matrix_ineq))deallocate(impurity_density_matrix_ineq)
    if(allocated(neigen_sector_ineq))deallocate(neigen_sector_ineq)
    if(allocated(neigen_total_ineq))deallocate(neigen_total_ineq)
    !
    Nineq = size(bath,1)
    if(bath_type=='replica'.AND..not.allocated(Hreplica_lambda_ineq))&
         stop "ERROR ed_init_solver: replica parameters lambda not defined for all sites"
    if(bath_type=='general'.AND..not.allocated(Hgeneral_lambda_ineq))&
         stop "ERROR ed_init_solver: general parameters lambda not defined for all sites"
    !
    allocate(dens_ineq(Nineq,Norb))
    allocate(docc_ineq(Nineq,Norb))
    allocate(mag_ineq(Nineq,3,Norb))
    allocate(phisc_ineq(Nineq,Norb,Norb))
    allocate(e_ineq(Nineq,4))
    allocate(dd_ineq(Nineq,4))
    allocate(single_particle_density_matrix_ineq(Nineq,Nspin,Nspin,Norb,Norb))
    allocate(impurity_density_matrix_ineq(Nineq,4**Norb,4**Norb))
    !
    do ilat=1,Nineq
       call ed_set_suffix(ilat)
       if(bath_type=='replica')call Hreplica_site(ilat)
       if(bath_type=='general')call Hgeneral_site(ilat)
       !set the ilat-th lambda vector basis for the replica bath
       call ed_init_solver(bath(ilat,:))
    enddo
#ifdef _MPI
    if(check_MPI())call Barrier_MPI()
#endif
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
  end subroutine ed_init_solver_lattice







  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################








  !+-----------------------------------------------------------------------------+!
  !PURPOSE: solve the impurity problems for a single or many independent
  ! lattice site using ED. 
  !+-----------------------------------------------------------------------------+!
  subroutine ed_solve_single(bath,sflag,fmpi)
    real(8),dimension(:),intent(in)     :: bath  !user bath input array
    logical,optional                    :: sflag !flag to calculate ( :code:`.true.` ) or not ( :code:`.false.` ) Green's functions and susceptibilities. Default :code:`.true.` . 
    logical,optional                    :: fmpi  !flag to solve the impurity problem parallely ( :code:`.true.` ) or not ( :code:`.false.` ). Default :code:`.true.` . 
    logical                             :: fmpi_
    logical                             :: check,iflag
    !
    fmpi_=.true.;if(present(fmpi))fmpi_=fmpi
    iflag=.true. ;if(present(sflag))iflag=sflag
    !
#ifdef _MPI    
    if(check_MPI().AND.fmpi_)call ed_set_MpiComm()
#endif
    !
    if(.not.allocated(impHloc))stop "ED_SOLVE ERROR: impHloc not allocated. Please call ed_set_Hloc first."
    !
    check   = check_bath_dimension(bath)
    if(.not.check)stop "ED_SOLVE_SINGLE Error: wrong bath dimensions"
    !  
    if(MpiMaster.and.fmpi_)call save_input_file(str(ed_input_file))
    !
    call allocate_dmft_bath(dmft_bath)
    call set_dmft_bath(bath,dmft_bath)
    call write_dmft_bath(dmft_bath)
    call save_dmft_bath(dmft_bath,used=.true.)
    !
    !SOLVE THE QUANTUM IMPURITY PROBLEM:
    call diagonalize_impurity()
    if(iflag)then
       call buildgf_impurity()
       call buildchi_impurity()
    endif
    call observables_impurity()
    call local_energy_impurity()
    call rdm_impurity()
    !
    call deallocate_dmft_bath(dmft_bath)
    call es_delete_espace(state_list)
    !
    !DELETE THE LOCAL MPI COMMUNICATOR:
#ifdef _MPI    
    if(check_MPI().AND.fmpi_)call ed_del_MpiComm()
#endif    
    !
    nullify(spHtimesV_cc)
    nullify(spHtimesV_p)
    write(Logfile,"(A)")""
  end subroutine ed_solve_single



  subroutine ed_solve_lattice(bath,mpi_lanc,Uloc_ii,Ust_ii,Jh_ii,Jp_ii,Jx_ii,iflag)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    real(8)                                         :: bath(:,:) !user bath input array
    logical,optional                                :: mpi_lanc  !parallelization strategy flag: if :code:`.false.` each core serially solves an inequivalent site, if :code:`.true.` all cores parallely solve each site in sequence. Default :code:`.false.` .
    real(8),optional,dimension(size(bath,1),Norb)   :: Uloc_ii !site-dependent values for :f:var:`uloc` , overriding the ones in the input file. It has dimension [ :f:var:`nlat` , :f:var:`norb` ].
    real(8),optional,dimension(size(bath,1))        :: Ust_ii  !site-dependent values for :f:var:`ust` , overriding the ones in the input file. It has dimension [ :f:var:`nlat`].
    real(8),optional,dimension(size(bath,1))        :: Jh_ii   !site-dependent values for :f:var:`jh` , overriding the ones in the input file. It has dimension [ :f:var:`nlat`].
    real(8),optional,dimension(size(bath,1))        :: Jp_ii   !site-dependent values for :f:var:`jp` , overriding the ones in the input file. It has dimension [ :f:var:`nlat`].
    real(8),optional,dimension(size(bath,1))        :: Jx_ii   !site-dependent values for :f:var:`jx` , overriding the ones in the input file. It has dimension [ :f:var:`nlat`].
    logical,optional                                :: iflag   !flag to calculate ( :code:`.true.` ) or not ( :code:`.false.` ) Green's functions and susceptibilities. Default :code:`.true.` . 
    !
    !MPI  auxiliary vars
    complex(8)                          :: Dmats_tmp(size(bath,1),0:Lmats)
    complex(8)                          :: Dreal_tmp(size(bath,1),Lreal)
    real(8)                             :: dens_tmp(size(bath,1),Norb)
    real(8)                             :: docc_tmp(size(bath,1),Norb)
    real(8)                             :: mag_tmp(size(bath,1),3,Norb)
    real(8)                             :: phisc_tmp(size(bath,1),Norb,Norb)
    real(8)                             :: e_tmp(size(bath,1),4)
    real(8)                             :: dd_tmp(size(bath,1),4)
    !    
    complex(8)       :: single_particle_density_matrix_tmp(size(bath,1),Nspin,Nspin,Norb,Norb)
    complex(8)       :: impurity_density_matrix_tmp(size(bath,1),4**Norb,4**Norb)
    !
    integer                             :: neigen_sectortmp(size(bath,1),Nsectors)
    integer                             :: neigen_totaltmp(size(bath,1))
    ! 
    integer                             :: i,j,ilat,iorb,jorb,ispin,jspin
    integer                             :: Nineq
    logical                             :: check_dim,mpi_lanc_,iflag_
    character(len=5)                    :: tmp_suffix
    !
    integer                             :: MPI_ID=0
    integer                             :: MPI_SIZE=1
    logical                             :: MPI_MASTER=.true.
    !
    integer                             :: mpi_err 
    !
#ifdef _MPI    
    if(check_MPI())then
       MPI_ID     = get_Rank_MPI()
       MPI_SIZE   = get_Size_MPI()
       MPI_MASTER = get_Master_MPI()
    endif
#endif
    !
    mpi_lanc_=.false.;if(present(mpi_lanc))mpi_lanc_=mpi_lanc
    !
    iflag_=.true.
    if(present(iflag)) iflag_=iflag
    !
    ! Check dimensions
    Nineq=size(bath,1)
    !
    if(size(neigen_sector_ineq,1)<Nineq)stop "ed_solve_lattice error: size(neigen_sectorii,1)<Nineq"
    if(size(neigen_total_ineq)<Nineq)stop "ed_solve_lattice error: size(neigen_totalii,1)<Nineq"
    !
    if(.not.allocated(Hloc_ineq))stop "ed_solve_lattice error: Hloc_ineq not allocated. Call ed_set_Hloc first."
    !Check the dimensions of the bath are ok.
    !This can always be done in parallel no issues with mpi_lanc
    do ilat=1+MPI_ID,Nineq,MPI_SIZE
       check_dim = check_bath_dimension(bath(ilat,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    !
    Dmats_ph_ineq = zero ; Dreal_ph_ineq = zero 
    dens_ineq     = 0d0  ; docc_ineq     = 0d0
    mag_ineq      = 0d0  ; phisc_ineq    = 0d0  
    e_ineq        = 0d0  ; dd_ineq       = 0d0 
    single_particle_density_matrix_ineq = zero
    impurity_density_matrix_ineq = zero
    !
    Dmats_tmp  = zero ; Dreal_tmp  = zero
    dens_tmp   = 0d0  ; docc_tmp   = 0d0
    mag_tmp    = 0d0  ; phisc_tmp  = 0d0
    e_tmp      = 0d0  ; dd_tmp     = 0d0
    neigen_sectortmp = 0
    neigen_totaltmp  = 0
    single_particle_density_matrix_tmp = zero
    impurity_density_matrix_tmp = zero
    !
    !
    select case(mpi_lanc_)
    case default              !mpi_lanc=False => solve sites with MPI
       if(MPI_MASTER)call start_timer(unit=LOGfile)
       do ilat = 1 + MPI_ID, Nineq, MPI_SIZE
          write(LOGfile,*)"CPU: "//str(MPI_ID)//" SOLVES INEQ SITE: "//str(ilat,Npad=4)
          !
          call ed_set_suffix(ilat)
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
          call ed_set_Hloc(Hloc_ineq(ilat,:,:,:,:))
          !
          if(MPI_MASTER)call save_input_file(str(ed_input_file))
          !
          call ed_solve_single(bath(ilat,:),sflag=iflag_,fmpi=mpi_lanc_)
          !
          neigen_sectortmp(ilat,:)   = neigen_sector(:)
          neigen_totaltmp(ilat)      = lanc_nstates_total
          dens_tmp(ilat,1:Norb)      = ed_dens(1:Norb)
          docc_tmp(ilat,1:Norb)      = ed_docc(1:Norb)
          mag_tmp(ilat,:,1:Norb)     = ed_mag(:,1:Norb)
          phisc_tmp(ilat,1:Norb,1:Norb)     = ed_phisc(1:Norb,1:Norb)
          e_tmp(ilat,:)              = [ed_Epot,ed_Eint,ed_Ehartree,ed_Eknot]
          dd_tmp(ilat,:)             = [ed_Dust,ed_Dund,ed_Dse,ed_Dph]
          single_particle_density_matrix_tmp(ilat,:,:,:,:) = single_particle_density_matrix(:,:,:,:)
          impurity_density_matrix_tmp(ilat,:,:) = impurity_density_matrix(:,:)
       enddo
#ifdef _MPI
       call Barrier_MPI()
#endif
       if(MPI_MASTER)call stop_timer
       call ed_reset_suffix
       !
#ifdef _MPI
       if(check_MPI())then
          call AllReduce_MPI(MPI_COMM_WORLD,dens_tmp,dens_ineq)
          call AllReduce_MPI(MPI_COMM_WORLD,docc_tmp,docc_ineq)
          call AllReduce_MPI(MPI_COMM_WORLD,mag_tmp,mag_ineq)
          call AllReduce_MPI(MPI_COMM_WORLD,phisc_tmp,phisc_ineq)
          call AllReduce_MPI(MPI_COMM_WORLD,e_tmp,e_ineq)
          call AllReduce_MPI(MPI_COMM_WORLD,dd_tmp,dd_ineq)
          call AllReduce_MPI(MPI_COMM_WORLD,impurity_density_matrix_tmp,impurity_density_matrix_ineq)
          neigen_sector_ineq=0
          neigen_total_ineq=0
          call AllReduce_MPI(MPI_COMM_WORLD,neigen_sectortmp,neigen_sector_ineq)
          call AllReduce_MPI(MPI_COMM_WORLD,neigen_totaltmp,neigen_total_ineq)
       else
          dens_ineq               = dens_tmp
          docc_ineq               = docc_tmp
          mag_ineq                = mag_tmp
          phisc_ineq              = phisc_tmp
          e_ineq                  = e_tmp
          dd_ineq                 = dd_tmp
          neigen_sector_ineq      = neigen_sectortmp
          neigen_total_ineq       = neigen_totaltmp
          impurity_density_matrix_ineq = impurity_density_matrix_tmp
       endif
#else
       dens_ineq               = dens_tmp
       docc_ineq               = docc_tmp
       mag_ineq                = mag_tmp
       phisc_ineq              = phisc_tmp
       e_ineq                  = e_tmp
       dd_ineq                 = dd_tmp
       neigen_sector_ineq      = neigen_sector_tmp
       neigen_total_ineq       = neigen_total_tmp
       impurity_density_matrix_ineq = impurity_density_matrix_tmp
#endif
       !
       !
       !
    case(.true.)                !solve sites serial, Lanczos with MPI
       if(MPI_MASTER)call start_timer(unit=LOGfile)
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
          call ed_set_Hloc(Hloc_ineq(ilat,:,:,:,:))
          call ed_solve_single(bath(ilat,:),sflag=iflag_,fmpi=mpi_lanc_)
          !
          neigen_sector_ineq(ilat,:)  = neigen_sector(:)
          neigen_total_ineq(ilat)     = lanc_nstates_total
          dens_ineq(ilat,1:Norb)      = ed_dens(1:Norb)
          docc_ineq(ilat,1:Norb)      = ed_docc(1:Norb)
          mag_ineq(ilat,:,1:Norb)     = ed_mag(:,1:Norb)
          phisc_ineq(ilat,1:Norb,1:Norb)     = ed_phisc(1:Norb,1:Norb)
          e_ineq(ilat,:)              = [ed_Epot,ed_Eint,ed_Ehartree,ed_Eknot]
          dd_ineq(ilat,:)             = [ed_Dust,ed_Dund,ed_Dse,ed_Dph]
          single_particle_density_matrix_ineq(ilat,:,:,:,:) = single_particle_density_matrix(:,:,:,:)
          impurity_density_matrix_ineq(ilat,:,:) = impurity_density_matrix(:,:)
       enddo
       if(MPI_MASTER)call stop_timer
       call ed_reset_suffix
    end select
    !
  end subroutine ed_solve_lattice



  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: deallocate and finalize one or multiple baths -+!
  !+-----------------------------------------------------------------------------+!
  !                              SINGLE SITE                                      !
  !+-----------------------------------------------------------------------------+!
  subroutine ed_finalize_solver_single()
    logical                            :: check
    integer                            :: i
    !
    !SET THE MPI FRAMEWORK:
#ifdef _MPI
    if(check_MPI())call ed_set_MpiComm()
#endif
    !
    write(LOGfile,"(A)")"FINALIZE SOLVER "
    !
    !just in case deallocate some objects
    call deallocate_dmft_bath(dmft_bath)
    call es_delete_espace(state_list)
    call deallocate_grids
    nullify(spHtimesV_cc)
    nullify(spHtimesV_p)
    Hreplica_status=.false.
    Hgeneral_status=.false.
    !
    !Delete ED Structure & memory
    call delete_ed_structure()
    !
    !Ready to be setup again
    isetup=.true.
    !
    !DELETE THE MPI FRAMEWORK:
#ifdef _MPI
    if(check_MPI())call ed_del_MpiComm()
#endif
    !
  end subroutine ed_finalize_solver_single

  !+-----------------------------------------------------------------------------+!
  !                              Multiple sites                                   !
  !+-----------------------------------------------------------------------------+!

  subroutine ed_finalize_solver_lattice(Nineq)
    integer                              :: ilat
    integer                              :: Nineq !number of inequivalent impurity sites for real-space DMFT
    !
    if(allocated(dens_ineq))deallocate(dens_ineq)
    if(allocated(docc_ineq))deallocate(docc_ineq)
    if(allocated(mag_ineq))deallocate(mag_ineq)
    if(allocated(phisc_ineq))deallocate(phisc_ineq)
    if(allocated(e_ineq))deallocate(e_ineq)
    if(allocated(dd_ineq))deallocate(dd_ineq)

    if(allocated(single_particle_density_matrix_ineq))deallocate(single_particle_density_matrix_ineq)
    if(allocated(impurity_density_matrix_ineq))deallocate(impurity_density_matrix_ineq)
    if(allocated(neigen_sector_ineq))deallocate(neigen_sector_ineq)
    if(allocated(neigen_total_ineq))deallocate(neigen_total_ineq)
    !
    do ilat=1,Nineq
       call ed_set_suffix(ilat)
       call ed_finalize_solver()
    enddo
    call ed_reset_suffix
    !
  end subroutine ed_finalize_solver_lattice










end module ED_MAIN




