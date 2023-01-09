MODULE ED_BATH_FIT
  USE SF_CONSTANTS
  USE SF_OPTIMIZE, only:fmin_cg,fmin_cgplus,fmin_cgminimize
  USE SF_LINALG,   only:eye,zeye,inv,inv_her,operator(.x.)
  USE SF_IOTOOLS,  only:reg,free_unit,txtfy
  USE SF_ARRAYS,   only:arange
  USE SF_MISC,     only:assert_shape
  !
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL  
  USE ED_AUX_FUNX
  USE ED_BATH
  USE ED_BATH_FUNCTIONS
  USE ED_FIT_COMMON
  USE ED_FIT_NORMAL
  USE ED_FIT_HYBRID
  USE ED_FIT_REPLICA
  USE ED_FIT_GENERAL
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif


  implicit none
  private

  interface ed_chi2_fitgf
     module procedure chi2_fitgf_single_normal
     module procedure chi2_fitgf_single_superc
     module procedure chi2_fitgf_lattice_normal
     module procedure chi2_fitgf_lattice_superc
  end interface ed_chi2_fitgf

  public :: ed_chi2_fitgf


contains


  !+----------------------------------------------------------------------+
  !PURPOSE  : Chi^2 fit of the G0/Delta 
  !
  ! - CHI2_FITGF_GENERIC_NORMAL interface for the normal case 
  !   * CHI2_FITGF_GENERIC_NORMAL_NOSPIN interface to fixed spin input
  !+----------------------------------------------------------------------+
  subroutine chi2_fitgf_single_normal(g,bath,ispin,iorb,fmpi)
    complex(8),dimension(..)                         :: g
    real(8),dimension(:)                             :: bath    
    integer,optional                                 :: ispin,iorb
    logical,optional                                 :: fmpi
    !
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: fg
    integer                                          :: ispin_,Liw
    logical                                          :: fmpi_
    !    
#ifdef _DEBUG
    write(Logfile,"(A)")""
    write(Logfile,"(A)")"DEBUG chi2_fitgf_generic_normal_mpi: Start Chi**2 fit"
#endif
    !
    ispin_=1;if(present(ispin))ispin_=ispin
    fmpi_=.true.;if(present(fmpi))fmpi_=fmpi
    !
#ifdef _MPI    
    if(check_MPI().AND.fmpi_)call ed_set_MpiComm()
#endif
    !
    select case(cg_method)
    case default
       stop "ED Error: cg_method > 1"
    case (0)
       if(ed_verbose>2)write(LOGfile,"(A,I1,A,A)")"Chi^2 fit with CG-nr and CG-weight: ",cg_weight," on: ",cg_scheme
    case (1)
       if(ed_verbose>2)write(LOGfile,"(A,I1,A,A)")"Chi^2 fit with CG-minimize and CG-weight: ",cg_weight," on: ",cg_scheme
    end select
    !
    if(MpiMaster)then
       !
       select rank(g)
       rank default
       stop "chi2_fitgf_generic_normal error: input G has a wrong rank: accepted 3 or 5"
       rank (3)
       call assert_shape(g,[Nspin*Norb,Nspin*Norb,Lmats],"chi2_fitgf_generic_normal","g")
       fg = so2nn_reshape(g,Nspin,Norb,Lmats)
       rank (5)
       call assert_shape(g,[Nspin,Nspin,Norb,Norb,Lmats],"chi2_fitgf_generic_normal","g")
       fg = g
       end select
       !
       select case(bath_type)
       case default
          select case(ed_mode)
          case ("normal")
             if(present(iorb))then
                call chi2_fitgf_normal_normal(fg(ispin_,ispin_,:,:,:),bath,ispin_,iorb)
             else
                call chi2_fitgf_normal_normal(fg(ispin_,ispin_,:,:,:),bath,ispin_)
             endif
          case ("nonsu2")
             call chi2_fitgf_normal_nonsu2(fg(:,:,:,:,:),bath)
          case default
             stop "chi2_fitgf ERROR: ed_mode!=normal/nonsu2 but only NORMAL component is provided"
          end select
          !
       case ("hybrid")
          select case(ed_mode)
          case ("normal")
             call chi2_fitgf_hybrid_normal(fg(ispin_,ispin_,:,:,:),bath,ispin_)
          case ("nonsu2")
             call chi2_fitgf_hybrid_nonsu2(fg(:,:,:,:,:),bath)
          case default
             stop "chi2_fitgf ERROR: ed_mode!=normal/nonsu2 but only NORMAL component is provided" 
          end select
          !
       case ("replica")
          select case(ed_mode)
          case ("normal","nonsu2")
             call chi2_fitgf_replica(fg,bath)
          case default
             stop "chi2_fitgf ERROR: ed_mode!=normal/nonsu2 but only NORMAL component is provided" 
          end select
       case ("general")
          select case(ed_mode)
          case ("normal","nonsu2")
             call chi2_fitgf_general(fg,bath)
          case default
             stop "chi2_fitgf ERROR: ed_mode!=normal/nonsu2 but only NORMAL component is provided" 
          end select
       end select
    endif
    !
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,bath)
       if(.not.MpiMaster)write(LOGfile,"(A)")"Bath received from master node"
    endif
#endif
    !
    !set trim_state_list to true after the first fit has been done: this 
    !marks the ends of the cycle of the 1st DMFT loop.
    trim_state_list=.true.
    !DELETE THE LOCAL MPI COMMUNICATOR:
#ifdef _MPI    
    if(check_MPI().AND.fmpi_)call ed_del_MpiComm()
#endif   
#ifdef _DEBUG
    write(Logfile,"(A)")""
#endif
  end subroutine chi2_fitgf_single_normal

  
  subroutine chi2_fitgf_single_superc(g,f,bath,ispin,iorb,fmpi)
    complex(8),dimension(..)                            :: g
    complex(8),dimension(..)                            :: f
    real(8),dimension(:)                                :: bath
    integer,optional                                    :: ispin,iorb
    logical,optional                                    :: fmpi
    !
    complex(8),dimension(2,Nspin,Nspin,Norb,Norb,Lmats) :: fg
    integer                                             :: ispin_
    logical                                             :: fmpi_
#ifdef _DEBUG
    write(Logfile,"(A)")""
    write(Logfile,"(A)")"DEBUG chi2_fitgf_generic_superc: Start Chi**2 fit"
#endif
    !
    !
    ispin_=1;if(present(ispin))ispin_=ispin
    fmpi_=.true.;if(present(fmpi))fmpi_=fmpi
    !
#ifdef _MPI    
    if(check_MPI().AND.fmpi_)call ed_set_MpiComm()
#endif
    !
    select case(cg_method)
    case default
       stop "ED Error: cg_method > 1"
    case (0)
       if(ed_verbose>2)write(LOGfile,"(A,I1,A,A)")"master: Chi^2 fit with CG-nr and CG-weight: ",cg_weight," on: ",cg_scheme
    case (1)
       if(ed_verbose>2)write(LOGfile,"(A,I1,A,A)")"master: Chi^2 fit with CG-minimize and CG-weight: ",cg_weight," on: ",cg_scheme
    end select
    !
    if(MpiMaster)then
       !
       select rank(g)
       rank default
       stop "chi2_fitgf_generic_superc error: input G has a wrong rank: accepted 3 or 5"
       rank (3)
       call assert_shape(g,[Nspin*Norb,Nspin*Norb,Lmats],"chi2_fitgf_generic_superc","g")
       fg(1,:,:,:,:,:) = so2nn_reshape(g,Nspin,Norb,Lmats)
       rank (5)
       call assert_shape(g,[Nspin,Nspin,Norb,Norb,Lmats],"chi2_fitgf_generic_superc","g")
       fg(1,:,:,:,:,:) = g
       end select
       !
       select rank(f)
       rank default
       stop "chi2_fitgf_generic_superc error: input F has a wrong rank: accepted 3 or 5"
       rank (3)
       call assert_shape(f,[Nspin*Norb,Nspin*Norb,Lmats],"chi2_fitgf_generic_superc","f")
       fg(2,:,:,:,:,:) = so2nn_reshape(f,Nspin,Norb,Lmats)
       rank (5)
       call assert_shape(f,[Nspin,Nspin,Norb,Norb,Lmats],"chi2_fitgf_generic_superc","f")
       fg(2,:,:,:,:,:) = f
       end select
       !       
       select case(bath_type)
       case default
          select case(ed_mode)
          case ("superc")
             if(present(iorb))then
                call chi2_fitgf_normal_superc(fg(:,ispin_,ispin_,:,:,:),bath,ispin_,iorb)
             else
                call chi2_fitgf_normal_superc(fg(:,ispin_,ispin_,:,:,:),bath,ispin_)
             endif
          case default
             write(LOGfile,"(A)") "chi2_fitgf WARNING: ed_mode=normal/nonsu2 but NORMAL & ANOMAL components provided."
             call chi2_fitgf_normal_normal(fg(1,ispin_,ispin_,:,:,:),bath,ispin_)          
          end select
       case ("hybrid")
          select case(ed_mode)
          case ("superc")
             call chi2_fitgf_hybrid_superc(fg(:,ispin_,ispin_,:,:,:),bath,ispin_)
          case default
             write(LOGfile,"(A)") "chi2_fitgf WARNING: ed_mode=normal/nonsu2 but NORMAL & ANOMAL components provided."
             call chi2_fitgf_hybrid_normal(fg(1,ispin_,ispin_,:,:,:),bath,ispin_)       
          end select
       case ("replica")
          select case(ed_mode)
          case ("superc")
             call chi2_fitgf_replica_superc(fg(:,:,:,:,:,:),bath)
          case default
             write(LOGfile,"(A)") "chi2_fitgf WARNING: ed_mode=normal/nonsu2 but NORMAL & ANOMAL components provided."
             call chi2_fitgf_replica(fg(1,:,:,:,:,:),bath)
          end select
       case ("general")
          select case(ed_mode)
          case ("superc")
             call chi2_fitgf_general_superc(fg(:,:,:,:,:,:),bath)
          case default
             write(LOGfile,"(A)") "chi2_fitgf WARNING: ed_mode=normal/nonsu2 but NORMAL & ANOMAL components provided."
             call chi2_fitgf_general(fg(1,:,:,:,:,:),bath)
          end select
       end select
    endif
    !
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,bath)
       if(.not.MpiMaster)write(LOGfile,"(A)")"Bath received from master node"
    endif
#endif
    !
    !set trim_state_list to true after the first fit has been done: this 
    !marks the ends of the cycle of the 1st DMFT loop.
    trim_state_list=.true.
#ifdef _MPI    
    if(check_MPI().AND.fmpi_)call ed_del_MpiComm()
#endif   
#ifdef _DEBUG
    write(Logfile,"(A)")""
#endif
  end subroutine chi2_fitgf_single_superc









  subroutine chi2_fitgf_lattice_normal(g,bath,ispin)
    real(8),intent(inout)    :: bath(:,:)
    complex(8),dimension(..) :: g
    integer,optional         :: ispin
    !
    complex(8)               :: fg(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    real(8)                  :: bath_tmp(size(bath,1),size(bath,2))
    integer                  :: ilat,i,iorb,ispin_
    integer                  :: Nsites
    logical                  :: check_dim
    character(len=5)         :: tmp_suffix
    integer                  :: MPI_ID=0
    integer                  :: MPI_SIZE=1
    logical                  :: MPI_MASTER=.true.
    !
    !
#ifdef _DEBUG
    write(Logfile,"(A)")""
#endif
    !
    ispin_=1;if(present(ispin))ispin_=ispin
    if(ispin_>Nspin)stop "ed_fit_bath_sites error: required spin index > Nspin"
    !
#ifdef _MPI    
    if(check_MPI())then
       MPI_ID     = get_Rank_MPI()
       MPI_SIZE   = get_Size_MPI()
       MPI_MASTER = get_Master_MPI()
    endif
#endif
    !
    ! Check dimensions !
    Nsites=size(bath,1)
    !
    do ilat=1+MPI_ID,Nsites,MPI_SIZE
       check_dim = check_bath_dimension(bath(ilat,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    !
    select rank(g)
    rank default
    stop "chi2_fitgf_generic_normal error: input G has a wrong rank: accepted 3 or 5"
    !
    rank (3)
    call assert_shape(g,[Nsites*Nspin*Norb,Nsites*Nspin*Norb,Lmats],'chi2_fitgf_generic_normal','g')
    fg = lso2nnn_reshape(g,Nsites,Nspin,Norb,Lmats)
    !
    rank (4)
    call assert_shape(g,[Nsites,Nspin*Norb,Nspin*Norb,Lmats],'chi2_fitgf_generic_normal','g')
    do ilat=1,Nsites
       fg(ilat,:,:,:,:,:)  = so2nn_reshape(g(ilat,:,:,:),Nspin,Norb,Lmats)
    enddo
    !
    rank (6)
    call assert_shape(g,[Nsites,Nspin,Nspin,Norb,Norb,Lmats],'chi2_fitgf_generic_normal','g')
    fg = g
    end select
    !
    bath_tmp=0d0
    !
    do ilat = 1+MPI_ID,Nsites,MPI_SIZE
#ifdef _DEBUG
       write(Logfile,"(A)")"DEBUG ed_fit_bath_sites_normal: Start Chi**2 fit for site:"//str(ilat)
#endif
       bath_tmp(ilat,:)=bath(ilat,:)
       call ed_set_Hloc(Hloc_ineq(ilat,:,:,:,:))
       ed_file_suffix=reg(ineq_site_suffix)//str(ilat,site_indx_padding)
       !
       if(present(ispin))then
          call ed_chi2_fitgf(fg(ilat,:,:,:,:,:),bath_tmp(ilat,:),ispin=ispin_,fmpi=.false.)
       else
          call ed_chi2_fitgf(fg(ilat,:,:,:,:,:),bath_tmp(ilat,:),fmpi=.false.)
       end if
    end do
    !
#ifdef _MPI
    if(check_MPI())then
       bath=0d0
       call AllReduce_MPI(MPI_COMM_WORLD,bath_tmp,bath)
    else
       bath = bath_tmp
    endif
#else
    bath = bath_tmp
#endif
    !
    call ed_reset_suffix
  end subroutine chi2_fitgf_lattice_normal

  subroutine chi2_fitgf_lattice_superc(g,f,bath,ispin)
    real(8),intent(inout)    :: bath(:,:)
    complex(8),dimension(..) :: g,f
    integer,optional         :: ispin
    !
    complex(8)               :: fg(2,size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    real(8)                  :: bath_tmp(size(bath,1),size(bath,2))
    integer                  :: ilat,i,iorb,ispin_
    integer                  :: Nsites
    logical                  :: check_dim
    character(len=5)         :: tmp_suffix
    integer                  :: MPI_ID=0
    integer                  :: MPI_SIZE=1
    logical                  :: MPI_MASTER=.true.
    !
#ifdef _DEBUG
    write(Logfile,"(A)")""
#endif
    ispin_=1;if(present(ispin))ispin_=ispin
    if(ispin_>Nspin)stop "ed_fit_bath_sites error: required spin index > Nspin"
    !
#ifdef _MPI    
    if(check_MPI())then
       MPI_ID     = get_Rank_MPI()
       MPI_SIZE   = get_Size_MPI()
       MPI_MASTER = get_Master_MPI()
    endif
#endif
    !
    ! Check dimensions !
    Nsites=size(bath,1)
    !
    do ilat = 1 + MPI_ID, Nsites, MPI_SIZE
       check_dim = check_bath_dimension(bath(ilat,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    !
    select rank(g)
    rank default
    stop "chi2_fitgf_generic_superc error: input G has a wrong rank: accepted 3 or 5"
    rank (3)
    call assert_shape(g,[Nsites*Nspin*Norb,Nsites*Nspin*Norb,Lmats],'chi2_fitgf_generic_superc','g')
    fg(1,:,:,:,:,:,:) = lso2nnn_reshape(g,Nsites,Nspin,Norb,Lmats)
    rank (4)
    call assert_shape(g,[Nsites,Nspin*Norb,Nspin*Norb,Lmats],'chi2_fitgf_generic_superc','g')
    do ilat=1,Nsites
       fg(1,ilat,:,:,:,:,:)  = so2nn_reshape(g(ilat,:,:,:),Nspin,Norb,Lmats)
    enddo
    rank (6)
    call assert_shape(g,[Nsites,Nspin,Nspin,Norb,Norb,Lmats],'chi2_fitgf_generic_superc','g')
    fg(1,:,:,:,:,:,:) = g
    end select
    !
    select rank(f)
    rank default
    stop "chi2_fitgf_generic_superc error: input F has a wrong rank: accepted 3 or 5"
    rank (3)
    call assert_shape(f,[Nsites*Nspin*Norb,Nsites*Nspin*Norb,Lmats],'chi2_fitgf_generic_superc','f')
    fg(2,:,:,:,:,:,:) = lso2nnn_reshape(f,Nsites,Nspin,Norb,Lmats)
    rank (4)
    call assert_shape(f,[Nsites,Nspin*Norb,Nspin*Norb,Lmats],'chi2_fitgf_generic_superc','f')
    do ilat=1,Nsites
       fg(2,ilat,:,:,:,:,:)  = so2nn_reshape(f(ilat,:,:,:),Nspin,Norb,Lmats)
    enddo
    rank (6)
    call assert_shape(f,[Nsites,Nspin,Nspin,Norb,Norb,Lmats],'chi2_fitgf_generic_superc','f')
    fg(2,:,:,:,:,:,:) = f
    end select
    !
    bath_tmp=0.d0
    !
    do ilat= 1 + MPI_ID, Nsites, MPI_SIZE
#ifdef _DEBUG
       write(Logfile,"(A)")"DEBUG ed_fit_bath_sites_superc: Start Chi**2 fit for site:"//str(ilat)
#endif
       bath_tmp(ilat,:) = bath(ilat,:)
       call ed_set_Hloc(Hloc_ineq(ilat,:,:,:,:))
       ed_file_suffix=reg(ineq_site_suffix)//str(ilat,site_indx_padding)
       !
       if(present(ispin))then
          call ed_chi2_fitgf(fg(:,ilat,:,:,:,:,:),bath_tmp(ilat,:),ispin_,fmpi=.false.)
       else
          do ispin_=1,Nspin
             call ed_chi2_fitgf(fg(:,ilat,:,:,:,:,:),bath_tmp(ilat,:),ispin_,fmpi=.false.)
          enddo
       endif
    end do
    !
#ifdef _MPI
    if(check_MPI())then
       bath=0d0
       call AllReduce_MPI(MPI_COMM_WORLD,bath_tmp,bath)
    else
       bath = bath_tmp
    endif
#else
    bath = bath_tmp
#endif
    !
    call ed_reset_suffix
  end subroutine chi2_fitgf_lattice_superc




  
  
end MODULE ED_BATH_FIT
















