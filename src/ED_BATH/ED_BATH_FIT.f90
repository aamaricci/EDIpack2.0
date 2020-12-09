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
  USE ED_HLOC_DECOMPOSITION
  USE ED_BATH
  USE ED_BATH_FUNCTIONS
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif


  implicit none
  private

  interface ed_chi2_fitgf
     !SINGLE_SITE
     module procedure chi2_fitgf_generic_normal
     module procedure chi2_fitgf_generic_superc
#ifdef _MPI
     module procedure chi2_fitgf_generic_normal_mpi
     module procedure chi2_fitgf_generic_superc_mpi
#endif
     !RDMFT_WRAPPER
     module procedure ed_fit_bath_sites_normal
     module procedure ed_fit_bath_sites_superc
#ifdef _MPI
     module procedure ed_fit_bath_sites_normal_mpi
     module procedure ed_fit_bath_sites_superc_mpi
#endif
  end interface ed_chi2_fitgf

  public :: ed_chi2_fitgf


  integer                               :: Ldelta
  complex(8),dimension(:,:),allocatable :: Gdelta
  complex(8),dimension(:,:),allocatable :: Fdelta
  real(8),dimension(:),allocatable      :: Xdelta,Wdelta
  integer                               :: totNorb,totNspin,totNso
  integer,dimension(:),allocatable      :: getIorb,getJorb,getIspin,getJspin
  integer                               :: Orb_indx,Spin_indx,Spin_mask
  type(effective_bath)                  :: chi2_bath
  integer                               :: cg_iter_count=0
  logical                               :: para_

  integer                               :: MPI_RANK=0
  integer                               :: MPI_SIZE=1
  logical                               :: MPI_MASTER=.true.
  integer                               :: MPI_IERR

  !This contains the number of the lambda expansion
  !for each replica of the impurity
  integer                              :: Nlambdas
  !
  !This is a dummy object which is used here to point
  !to the replica bath lambdas, i.e. the coefficients
  !of the bath item-th Hamiltonian expansion 
  type nsymm_vector
     real(8),dimension(:),allocatable   :: element          
  end type nsymm_vector

contains


  !+----------------------------------------------------------------------+
  !PURPOSE  : Chi^2 fit of the G0/Delta 
  !
  ! - CHI2_FITGF_GENERIC_NORMAL interface for the normal case 
  !   * CHI2_FITGF_GENERIC_NORMAL_NOSPIN interface to fixed spin input
  !+----------------------------------------------------------------------+
  subroutine chi2_fitgf_generic_normal(fg,bath,ispin,iorb)
    complex(8),dimension(:,:,:,:,:) :: fg ![Nspin][Nspin][Norb][Norb][Niw] 
    real(8),dimension(:)            :: bath
    integer,optional                :: ispin,iorb
    integer                         :: ispin_
    ispin_=1;if(present(ispin))ispin_=ispin
    call assert_shape(fg,[Nspin,Nspin,Norb,Norb,size(fg,5)],"chi2_fitgf_generic_normal","fg")
    !
    select case(cg_method)
    case default
       stop "ED Error: cg_method > 1"
    case (0)
       if(ed_verbose>2)write(LOGfile,"(A,I1,A,A)")"\Chi2 fit with CG-nr and CG-weight: ",cg_weight," on: ",cg_scheme
    case (1)
       if(ed_verbose>2)write(LOGfile,"(A,I1,A,A)")"\Chi2 fit with CG-minimize and CG-weight: ",cg_weight," on: ",cg_scheme
    end select
    !
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
       call chi2_fitgf_replica(fg,bath)
       !
    end select
    !
    !set trim_state_list to true after the first fit has been done: this 
    !marks the ends of the cycle of the 1st DMFT loop.
    trim_state_list=.true.
  end subroutine chi2_fitgf_generic_normal

  subroutine chi2_fitgf_generic_superc(fg,bath,ispin,iorb)
    complex(8),dimension(:,:,:,:,:,:) :: fg ![2][Nspin][Nspin][Norb][Norb][Niw]
    real(8),dimension(:)            :: bath
    integer,optional                :: ispin,iorb
    integer                         :: ispin_
    ispin_=1;if(present(ispin))ispin_=ispin
    call assert_shape(fg,[2,Nspin,Nspin,Norb,Norb,size(fg,6)],"chi2_fitgf_generic_superc","fg")
    !
    select case(cg_method)
    case default
       stop "ED Error: cg_method > 1"
    case (0)
       if(ed_verbose>2)write(LOGfile,"(A,I1,A,A)")"\Chi2 fit with CG-nr and CG-weight: ",cg_weight," on: ",cg_scheme
    case (1)
       if(ed_verbose>2)write(LOGfile,"(A,I1,A,A)")"\Chi2 fit with CG-minimize and CG-weight: ",cg_weight," on: ",cg_scheme
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
          call sleep(1)
          call chi2_fitgf_normal_normal(fg(1,ispin_,ispin_,:,:,:),bath,ispin_)          
       end select
    case ("hybrid")
       select case(ed_mode)
       case ("superc")
          call chi2_fitgf_hybrid_superc(fg(:,ispin_,ispin_,:,:,:),bath,ispin_)
       case default
          write(LOGfile,"(A)") "chi2_fitgf WARNING: ed_mode=normal/nonsu2 but NORMAL & ANOMAL components provided."
          call sleep(1)
          call chi2_fitgf_hybrid_normal(fg(1,ispin_,ispin_,:,:,:),bath,ispin_)          
       end select
    end select
    !set trim_state_list to true after the first fit has been done: this 
    !marks the ends of the cycle of the 1st DMFT loop.
    trim_state_list=.true.
  end subroutine chi2_fitgf_generic_superc







#ifdef _MPI
  subroutine chi2_fitgf_generic_normal_mpi(comm,fg,bath,ispin,iorb)
    integer                         :: comm
    complex(8),dimension(:,:,:,:,:) :: fg ![Nspin][Nspin][Norb][Norb][Niw] 
    real(8),dimension(:)            :: bath
    integer,optional                :: ispin,iorb
    integer                         :: ispin_
    !
    MPI_MASTER=get_Master_MPI(comm)
    !
    ispin_=1;if(present(ispin))ispin_=ispin
    !
    call assert_shape(fg,[Nspin,Nspin,Norb,Norb,size(fg,5)],"chi2_fitgf_generic_normal","fg")
    select case(cg_method)
    case default
       stop "ED Error: cg_method > 1"
    case (0)
       if(ed_verbose>2)write(LOGfile,"(A,I1,A,A)")"master: Chi^2 fit with CG-nr and CG-weight: ",cg_weight," on: ",cg_scheme
    case (1)
       if(ed_verbose>2)write(LOGfile,"(A,I1,A,A)")"master: Chi^2 fit with CG-minimize and CG-weight: ",cg_weight," on: ",cg_scheme
    end select
    !
    if(MPI_MASTER)then
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
          call chi2_fitgf_replica(fg,bath)
       end select
    endif
    !
    call Bcast_MPI(comm,bath)
    if(.not.MPI_MASTER)write(LOGfile,"(A)")"Bath received from master node"
    !
    !set trim_state_list to true after the first fit has been done: this 
    !marks the ends of the cycle of the 1st DMFT loop.
    trim_state_list=.true.
  end subroutine chi2_fitgf_generic_normal_mpi

  subroutine chi2_fitgf_generic_superc_mpi(comm,fg,bath,ispin,iorb)
    integer                           :: comm
    complex(8),dimension(:,:,:,:,:,:) :: fg ![2][Nspin][Nspin][Norb][Norb][Niw]
    real(8),dimension(:)              :: bath
    integer,optional                  :: ispin,iorb
    integer                           :: ispin_
    !
    MPI_MASTER=get_Master_MPI(comm)
    !
    ispin_=1;if(present(ispin))ispin_=ispin
    call assert_shape(fg,[2,Nspin,Nspin,Norb,Norb,size(fg,6)],"chi2_fitgf_generic_superc","fg")
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
    if(MPI_MASTER)then
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
             call sleep(1)
             call chi2_fitgf_normal_normal(fg(1,ispin_,ispin_,:,:,:),bath,ispin_)          
          end select
       case ("hybrid")
          select case(ed_mode)
          case ("superc")
             call chi2_fitgf_hybrid_superc(fg(:,ispin_,ispin_,:,:,:),bath,ispin_)
          case default
             write(LOGfile,"(A)") "chi2_fitgf WARNING: ed_mode=normal/nonsu2 but NORMAL & ANOMAL components provided."
             call sleep(1)
             call chi2_fitgf_hybrid_normal(fg(1,ispin_,ispin_,:,:,:),bath,ispin_)          
          end select
       end select
    endif
    !
    call Bcast_MPI(comm,bath)
    if(.not.MPI_MASTER)write(LOGfile,"(A)")"Bath received from master node"
    !
    !set trim_state_list to true after the first fit has been done: this 
    !marks the ends of the cycle of the 1st DMFT loop.
    trim_state_list=.true.
  end subroutine chi2_fitgf_generic_superc_mpi
#endif




  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************
  !NORMAL bath
  include "fitgf_normal_normal.f90"
  include "fitgf_normal_superc.f90"
  include "fitgf_normal_nonsu2.f90"

  !HYBRID bath
  include "fitgf_hybrid_normal.f90"  
  include "fitgf_hybrid_superc.f90"  
  include "fitgf_hybrid_nonsu2.f90"

  !REPLICA bath
  include "fitgf_replica.f90"
  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************










  !+----------------------------------------------------------------------!
  ! PURPOSE: given a number of independent baths, evaluate N independent
  ! Delta/G0 functions and fit them to update the effective baths for ED.
  !+----------------------------------------------------------------------!
  !RDMFT WRAPPER:
  subroutine ed_fit_bath_sites_normal(bath,Delta,Hloc,ispin)
    real(8),intent(inout)    :: bath(:,:)
    complex(8),intent(inout) :: Delta(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)               :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    integer,optional         :: ispin
    !MPI auxiliary vars
    real(8)                  :: bath_tmp(size(bath,1),size(bath,2))
    integer                  :: ilat,i,iorb,ispin_
    integer                  :: Nsites
    logical                  :: check_dim
    character(len=5)         :: tmp_suffix
    !
    ! Check dimensions !
    Nsites=size(bath,1)
    !
    do ilat=1,Nsites
       check_dim = check_bath_dimension(bath(ilat,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    !
    bath_tmp=0d0
    do ilat = 1, Nsites
       bath_tmp(ilat,:)=bath(ilat,:)
       call set_Hloc(Hloc(ilat,:,:,:,:))
       !
       ed_file_suffix=reg(ineq_site_suffix)//str(ilat,site_indx_padding)
       !
       if(present(ispin))then
          ispin_=ispin
          if(ispin_>Nspin)stop "ed_fit_bath_sites error: required spin index > Nspin"
          call ed_chi2_fitgf(Delta(ilat,:,:,:,:,:),bath_tmp(ilat,:),ispin=ispin_)
       else
          call ed_chi2_fitgf(Delta(ilat,:,:,:,:,:),bath_tmp(ilat,:))
       end if
    end do
    !
    bath = bath_tmp
    !
    ed_file_suffix=""
  end subroutine ed_fit_bath_sites_normal

  subroutine ed_fit_bath_sites_superc(bath,Delta,Hloc,ispin)
    integer                  :: comm
    real(8),intent(inout)    :: bath(:,:)
    complex(8),intent(inout) :: Delta(2,size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)               :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    integer,optional         :: ispin
    !MPI auxiliary
    real(8)                  :: bath_tmp(size(bath,1),size(bath,2))
    integer                  :: ilat,i,iorb,ispin_,Nsites
    logical                  :: check_dim
    character(len=5)         :: tmp_suffix
    !
    Nsites=size(bath,1)
    !
    do ilat = 1,Nsites
       check_dim = check_bath_dimension(bath(ilat,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    !
    bath_tmp=0.d0
    !
    do ilat=1, Nsites
       !
       bath_tmp(ilat,:) = bath(ilat,:)
       !
       call set_Hloc(Hloc(ilat,:,:,:,:))
       !
       ed_file_suffix=reg(ineq_site_suffix)//str(ilat,site_indx_padding)
       !
       if(present(ispin))then
          ispin_=ispin
          if(ispin_>Nspin)stop "ed_fit_bath_sites error: required spin index > Nspin"
          call ed_chi2_fitgf(Delta(:,ilat,:,:,:,:,:),bath_tmp(ilat,:),ispin_)
       else
          do ispin_=1,Nspin
             call ed_chi2_fitgf(Delta(:,ilat,:,:,:,:,:),bath_tmp(ilat,:),ispin_)
          enddo
       endif
    end do
    bath = bath_tmp
    ed_file_suffix=""
  end subroutine ed_fit_bath_sites_superc




#ifdef _MPI
  subroutine ed_fit_bath_sites_normal_mpi(comm,bath,Delta,Hloc,ispin)
    integer                  :: comm
    real(8),intent(inout)    :: bath(:,:)
    complex(8),intent(inout) :: Delta(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)               :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    integer,optional         :: ispin
    !MPI auxiliary vars
    real(8)                  :: bath_tmp(size(bath,1),size(bath,2))
    integer                  :: ilat,i,iorb,ispin_
    integer                  :: Nsites
    logical                  :: check_dim
    character(len=5)         :: tmp_suffix
    !
    MPI_RANK = get_Rank_MPI(comm)
    MPI_SIZE = get_Size_MPI(comm)
    MPI_MASTER=get_Master_MPI(comm)
    !
    ! Check dimensions !
    Nsites=size(bath,1)
    !
    do ilat=1+MPI_RANK,Nsites,MPI_SIZE
       check_dim = check_bath_dimension(bath(ilat,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    !
    bath_tmp=0d0
    do ilat = 1+MPI_RANK,Nsites,MPI_SIZE
       bath_tmp(ilat,:)=bath(ilat,:)
       call set_Hloc(Hloc(ilat,:,:,:,:))
       !
       ed_file_suffix=reg(ineq_site_suffix)//str(ilat,site_indx_padding)
       !
       if(present(ispin))then
          ispin_=ispin
          if(ispin_>Nspin)stop "ed_fit_bath_sites error: required spin index > Nspin"
          call ed_chi2_fitgf(Delta(ilat,:,:,:,:,:),bath_tmp(ilat,:),ispin=ispin_)
       else
          call ed_chi2_fitgf(Delta(ilat,:,:,:,:,:),bath_tmp(ilat,:))
       end if
    end do
    !
    bath=0d0
    call MPI_ALLREDUCE(bath_tmp(:,:),bath(:,:),size(bath),MPI_DOUBLE_PRECISION,MPI_SUM,comm,MPI_IERR)
    !
    ed_file_suffix=""
  end subroutine ed_fit_bath_sites_normal_mpi

  subroutine ed_fit_bath_sites_superc_mpi(comm,bath,Delta,Hloc,ispin)
    integer                  :: comm
    real(8),intent(inout)    :: bath(:,:)
    complex(8),intent(inout) :: Delta(2,size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)               :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    integer,optional         :: ispin
    !MPI auxiliary
    real(8)                  :: bath_tmp(size(bath,1),size(bath,2))
    integer                  :: ilat,i,iorb,ispin_,Nsites
    logical                  :: check_dim
    character(len=5)         :: tmp_suffix
    !
    MPI_RANK = get_Rank_MPI(comm)
    MPI_SIZE = get_Size_MPI(comm)
    MPI_MASTER=get_Master_MPI(comm)
    !
    Nsites=size(bath,1)
    !
    do ilat = 1 + MPI_RANK, Nsites, MPI_SIZE
       check_dim = check_bath_dimension(bath(ilat,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    !
    bath_tmp=0.d0
    !
    do ilat= 1 + MPI_RANK, Nsites, MPI_SIZE
       !
       bath_tmp(ilat,:) = bath(ilat,:)
       !
       call set_Hloc(Hloc(ilat,:,:,:,:))
       !
       ed_file_suffix=reg(ineq_site_suffix)//str(ilat,site_indx_padding)
       !
       if(present(ispin))then
          ispin_=ispin
          if(ispin_>Nspin)stop "ed_fit_bath_sites error: required spin index > Nspin"
          call ed_chi2_fitgf(Delta(:,ilat,:,:,:,:,:),bath_tmp(ilat,:),ispin_)
       else
          do ispin_=1,Nspin
             call ed_chi2_fitgf(Delta(:,ilat,:,:,:,:,:),bath_tmp(ilat,:),ispin_)
          enddo
       endif
    end do
    !
    bath=0d0
    call MPI_ALLREDUCE(bath_tmp,bath,size(bath),MPI_DOUBLE_PRECISION,MPI_SUM,comm,MPI_IERR)
    !
    ed_file_suffix=""
  end subroutine ed_fit_bath_sites_superc_mpi
#endif

end MODULE ED_BATH_FIT
















