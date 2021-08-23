MODULE ED_HAMILTONIAN_SUPERC
  USE ED_HAMILTONIAN_SUPERC_COMMON
  USE ED_HAMILTONIAN_SUPERC_STORED_HxV
  USE ED_HAMILTONIAN_SUPERC_DIRECT_HxV
  !
  implicit none
  private

  !>Build sparse hamiltonian of the sector
  public  :: build_Hv_sector_superc
  public  :: delete_Hv_sector_superc
  public  :: vecDim_Hv_sector_superc


  !> Tridiag sparse Hamiltonian of the sector
  public  :: tridiag_Hv_sector_superc





contains






  !####################################################################
  !                 MAIN ROUTINES: BUILD/DELETE SECTOR
  !####################################################################
  subroutine build_Hv_sector_superc(isector,Hmat)
    integer                            :: isector,SectorDim
    complex(8),dimension(:,:),optional :: Hmat
    integer                            :: irank
    integer                            :: i,j,Dim
    !
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")&
         "DEBUG build_Hv_sector_SUPERC: build H*v info. present(Hmat):"//str(present(Hmat))//&
         ", total Nup,Ndw:"//str(ed_total_ud)//", using sparse H:"//str(ed_sparse_H)
#endif
    !
    call build_sector(isector,Hsector)
    !
    Dim    = Hsector%Dim
    !
    !#################################
    !          MPI SETUP
    !#################################
    mpiAllThreads=.true.
    MpiQ = Dim/MpiSize
    MpiR = 0
    if(MpiRank==(MpiSize-1))MpiR=mod(Dim,MpiSize)
    !
    MpiIshift = MpiRank*mpiQ
    MpiIstart = MpiRank*mpiQ + 1
    MpiIend   = (MpiRank+1)*mpiQ + mpiR
    !
#ifdef _MPI
#ifdef _DEBU
    if(MpiStatus.AND.ed_verbose>4)then
       write(LOGfile,*)&
            "         mpiRank,   mpi_Q,   mpi_R,   mpi_Istart,   mpi_Iend,   mpi_Iend-mpi_Istart"
       do irank=0,MpiSize-1
          call Barrier_MPI(MpiComm)
          if(MpiRank==irank)then
             write(LOGfile,*)MpiRank,MpiQ,MpiR,MpiIstart,MpiIend,MpiIend-MpiIstart+1
          endif
       enddo
       call Barrier_MPI(MpiComm)
    endif
#endif
#endif
    !
    !
    !#################################
    !          HxV SETUP
    !#################################
    if(present(Hmat))then
       spHtimesV_cc => null()
       call ed_buildh_superc_main(isector,Hmat)
       return
    endif
    !
    select case (ed_sparse_H)
    case (.true.)
       spHtimesV_cc => spMatVec_superc_main
#ifdef _MPI
       if(MpiStatus)spHtimesV_cc => spMatVec_MPI_superc_main
#endif
       call ed_buildh_superc_main(isector)
       !
       !
    case (.false.)
#ifdef _DEBUG
       if(ed_verbose>2)write(Logfile,"(A)")"DEBUG ed_build_Hv_sector SUPERC: direct H*v product, no further debug info..."
#endif
       spHtimesV_cc => directMatVec_superc_main
#ifdef _MPI
       if(MpiStatus)spHtimesV_cc => directMatVec_MPI_superc_main
#endif
    end select
    !
  end subroutine build_Hv_sector_superc





  
  subroutine delete_Hv_sector_superc()
    integer :: iud
    !
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")"DEBUG delete_Hv_sector_SUPERC: delete H*v info"
#endif
    !
    call delete_sector(Hsector)
    Dim    = 0
#ifdef _MPI
    if(MpiStatus)then
       call sp_delete_matrix(MpiComm,spH0)
    else
       call sp_delete_matrix(spH0)
    endif
#else
    call sp_delete_matrix(spH0)
#endif
    !
    spHtimesV_cc => null()
    !
#ifdef _MPI
    if(MpiStatus)then
       MpiComm = MpiComm_Global
       MpiSize = get_Size_MPI(MpiComm_Global)
       MpiRank = get_Rank_MPI(MpiComm_Global)
       mpiQ=0
       mpiR=0
       mpiIstart=0
       mpiIend=0
       mpiIshift=0
    endif
#endif
    !
  end subroutine delete_Hv_sector_superc






  

  function vecDim_Hv_sector_superc(isector) result(vecDim)
    integer :: isector
    integer :: Dim
    integer :: vecDim
    !
    Dim  = getdim(isector)
    !
#ifdef _MPI
    if(MpiStatus)then
       MpiQ = Dim/MpiSize
       MpiR = 0
       if(MpiRank==(MpiSize-1))MpiR=mod(Dim,MpiSize)
    else
       MpiQ = Dim
       MpiR = 0
    endif
#else
    MpiQ = Dim
    MpiR = 0
#endif
    !
    vecDim=MpiQ + MpiR
    !
  end function vecDim_Hv_sector_superc






  
  subroutine tridiag_Hv_sector_superc(isector,vvinit,alanc,blanc,norm2)
    integer                             :: isector
    complex(8),dimension(:)             :: vvinit
    real(8),dimension(:),allocatable    :: alanc,blanc
    real(8)                             :: norm2
    complex(8),dimension(:),allocatable :: vvloc
    integer                             :: vecDim
    !
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A)")&
         "DEBUG tridiag_Hv_sector SUPERC: start tridiag of H sector:"//str(isector)
#endif
    !
    if(MpiMaster)then
       norm2=dot_product(vvinit,vvinit)
       vvinit=vvinit/sqrt(norm2)
    endif
#ifdef _MPI
    if(MpiStatus)call bcast_MPI(MpiComm,norm2)
#endif
    call build_Hv_sector_superc(isector)
    allocate(alanc(Hsector%Nlanc),blanc(Hsector%Nlanc))
    alanc=0d0 ; blanc=0d0
    if(norm2/=0d0)then
#ifdef _MPI
       if(MpiStatus)then
          vecDim = vecDim_Hv_sector_superc(isector)
          allocate(vvloc(vecDim))
          call scatter_vector_MPI(MpiComm,vvinit,vvloc)
          call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvloc,alanc,blanc)
       else
          call sp_lanc_tridiag(spHtimesV_cc,vvinit,alanc,blanc)
       endif
#else
       call sp_lanc_tridiag(spHtimesV_cc,vvinit,alanc,blanc)
#endif
    endif
    call delete_Hv_sector_superc()
  end subroutine tridiag_Hv_sector_superc


end MODULE ED_HAMILTONIAN_SUPERC
