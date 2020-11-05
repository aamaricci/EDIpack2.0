MODULE ED_HAMILTONIAN_NONSU2
  USE ED_HAMILTONIAN_NONSU2_COMMON
  USE ED_HAMILTONIAN_NONSU2_STORED_HxV
  USE ED_HAMILTONIAN_NONSU2_DIRECT_HxV
  !
  implicit none
  private

  !>Build sparse hamiltonian of the sector
  public  :: build_Hv_sector_nonsu2
  public  :: delete_Hv_sector_nonsu2
  public  :: vecDim_Hv_sector_nonsu2

  !> Tridiag sparse Hamiltonian of the sector
  public  :: tridiag_Hv_sector_nonsu2




contains






  !####################################################################
  !                 MAIN ROUTINES: BUILD/DELETE SECTOR
  !####################################################################
  subroutine build_Hv_sector_nonsu2(isector,Hmat)
    integer                            :: isector,SectorDim
    complex(8),dimension(:,:),optional :: Hmat
    integer                            :: irank
    integer                            :: i,j,Dim
    !
    call build_sector(isector,Hsector)
    !
    Dim    = Hsector%Dim
    !
    !Total Split:
    MpiQ = Dim/MpiSize
    MpiR = 0
    if(MpiRank==(MpiSize-1))MpiR=mod(Dim,MpiSize)
    !
    MpiIshift = MpiRank*mpiQ
    MpiIstart = MpiRank*mpiQ + 1
    MpiIend   = (MpiRank+1)*mpiQ + mpiR
    !
#ifdef _MPI
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
    !
    !
    if(present(Hmat))then
       spHtimesV_cc => null()
       call ed_buildh_nonsu2_main(isector,Hmat)
       return
    endif
    !
    select case (ed_sparse_H)
    case (.true.)
       spHtimesV_cc => spMatVec_nonsu2_main
#ifdef _MPI
       if(MpiStatus)spHtimesV_cc => spMatVec_MPI_nonsu2_main
#endif
       call ed_buildh_nonsu2_main(isector)
       !
       !
    case (.false.)
       spHtimesV_cc => directMatVec_nonsu2_main
#ifdef _MPI
       if(MpiStatus)spHtimesV_cc => directMatVec_MPI_nonsu2_main
#endif
    end select
    !
  end subroutine build_Hv_sector_nonsu2


  subroutine delete_Hv_sector_nonsu2()
    integer :: iud
    call delete_sector(Hsector)
    Dim = 0
    !
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
  end subroutine delete_Hv_sector_nonsu2


  function vecDim_Hv_sector_nonsu2(isector) result(vecDim)
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
  end function vecDim_Hv_sector_nonsu2




  subroutine tridiag_Hv_sector_nonsu2(isector,vvinit,alanc,blanc,norm2)
    integer                             :: isector
    complex(8),dimension(:)             :: vvinit
    real(8),dimension(:),allocatable    :: alanc,blanc
    real(8)                             :: norm2
    complex(8),dimension(:),allocatable :: vvloc
    integer                             :: vecDim
    !
    if(MpiMaster)then
       norm2=dot_product(vvinit,vvinit)
       vvinit=vvinit/sqrt(norm2)
    endif
#ifdef _MPI
    if(MpiStatus)call bcast_MPI(MpiComm,norm2)
#endif
    call build_Hv_sector_nonsu2(isector)
    allocate(alanc(Hsector%Nlanc),blanc(Hsector%Nlanc))
    alanc=0d0 ; blanc=0d0
    if(norm2/=0d0)then
#ifdef _MPI
       if(MpiStatus)then
          vecDim = vecDim_Hv_sector_nonsu2(isector)
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
    call delete_Hv_sector_nonsu2()
  end subroutine tridiag_Hv_sector_nonsu2


end MODULE ED_HAMILTONIAN_NONSU2
