! > SPARSE MAT-VEC DIRECT ON-THE-FLY PRODUCT 
MODULE ED_HAMILTONIAN_NONSU2_DIRECT_HxV
  USE ED_HAMILTONIAN_NONSU2_COMMON
  implicit none
  private


  !>Sparse Mat-Vec direct on-the-fly product 
  public  :: directMatVec_nonsu2_main
#ifdef _MPI
  public  :: directMatVec_MPI_nonsu2_main
#endif




contains



  subroutine directMatVec_nonsu2_main(Nloc,vin,Hv)
    integer                                           :: Nloc
    complex(8),dimension(Nloc)                        :: vin
    complex(8),dimension(Nloc)                        :: Hv
    integer                                           :: isector
    integer,dimension(Nlevels)                        :: ib
    integer,dimension(Ns)                             :: ibup,ibdw
    real(8),dimension(Norb)                           :: nup,ndw
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Nbath) :: Hbath_tmp
    integer                                           :: first_state,last_state
    integer                                           :: first_state_up,last_state_up
    integer                                           :: first_state_dw,last_state_dw
    !
    if(.not.Hsector%status)stop "directMatVec_cc ERROR: Hsector NOT allocated"
    isector=Hsector%index
    !
    Dim = getdim(isector)
    !
    if(Nloc/=dim)stop "directMatVec_cc ERROR: Nloc != dim(isector)"
    !
    !Get diagonal hybridization, bath energy
    if(allocated(diag_hybr))deallocate(diag_hybr)
    if(allocated(bath_diag))deallocate(bath_diag)
    select case (bath_type)
    case default
       Nfoo = size(dmft_bath%e,2)
       allocate(diag_hybr(Nspin,Norb,Nbath));diag_hybr=0d0
       allocate(bath_diag(Nspin,Nfoo,Nbath));bath_diag=0d0       
       do ibath=1,Nbath
          do ispin=1,Nspin             
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%v(ispin,iorb,ibath)
             enddo
             do iorb=1,Nfoo
                bath_diag(ispin,iorb,ibath)=dmft_bath%e(ispin,iorb,ibath)
             enddo
          enddo
       enddo
    case ("replica")
       allocate(diag_hybr(Nspin,Norb,Nbath));diag_hybr=0d0
       allocate(bath_diag(Nspin,Norb,Nbath));bath_diag=0d0
       do ibath=1,Nbath
          Hbath_tmp(:,:,:,:,ibath) = Hreplica_build(dmft_bath%item(ibath)%lambda)
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%item(ibath)%v
                bath_diag(ispin,iorb,ibath)=Hbath_tmp(ispin,ispin,iorb,iorb,ibath)
             enddo
          enddo
       enddo
    end select
    !
    Hv=zero
    !-----------------------------------------------!
    states: do j=MpiIstart,MpiIend
       m    = Hsector%H(1)%map(j)
       ib   = bdecomp(m,2*Ns)
       !
       do iorb=1,Norb
          nup(iorb)=dble(ib(iorb))
          ndw(iorb)=dble(ib(iorb+Ns))
       enddo
       !
       !
       !IMPURITY  HAMILTONIAN
       include "direct/HxVimp.f90"
       !
       !LOCAL INTERACTION
       include "direct/HxVint.f90"
       !
       !BATH HAMILTONIAN
       include "direct/HxVbath.f90"
       !
       !IMPURITY- BATH HYBRIDIZATION
       include "direct/HxVimp_bath.f90"
    enddo states
    !-----------------------------------------------!
    !
  end subroutine directMatVec_nonsu2_main



#ifdef _MPI
  subroutine directMatVec_MPI_nonsu2_main(Nloc,v,Hv)
    integer                                           :: Nloc
    complex(8),dimension(Nloc)                        :: v
    complex(8),dimension(Nloc)                        :: Hv
    integer                                           :: N
    complex(8),dimension(:),allocatable               :: vin
    integer,allocatable,dimension(:)                  :: Counts,Offset
    integer                                           :: isector
    integer,dimension(Nlevels)                        :: ib
    integer,dimension(Ns)                             :: ibup,ibdw
    real(8),dimension(Norb)                           :: nup,ndw
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Nbath) :: Hbath_tmp
    integer                                           :: first_state,last_state
    integer                                           :: first_state_up,last_state_up
    integer                                           :: first_state_dw,last_state_dw
    integer                                           :: mpiIerr
    !
    if(MpiComm==MPI_UNDEFINED)stop "directMatVec_MPI_cc ERRROR: MpiComm = MPI_UNDEFINED"
    if(.not.MpiStatus)stop "directMatVec_MPI_cc ERROR: MpiStatus = F"
    !
    if(.not.Hsector%status)stop "directMatVec_cc ERROR: Hsector NOT allocated"
    isector=Hsector%index
    !
    Dim = getdim(isector)
    !
    !
    !Get diagonal hybridization, bath energy
    if(allocated(diag_hybr))deallocate(diag_hybr)
    if(allocated(bath_diag))deallocate(bath_diag)
    select case (bath_type)
    case default
       Nfoo = size(dmft_bath%e,2)
       allocate(diag_hybr(Nspin,Norb,Nbath));diag_hybr=0d0
       allocate(bath_diag(Nspin,Nfoo,Nbath));bath_diag=0d0       
       do ibath=1,Nbath
          do ispin=1,Nspin             
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%v(ispin,iorb,ibath)
             enddo
             do iorb=1,Nfoo
                bath_diag(ispin,iorb,ibath)=dmft_bath%e(ispin,iorb,ibath)
             enddo
          enddo
       enddo
    case ("replica")
       allocate(diag_hybr(Nspin,Norb,Nbath));diag_hybr=0d0
       allocate(bath_diag(Nspin,Norb,Nbath));bath_diag=0d0
       do ibath=1,Nbath
          Hbath_tmp(:,:,:,:,ibath) = Hreplica_build(dmft_bath%item(ibath)%lambda)
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%item(ibath)%v
                bath_diag(ispin,iorb,ibath)=Hbath_tmp(ispin,ispin,iorb,iorb,ibath)
             enddo
          enddo
       enddo
    end select
    !
    N=0
    call AllReduce_MPI(MpiComm,Nloc,N)
    !
    !Reconstruct Vin and get the displacements for AllGatherV call
    allocate(Counts(0:MpiSize-1)) ; Counts(0:)=0
    allocate(Offset(0:MpiSize-1)) ; Offset(0:)=0
    !
    Counts(0:)        = N/MpiSize
    Counts(MpiSize-1) = N/MpiSize+mod(N,MpiSize)
    !
    do i=1,MpiSize-1
       Offset(i) = Counts(i-1) + Offset(i-1)
    enddo
    !
    allocate(vin(N)); vin  = zero
    call MPI_Allgatherv(&
         v(1:Nloc),Nloc,MPI_Double_Complex,&
         vin,Counts,Offset,MPI_Double_Complex,&
         MpiComm,MpiIerr)
    !
    Hv=zero
    !
    !-----------------------------------------------!
    states: do j=MpiIstart,MpiIend
       m  = Hsector%H(1)%map(j)
       ib = bdecomp(m,2*Ns)
       !
       do iorb=1,Norb
          nup(iorb)=dble(ib(iorb))
          ndw(iorb)=dble(ib(iorb+Ns))
       enddo
       !
       !
       !IMPURITY  HAMILTONIAN
       include "direct/HxVimp.f90"
       !
       !LOCAL INTERACTION
       include "direct/HxVint.f90"
       !
       !BATH HAMILTONIAN
       include "direct/HxVbath.f90"
       !
       !IMPURITY- BATH HYBRIDIZATION
       include "direct/HxVimp_bath.f90"
    enddo states
    !-----------------------------------------------!
    !
  end subroutine directMatVec_MPI_nonsu2_main
#endif





end MODULE ED_HAMILTONIAN_NONSU2_DIRECT_HXV






