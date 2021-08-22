!-------------------------------------------------------------------!
! PURPOSE: INITIALIZE INTERNAL Hreplica STRUCTURES
!-------------------------------------------------------------------!

!allocate GLOBAL basis for H (used for impHloc and bath) and vectors coefficient
subroutine allocate_hreplica(N)
  integer          :: N
  integer          :: isym
  !
#ifdef _DEBUG
  if(ed_verbose>3)write(Logfile,"(A)")"DEBUG allocate_Hreplica"
#endif
  if(allocated(Hreplica_basis))deallocate(Hreplica_basis)
  if(allocated(Hreplica_lambda))deallocate(Hreplica_lambda)
  !
  allocate(Hreplica_basis(N))
  allocate(Hreplica_lambda(N))
  do isym=1,N
     allocate(Hreplica_basis(isym)%O(Nspin,Nspin,Norb,Norb))
     Hreplica_basis(isym)%O=0d0
     Hreplica_lambda(isym)=0d0
  enddo
  Hreplica_status=.true.
end subroutine allocate_hreplica


!deallocate GLOBAL basis for H (used for impHloc and bath) and vectors coefficient
subroutine deallocate_hreplica()
  integer              :: isym
  !
#ifdef _DEBUG
  if(ed_verbose>3)write(Logfile,"(A)")"DEBUG deallocate_Hreplica"
#endif
  do isym=1,size(Hreplica_basis)
     deallocate(Hreplica_basis(isym)%O)
  enddo
  deallocate(Hreplica_basis)
  deallocate(Hreplica_lambda)
  Hreplica_status=.false.
end subroutine deallocate_hreplica


!+------------------------------------------------------------------+
!PURPOSE  : Set Hreplica from user defined Hloc
!1: [Nspin,Nspin,Norb,Norb]
!2: [Nspin*Norb,Nspin*Norb]
!+------------------------------------------------------------------+
subroutine init_Hreplica_direct_nn(Hloc)
  integer                                     :: ispin,jspin,iorb,jorb,counter,io,jo,Nsym
  complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hloc
  logical(8),dimension(Nspin,Nspin,Norb,Norb) :: Hmask
  !
#ifdef _DEBUG
  if(ed_verbose>3)write(Logfile,"(A)")"DEBUG init_Hreplica_direct_nn: from Hloc[:,:,:,:]"
#endif
  !
  Hmask=.false.
  !
  counter=0
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              io=index_stride_so(ispin,iorb)
              jo=index_stride_so(jspin,jorb)
              if(io > jo )cycle
              if(Hloc(ispin,jspin,iorb,jorb)/=zero)counter=counter+1
           enddo
        enddo
     enddo
  enddo
  !
  call allocate_hreplica(counter)
  !
  counter=0
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              io=index_stride_so(ispin,iorb)
              jo=index_stride_so(jspin,jorb)
              if(io > jo )cycle
              if(Hloc(ispin,jspin,iorb,jorb)/=zero)then
                 counter=counter+1
                 Hreplica_basis(counter)%O(ispin,jspin,iorb,jorb)=1d0
                 Hreplica_basis(counter)%O(ispin,jspin,jorb,iorb)=1d0
                 Hreplica_lambda(counter)=Hloc(ispin,ispin,iorb,jorb)
              endif
           enddo
        enddo
     enddo
  enddo
  !
end subroutine init_Hreplica_direct_nn

subroutine init_Hreplica_direct_so(Hloc)
  complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hloc
#ifdef _DEBUG
  if(ed_verbose>3)write(Logfile,"(A)")"DEBUG init_Hreplica_direct_so: from Hloc[:,:]"
#endif
  call init_Hreplica_direct_nn(so2nn_reshape(Hloc,Nspin,Norb))
end subroutine init_Hreplica_direct_so


subroutine init_Hreplica_symmetries_site(Hvec,lambdavec)
  complex(8),dimension(:,:,:,:,:) :: Hvec
  real(8),dimension(:)            :: lambdavec
  integer                         :: isym,N
  !
#ifdef _DEBUG
  if(ed_verbose>3)write(Logfile,"(A)")"DEBUG init_Hreplica_direct_nn: from {[Hs,Lam]}_b"
#endif
  !
  N=size(lambdavec)
  call assert_shape(Hvec,[Nspin,Nspin,Norb,Norb,N],"init_Hreplica_symmetries","Hvec")
  !
  call allocate_hreplica(N)
  !
  do isym=1,N
     Hreplica_lambda(isym)  = lambdavec(isym)
     Hreplica_basis(isym)%O = Hvec(:,:,:,:,isym)
  enddo
  !
  if(ed_verbose>2)call print_hloc(Hreplica_build(Hreplica_lambda))
end subroutine init_Hreplica_symmetries_site



subroutine init_Hreplica_symmetries_lattice(Hvec,lambdavec)
  complex(8),dimension(:,:,:,:,:) :: Hvec
  real(8),dimension(:,:)          :: lambdavec ![Nlat,Nsym]
  integer                         :: isym,ilat,N,Nlat
  !
#ifdef _DEBUG
  if(ed_verbose>3)write(Logfile,"(A)")"DEBUG init_Hreplica_direct_nn: from ({[Hs,Lam]}_b)_site"
#endif
  !
  Nlat=size(lambdavec,1)
  N   =size(lambdavec,2)
  call assert_shape(Hvec,[Nspin,Nspin,Norb,Norb,N],"init_Hreplica_symmetries","Hvec")
  !
  if(allocated(Hreplica_lambda_ineq))deallocate(Hreplica_lambda_ineq)
  allocate(Hreplica_lambda_ineq(Nlat,N))
  call allocate_hreplica(N)
  !
  do isym=1,N
     Hreplica_lambda_ineq(:,isym)  = lambdavec(:,isym)
     Hreplica_basis(isym)%O = Hvec(:,:,:,:,isym)
  enddo
  !
  if(ed_verbose>2)then
     do ilat=1,Nlat
        call print_hloc(Hreplica_build(Hreplica_lambda_ineq(ilat,:)))
     enddo
  endif
end subroutine init_Hreplica_symmetries_lattice











