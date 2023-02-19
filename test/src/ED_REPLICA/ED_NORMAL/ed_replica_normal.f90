program ed_replica_normal
  USE DMFT_ED
  USE SCIFOR
  USE MPI
  USE SF_MPI
  USE ASSERTING
  implicit none
  integer                                     :: i,iw,jo,js,Nso,Nsymm
  integer                                     :: unit,unit_
  real(8)                                     :: w,Re,Im
  !Bath:
  integer                                     :: Nb,iorb,jorb,ispin,jspin,inso,print_mode
  real(8),allocatable                         :: Bath(:)
  !GFs and Sigma:
  complex(8),allocatable                      :: Weiss(:,:,:,:,:,:)
  complex(8),allocatable                      :: Smats(:,:,:,:,:)
  !hamiltonian input:
  complex(8),allocatable                      :: Hloc(:,:,:,:)
  real(8),dimension(:),allocatable            :: H0     ![Nso]
  !variables for the model:
  real(8)                                     :: Delta
  character(len=16)                           :: finput
  !Replica variables:
  real(8),allocatable                         :: dens(:),docc(:),exciton(:),energy(:),imp(:)
  complex(8),allocatable                      :: Smats11(:), Smats12(:)
  !CHECK variables
  real(8),allocatable                         :: dens_(:),docc_(:),exciton_(:),energy_(:),imp_(:)
  complex(8),allocatable                      :: Smats11_(:),Smats12_(:)
  !
  complex(8),dimension(4,4)                   :: GammaN,GammaE0
  real(8),dimension(:,:),allocatable          :: lambdasym_vector
  complex(8),dimension(:,:,:,:,:),allocatable :: Hsym_basis
  !MPI Vars:
  integer                                     :: irank,comm,rank,size2,ierr
  logical                                     :: master
  !
  ! MPI initialization
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  size2 = get_Size_MPI(comm)
  master = get_Master_MPI(comm)
  !
  !Parse additional variables && read Input
  call parse_cmd_variable(finput,"FINPUT",default="inputED.in")
  call parse_input_variable(delta,"DELTA",finput,default=0.d0)
  !
  !
  call ed_read_input(trim(finput))
  !
  if(bath_type/="replica")stop "Wrong setup from input file: non replica bath"
  if(Norb/=2)stop "Wrong setup from input file: Norb!=2"
  if(Nspin/=1 )stop "Wrong setup from input file: Nspin/=1"
  Nso=Nspin*Norb
  !Allocate Weiss Field:
  allocate(Weiss(1,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  !
  ! Matrices for replica hamiltonian
  gammaN =kron_pauli( pauli_sigma_0, pauli_tau_0)
  gammaE0=kron_pauli( pauli_sigma_0, pauli_tau_x )
  !
  allocate(Hloc(Nspin,Nspin,Norb,Norb))
  allocate(H0(Nso))
  Hloc = zero
  H0   = zero
  do js=1,Nspin
     Hloc(js,js,:,:)= Delta*pauli_sigma_z
     do jo=1,Norb
        H0(jo+2*(js-1)) =Hloc(js,js,jo,jo)
     end do
  end do
  !
  print_mode=3
  !
  ! Set up replica hamiltonian
  Nsymm=2
  allocate(lambdasym_vector(Nbath,Nsymm))
  allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,Nsymm))
  !
  ! N
  Hsym_basis(:,:,:,:,1)=j2so(GammaN(:Nso,:Nso))
  do i=1,Nbath
     lambdasym_vector(i,1) = -1.0 + 2.0*dble(i-1)/dble(Nbath-1)
  end do
  !
  ! E0
  Hsym_basis(:,:,:,:,2)=j2so(GammaE0(:Nso,:Nso))
  lambdasym_vector(:,2)=0.1d0
  !
  call ed_set_Hreplica(Hsym_basis,lambdasym_vector)
  Nb=ed_get_bath_dimension(Nsymm)
  allocate(Bath(Nb))
  call ed_init_solver(bath)
  !
  !
  !set Hloc
  call ed_set_Hloc(hloc)
  !
  !Solve the IMPURITY PROBLEM
  call ed_solve(bath)
  call ed_get_sigma(Smats,axis="m",type="n")
  !
  !
  ! Check observables
  allocate(dens(Norb),dens_(Norb))
  allocate(docc(Norb),docc_(Norb))
  allocate(exciton(2),exciton_(2))
  allocate(energy(8),energy_(8))
  allocate(imp(4),imp_(4))
  allocate(Smats11(size(Smats,1)),Smats11_(size(Smats,1)))
  allocate(Smats12(size(Smats,1)),Smats12_(size(Smats,1)))
  ! density
  unit =free_unit()
  unit_=free_unit()
  open(unit,file="dens_last.ed")
  read(unit,*) dens(:)
  close(unit)
  open(unit_,file="dens_last.check")
  read(unit_,*) dens_(:)
  close(unit_)
  call assert(dens,dens_,"REPLICA_NORMAL dens(:)")
  
  open(unit,file="docc_last.ed")
  read(unit,*) docc(:)
  close(unit)
  open(unit_,file="docc_last.check")
  read(unit_,*) docc_(:)
  close(unit_)
  call assert(docc,docc_,"REPLICA_NORMAL docc(:)")
  
  open(unit,file="exciton_last.ed")
  read(unit,*) exciton(:)
  close(unit)
  open(unit_,file="exciton_last.check")
  read(unit_,*) exciton_(:)
  close(unit_)
  call assert(exciton,exciton_,"REPLICA_NORMAL exciton(:)")

  open(unit,file="energy_last.ed")
  read(unit,*) energy(:)
  close(unit)
  open(unit_,file="energy_last.check")
  read(unit_,*) energy_(:)
  close(unit_)
  call assert(energy,energy_,"REPLICA_NORMAL energy(:)")
  
  open(unit,file="imp_last.ed")
  read(unit,*) imp(:)
  close(unit)
  open(unit_,file="imp_last.check")
  read(unit_,*) imp_(:)
  close(unit_)
  call assert(imp,imp_,"REPLICA_NORMAL imp(:)")
  
  open(unit,file="impSigma_l11_s1_iw.ed")
  do iw=1,size(Smats,1)
     read(unit,*) w, Im, Re
     Smats11(iw) = Re+xi*Im
  end do
  close(unit)
  open(unit_,file="impSigma_l11_s1_iw.check")
  do iw=1,size(Smats,1)
     read(unit_,*) w, Im, Re
     Smats11_(iw) = Re+xi*Im
  end do
  close(unit_)
  call assert(Smats11,Smats11_,"REPLICA_NORMAL Smats11(:)",tol=1.0d-8)

  open(unit,file="impSigma_l12_s1_iw.ed")
  do iw=1,size(Smats,1)
     read(unit,*) w, Im, Re
     Smats12(iw) = Re+xi*Im
  end do
  close(unit)
  open(unit_,file="impSigma_l12_s1_iw.check")
  do iw=1,size(Smats,1)
     read(unit_,*) w, Im, Re
     Smats12_(iw) = Re+xi*Im
  end do
  close(unit_)
  call assert(Smats12,Smats12_,"REPLICA_NORMAL Smats12(:)",tol=1.0d-8)
  
  call finalize_MPI()



contains


  function so2j_index(ispin,iorb) result(isporb)
    integer :: ispin,iorb
    integer :: isporb
    if(iorb>Norb)stop "error so2j_index: iorb>Norb"
    if(ispin>Nspin)stop "error so2j_index: ispin>Nspin"
    isporb=(ispin-1)*Nspin + iorb
  end function so2j_index


  function so2j(fg) result(g)
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: fg
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: g
    integer                                     :: i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i=so2j_index(ispin,iorb)
                j=so2j_index(jspin,jorb)
                g(i,j) = fg(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function so2j

  function j2so(fg) result(g)
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: fg
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: g
    integer                                     :: i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i=so2j_index(ispin,iorb)
                j=so2j_index(jspin,jorb)
                g(ispin,jspin,iorb,jorb)  = fg(i,j)
             enddo
          enddo
       enddo
    enddo
  end function j2so

end program ed_replica_normal



