program ed_replica_normal
  USE EDIPACK2
  USE SCIFOR
  USE MPI
  USE SF_MPI
  USE ASSERTING
  implicit none
  integer                                     :: i,iw,jo,js,Nso,Nsymm,Nmomenta
  integer                                     :: unit,unit_
  real(8)                                     :: w,Re,Im
  !Bath:
  integer                                     :: Nb,iorb,jorb,ispin,jspin,inso,print_mode
  real(8),allocatable                         :: Bath(:),Wlist(:)
  !GFs and Sigma:
  complex(8),allocatable                      :: Weiss(:,:,:,:,:)
  complex(8),allocatable                      :: Smats(:,:,:,:,:)
  !hamiltonian input:
  complex(8),allocatable                      :: Hloc(:,:,:,:)
  !variables for the model:
  real(8)                                     :: mh,lambda
  character(len=30)                           :: Params
  character(len=16)                           :: finput
  !Replica variables:
  real(8),allocatable                         :: dens(:),docc(:),magX(:),energy(:),imp(:)
  real(8),allocatable                         :: Sab11s11mom(:),Sab12s11mom(:),Sab11s12mom(:),Sab12s12mom(:)
  !CHECK variables
  real(8),allocatable                         :: dens_(:),docc_(:),magX_(:),energy_(:),imp_(:)
  real(8),allocatable                         :: Sab11s11mom_(:),Sab12s11mom_(:),Sab11s12mom_(:),Sab12s12mom_(:)
  !
  complex(8),dimension(4,4)                   :: Gamma1,Gamma2,Gamma5,GammaN,GammaS
  complex(8),dimension(4,4)                   :: GammaE0,GammaEx,GammaEy,GammaEz
  real(8),dimension(:),allocatable            :: lambdasym_vector
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
  call parse_input_variable(mh,"MH",finput,default=1.d0)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.3d0)
  !
  !
  call ed_read_input(trim(finput))
  !
  if(ed_mode/='nonsu2')stop "Wrong setup from input file: ed_mode != nonsu2"
  if(bath_type/="normal")stop "Wrong setup from input file: non replica bath"
  if(Norb/=2)stop "Wrong setup from input file: Norb!=2"
  if(Nspin/=2 )stop "Wrong setup from input file: Nspin/=2"
  Nso=Nspin*Norb
  Nmomenta=4
  !Allocate Weiss Field:
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  !
  ! Matrices for replica hamiltonian
  gamma1=kron( pauli_sigma_z, pauli_tau_x)
  gamma2=kron( pauli_sigma_0,-pauli_tau_y)
  gamma5=kron( pauli_sigma_0, pauli_tau_z)
  gammaN=kron( pauli_sigma_0, pauli_tau_0)
  !
  gammaE0=kron( pauli_sigma_0, pauli_tau_x )
  gammaEx=kron( pauli_sigma_x, pauli_tau_x )
  gammaEy=kron( pauli_sigma_y, pauli_tau_x )
  gammaEz=kron( pauli_sigma_z, pauli_tau_x )


  Nb=ed_get_bath_dimension()
  !
  allocate(Bath(Nb))
  call ed_init_solver(bath)
  !
  !set Hloc
  allocate(Hloc(Nspin,Nspin,Norb,Norb))
  Hloc = j2so(Mh*Gamma5)
  call ed_set_Hloc(Hloc)
  !
  !Solve the IMPURITY PROBLEM
  call ed_solve(bath)
  call ed_get_sigma(Smats,axis="m")
  !


  ! !
  ! ! Check observables
  allocate(dens(Norb),dens_(Norb))
  allocate(docc(Norb),docc_(Norb))
  allocate(magX(Norb),magX_(Norb))
  allocate(energy(8),energy_(8))
  allocate(imp(2),imp_(2))
  allocate(Wlist(size(Smats,5)))
  allocate(Sab11s11mom(Nmomenta),Sab11s11mom_(Nmomenta))
  allocate(Sab11s12mom(Nmomenta),Sab11s12mom_(Nmomenta))
  write(*,*) ""
  write(*,*) "ED_MODE = NONSU2   |   BATH_TYPE = NORMAL"
  write(*,*) "Checking..."
  unit =free_unit()
  unit_=free_unit()


  ! density
  call ed_get_dens(dens)
  open(unit_,file="dens_last.check")
  read(unit_,*) dens_(:)
  close(unit_)
  call assert(dens,dens_,"dens(:)")
  !double occupancy
  call ed_get_docc(docc)
  open(unit_,file="docc_last.check")
  read(unit_,*) docc_(:)
  close(unit_)
  call assert(docc,docc_,"docc(:)")
  !magXYZ
  open(unit,file="magX_last.ed")
  read(unit,*) magX(:)
  close(unit)
  open(unit,file="magX_last.check")
  read(unit,*) magX_(:)
  close(unit)
  call assert(magX,magX_,"magX(:)")
  !energies
  open(unit,file="energy_last.ed")
  read(unit,*) energy(:)
  close(unit)
  open(unit_,file="energy_last.check")
  read(unit_,*) energy_(:)
  close(unit_)
  call assert(energy,energy_,"energy(:)")
  !impurity
  call ed_get_imp_info(imp)
  open(unit_,file="imp_last.check")
  read(unit_,*) imp_(:)
  close(unit_)
  call assert(imp,imp_,"imp(:)")


  !Self-Energies
  open(unit,file="impSigma_l11_s11_iw.ed")
  do iw=1,size(Smats,5)
     read(unit,*)Wlist(iw)
  end do
  close(unit)
  !
  ! Get momenta
  do i=1,Nmomenta
     call compute_momentum(Wlist,Smats(1,1,1,1,:),i,Sab11s11mom(i))
     call compute_momentum(Wlist,Smats(1,2,1,1,:),i,Sab11s12mom(i))
  enddo
  ! Write new momenta
  open(unit_,file="impSigma_l11_s11_iw.momenta.new")
  do i=1,Nmomenta
     write(unit_,*) i, Sab11s11mom(i)
  enddo
  close(unit_)
  open(unit_,file="impSigma_l11_s12_iw.momenta.new")
  do i=1,Nmomenta
     write(unit_,*) i, Sab11s12mom(i)
  enddo
  close(unit_)
  ! Read check momenta
  open(unit_,file="impSigma_l11_s11_iw.momenta.check")
  do i=1,Nmomenta
     read(unit_,*) iw, Sab11s11mom_(i)
  enddo
  close(unit_)
  open(unit_,file="impSigma_l11_s12_iw.momenta.check")
  do i=1,Nmomenta
     read(unit_,*) iw, Sab11s12mom_(i)
  enddo
  close(unit_)
  call assert(Sab11s11mom/Sab11s11mom_,dble(ones(Nmomenta)),"Sigma_matsubara_l11_s11(:)",tol=1.0d-8)
  call assert(Sab11s12mom/Sab11s12mom_,dble(ones(Nmomenta)),"Sigma_matsubara_l11_s12(:)",tol=1.0d-8)


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

  ! Subroutine to compute momenta
  ! 
  ! ( sum_w abs(F(w))*w**n ) / ( sum_w abs(F(w)) )
  subroutine compute_momentum(x,Fx,n,momentum)
    real(8)   ,dimension(:),intent(in)       :: x
    complex(8),dimension(:),intent(in)       :: Fx
    integer   ,intent(in)                    :: n
    real(8)   ,intent(out)                   :: momentum
    !
    integer                                  :: iw
    real(8)                                  :: num,den
    !
    num=0.0;    den=0.0
    ! num = sum(abs(Fx*x**n))
    ! den = sum(abs(Fx))
    do iw=1,size(x)
       num = num + (abs(Fx(iw)) )*x(iw)**n
       den = den + (abs(Fx(iw)) )
    enddo
    momentum=num/den
  end subroutine compute_momentum

end program ed_replica_normal



