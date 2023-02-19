program ed_normal_normal
  USE DMFT_ED
  USE SCIFOR
  USE MPI
  USE SF_MPI
  USE ASSERTING
  implicit none
  integer                                     :: i,iw,jo,js,Nso
  integer                                     :: unit,unit_
  real(8)                                     :: w, Im, Re
  !Bath:
  integer                                     :: Nb,iorb,jorb,ispin,jspin,inso,print_mode
  real(8),allocatable                         :: Bath(:)
  !GFs and Sigma:
  complex(8),allocatable                      :: Weiss(:,:,:,:,:,:)
  complex(8),allocatable                      :: Smats(:,:,:,:,:,:)
  !hamiltonian input:
  complex(8),allocatable                      :: Hloc(:,:,:,:)
  real(8),dimension(:),allocatable            :: H0     ![Nso]
  !variables for the model:
  real(8)                                     :: Delta
  character(len=16)                           :: finput
  !NORMAL variables:
  real(8),allocatable                         :: dens(:), docc(:), energy(:), phisc(:), imp(:)
  complex(8),allocatable                      :: Smats11(:), Asmats11(:)
  !CHECK variables:
  real(8),allocatable                         :: dens_(:),docc_(:),energy_(:),phisc_(:),imp_(:)
  complex(8),allocatable                      :: Smats11_(:),ASmats11_(:)
  !
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
  if(bath_type/="normal")stop "Wrong setup from input file: non normal bath"
  if(ed_mode/="superc")stop "Wrong setup from input file: non superc ed_mode"
  if(Norb/=2)stop "Wrong setup from input file: Norb!=2"
  if(Nspin/=1 )stop "Wrong setup from input file: Nspin/=1"
  Nso=Nspin*Norb
  !Allocate Weiss Field:
  allocate(Weiss(2,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(2,Nspin,Nspin,Norb,Norb,Lmats))
  !
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
  Nb=ed_get_bath_dimension()
  allocate(Bath(Nb))
  call ed_init_solver(bath)
  !
  !
  !set Hloc
  call ed_set_Hloc(hloc)
  !
  !Solve the IMPURITY PROBLEM
  call ed_solve(bath)
  call ed_get_sigma(Smats(1,:,:,:,:,:),axis="m",type="n")
  call ed_get_sigma(Smats(2,:,:,:,:,:),axis="m",type="a")
  !
  !
  ! Check observables
  allocate(dens(Norb),dens_(Norb))
  allocate(docc(Norb),docc_(Norb))
  allocate(energy(8),energy_(8))
  allocate(phisc(4),phisc_(4))
  allocate(imp(4),imp_(4))
  allocate(Smats11(size(Smats,6)),Smats11_(size(Smats,6)))
  allocate(ASmats11(size(Smats,6)),ASmats11_(size(Smats,6)))
  ! density
  unit =free_unit()
  unit_=free_unit()
  open(unit,file="dens_last.ed")
  read(unit,*) dens(:)
  close(unit)
  open(unit_,file="dens_last.check")
  read(unit_,*) dens_(:)
  close(unit_)
  call assert(dens,dens_,"NORMAL_SUPERC dens(:)")

  open(unit,file="docc_last.ed")
  read(unit,*) docc(:)
  close(unit)
  open(unit_,file="docc_last.check")
  read(unit_,*) docc_(:)
  close(unit_)
  call assert(docc,docc_,"NORMAL_SUPERC docc(:)")

  open(unit,file="energy_last.ed")
  read(unit,*) energy(:)
  close(unit)
  open(unit_,file="energy_last.check")
  read(unit_,*) energy_(:)
  close(unit_)
  call assert(energy,energy_,"NORMAL_SUPERC energy(:)")

  open(unit,file="imp_last.ed")
  read(unit,*) imp(:)
  close(unit)
  open(unit_,file="imp_last.check")
  read(unit_,*) imp_(:)
  close(unit_)
  call assert(imp,imp_,"NORMAL_SUPERC imp(:)")

  open(unit,file="phisc_last.ed")
  read(unit,*) phisc(:)
  close(unit)
  open(unit_,file="phisc_last.check")
  read(unit_,*) phisc_(:)
  close(unit_)
  call assert(phisc,phisc_,"NORMAL_SUPERC phisc(:)")

  open(unit,file="impSigma_l11_s1_iw.ed")
  do iw=1,size(Smats,6)
     read(unit,*) w, Im, Re
     Smats11(iw)=Re+xi*Im
  enddo
  close(unit)
  open(unit_,file="impSigma_l11_s1_iw.check")
  do iw=1,size(Smats,6)
     read(unit_,*) w, Im, Re
     Smats11_(iw)=Re+xi*Im
  end do
  close(unit_)
  call assert(Smats11,Smats11_,"NORMAL_SUPERC Smats11(:)",tol=1.0d-8)

  open(unit,file="impSelf_l11_s1_iw.ed")
  do iw=1,size(Smats,6)
     read(unit,*) w, Im, Re
     ASmats11(iw)=Re+xi*Im
  enddo
  close(unit)
  open(unit_,file="impSelf_l11_s1_iw.check")
  do iw=1,size(Smats,6)
     read(unit_,*) w, Im, Re
     ASmats11_(iw)=Re+xi*Im
  end do
  close(unit_)
  call assert(ASmats11,ASmats11_,"NORMAL_SUPERC ASmats11(:)",tol=1.0d-8)

  
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

end program ed_normal_normal



