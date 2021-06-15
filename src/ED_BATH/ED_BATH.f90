MODULE ED_BATH
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,txtfy
  USE SF_LINALG, only: eye,inv
  USE SF_ARRAYS, only:linspace
  USE SF_MISC, only: assert_shape
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  implicit none

  private

  !##################################################################
  !
  !     USER BATH ROUTINES:
  !
  !##################################################################
  !Interface for user bath I/O operations: get,set,copy
  interface get_bath_dimension
     module procedure ::  get_bath_dimension_direct
     module procedure ::  get_bath_dimension_symmetries
  end interface get_bath_dimension

  interface set_Hreplica
     module procedure init_Hreplica_direct_so
     module procedure init_Hreplica_direct_nn
     module procedure init_Hreplica_symmetries_site
     module procedure init_Hreplica_symmetries_lattice
  end interface set_Hreplica

  !explicit symmetries:
  interface break_symmetry_bath
     module procedure break_symmetry_bath_site
     module procedure break_symmetry_bath_lattice
  end interface break_symmetry_bath

  interface spin_symmetrize_bath
     module procedure spin_symmetrize_bath_site
     module procedure spin_symmetrize_bath_lattice
  end interface spin_symmetrize_bath

  interface orb_symmetrize_bath
     module procedure orb_symmetrize_bath_site
     module procedure orb_symmetrize_bath_lattice
  end interface orb_symmetrize_bath

  interface orb_equality_bath
     module procedure orb_equality_bath_site
     module procedure orb_equality_bath_lattice
  end interface orb_equality_bath

  interface ph_symmetrize_bath
     module procedure ph_symmetrize_bath_site
     module procedure ph_symmetrize_bath_lattice
  end interface ph_symmetrize_bath

  interface ph_trans_bath
     module procedure ph_trans_bath_site
     module procedure ph_trans_bath_lattice
  end interface ph_trans_bath

  interface enforce_normal_bath
     module procedure enforce_normal_bath_site
     module procedure enforce_normal_bath_lattice
  end interface enforce_normal_bath

  !Aux:
  interface get_Whyb_matrix
     module procedure get_Whyb_matrix_1orb
     module procedure get_Whyb_matrix_Aorb
     module procedure get_Whyb_matrix_dmft_bath
  end interface get_Whyb_matrix

  interface is_identity
     module procedure ::  is_identity_so
     module procedure ::  is_identity_nn
  end interface is_identity

  interface is_diagonal
     module procedure ::  is_diagonal_so
     module procedure ::  is_diagonal_nn
  end interface is_diagonal






  !##################################################################
  !
  !     USER BATH ROUTINES:
  !
  !##################################################################
  public :: get_bath_dimension
  public :: check_bath_dimension
  !
  public :: get_bath_component_dimension
  public :: get_bath_component
  public :: set_bath_component
  public :: copy_bath_component
  !

  public :: break_symmetry_bath              
  public :: spin_symmetrize_bath
  public :: orb_symmetrize_bath
  public :: orb_equality_bath
  public :: ph_symmetrize_bath
  public :: ph_trans_bath
  public :: enforce_normal_bath
  public :: get_Whyb_matrix
  public :: impose_equal_lambda
  public :: impose_bath_offset
  !
  public :: set_Hreplica


  !##################################################################
  !
  !     DMFT BATH ROUTINES:
  !
  !##################################################################
  public :: allocate_dmft_bath               !INTERNAL (for effective_bath)
  public :: deallocate_dmft_bath             !INTERNAL (for effective_bath)
  public :: init_dmft_bath                   !INTERNAL (for effective_bath)
  public :: write_dmft_bath                  !INTERNAL (for effective_bath)
  public :: save_dmft_bath                   !INTERNAL (for effective_bath)
  public :: set_dmft_bath                    !INTERNAL (for effective_bath)
  public :: get_dmft_bath                    !INTERNAL (for effective_bath)
  !
  public :: hreplica_build                   !INTERNAL (for effective_bath)
  public :: hreplica_mask                    !INTERNAL (for effective_bath)
  public :: hreplica_site                    !INTERNAL (for effective_bath)


  integer :: ibath,ilat,iorb



contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : Inquire the correct bath size to allocate the 
  ! the bath array in the calling program.
  !+-------------------------------------------------------------------+
  function get_bath_dimension_direct(Hloc_nn) result(bath_size)
    complex(8),optional,intent(in) :: Hloc_nn(:,:,:,:)
    integer                        :: bath_size,ndx,ispin,iorb,jspin,jorb,io,jo,Maxspin
    complex(8),allocatable         :: Hloc(:,:,:,:)
    !
    select case(bath_type)
       !
    case default
       select case(ed_mode)
       case default
          bath_size = Norb*Nbath + Norb*Nbath
          !( e [Nspin][Norb][Nbath] + v [Nspin][Norb][Nbath] )
       case ("superc")
          bath_size = Norb*Nbath + Norb*Nbath + Norb*Nbath
          !( e [Nspin][Norb][Nbath] + d [Nspin][Norb][Nbath] + v [Nspin][Norb][Nbath] )
       case ("nonsu2")
          bath_size = Norb*Nbath + Norb*Nbath + Norb*Nbath
          !( e [Nspin][Norb][Nbath] + v [Nspin][Norb][Nbath] + u [Nspin][Norb][Nbath] )
       end select
       bath_size=Nspin*bath_size
       !
    case('hybrid')
       select case(ed_mode)
       case default          
          bath_size = Nbath + Norb*Nbath
          !(e [Nspin][1][Nbath] + v [Nspin][Norb][Nbath] )
       case ("superc")
          bath_size = Nbath + Nbath + Norb*Nbath
          !(e [Nspin][1][Nbath] + d [Nspin][1][Nbath] + v [Nspin][Norb][Nbath] )
       case ("nonsu2")
          bath_size = Nbath + Norb*Nbath + Norb*Nbath
          !(e [Nspin][1][Nbath] + v [Nspin][Norb][Nbath] + u [Nspin][Norb][Nbath] )
       end select
       bath_size=Nspin*bath_size
       !
    case('replica')
       allocate(Hloc(Nspin,Nspin,Norb,Norb))
       if(present(Hloc_nn))then              !User defined Hloc_nn 
          Hloc=Hloc_nn
       elseif(Hreplica_status)then !User defined Hreplica_basis
          Hloc=Hreplica_build(Hreplica_lambda)
       else                                  !Error:
          deallocate(Hloc)
          stop "ERROR get_bath_dimension_direct: ed_mode=replica neither Hloc_nn present nor Hreplica_basis defined"
       endif
       !
       !Check Hermiticity:
       do ispin=1,Nspin
          do iorb=1,Norb
             if(abs(dimag(Hloc(ispin,ispin,iorb,iorb))).gt.1d-6)stop "Hloc is not Hermitian"
          enddo
       enddo
       !
       !Re/Im off-diagonal non-vanishing elements
       ndx=0
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   if(io > jo)cycle
                   if(dreal(Hloc(ispin,jspin,iorb,jorb)) /= 0d0)ndx=ndx+1
                   if(dimag(Hloc(ispin,jspin,iorb,jorb)) /= 0d0)ndx=ndx+1
                enddo
             enddo
          enddo
       enddo
       !
       ndx = ndx * Nbath !number of non vanishing elements for each replica
       ndx = ndx + Nbath !diagonal hybridizations: Vs (different per spin)
       ndx = ndx + 1     !we also print Nbasis
       bath_size = ndx
    end select
  end function get_bath_dimension_direct

  function get_bath_dimension_symmetries(Hloc_nn) result(bath_size)
    complex(8),dimension(:,:,:,:,:),intent(in) :: Hloc_nn
    integer                                    :: bath_size,ndx,isym,Nsym
    !
    !number of symmetries
    Nsym=size(Hloc_nn,5)
    if(Nsym/=size(Hreplica_lambda))stop "ERROR get_bath_dimension_symmetries:  neither Hloc_nn present nor Hreplica_basis defined"
    !
    ndx=Nsym
    !
    !number of replicas
    ndx = ndx * Nbath
    !diagonal hybridizations: Vs
    ndx = ndx + Nbath
    !
    !include Nbasis
    ndx=ndx+1
    !
    bath_size = ndx
    !
  end function get_bath_dimension_symmetries



  !+-------------------------------------------------------------------+
  !PURPOSE  : Check if the dimension of the bath array are consistent
  !+-------------------------------------------------------------------+
  function check_bath_dimension(bath_) result(bool)
    real(8),dimension(:)           :: bath_
    integer                        :: Ntrue,i
    logical                        :: bool
    complex(8),allocatable         :: Hreplica(:,:,:,:,:)![Nspin][:][Norb][:][Nsym]
    select case (bath_type)
    case default
       Ntrue = get_bath_dimension()
    case ('replica')
       if(.not.Hreplica_status)STOP "check_bath_dimension: Hreplica_basis not allocated"
       if(.not.allocated(Hreplica))allocate(Hreplica(Nspin,Nspin,Norb,Norb,size(Hreplica_basis)))
       do i=1,size(Hreplica_basis)
          Hreplica(:,:,:,:,i)=Hreplica_basis(i)%O
       enddo
       Ntrue   = get_bath_dimension_symmetries(Hreplica)
    end select
    bool  = ( size(bath_) == Ntrue )
  end function check_bath_dimension





  !##################################################################
  !
  !     USER BATH ROUTINES:
  !
  !##################################################################
  include 'user_aux.f90'



  !##################################################################
  !
  !     DMFT BATH ROUTINES:
  !
  !##################################################################
  include 'dmft_aux.f90'


  !##################################################################
  !
  !     H_REPLICA ROUTINES:
  !
  !##################################################################
  include 'hreplica_setup.f90'




  !##################################################################
  !
  !     AUX FUNCTIONS:
  !
  !##################################################################
  subroutine Hreplica_site(site)
    integer :: site
    if(site<1.OR.site>size(Hreplica_lambda_ineq,1))stop "ERROR Hreplica_site: site not in [1,Nlat]"
    if(.not.allocated(Hreplica_lambda_ineq))stop "ERROR Hreplica_site: Hreplica_lambda_ineq not allocated"
    Hreplica_lambda(:)  = Hreplica_lambda_ineq(site,:)
  end subroutine Hreplica_site


  !reconstruct [Nspin,Nspin,Norb,Norb] hamiltonian from basis expansion given [lambda]
  function Hreplica_build(lambdavec) result(H)
    real(8),dimension(:)                        :: lambdavec
    integer                                     :: isym
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: H
    !
    if(.not.Hreplica_status)STOP "ERROR Hreplica_build: Hreplica_basis is not setup"
    if(size(lambdavec)/=size(Hreplica_basis)) STOP "ERROR Hreplica_build: Wrong coefficient vector size"
    H=zero
    do isym=1,size(lambdavec)
       H=H+lambdavec(isym)*Hreplica_basis(isym)%O
    enddo
  end function Hreplica_build



  !+-------------------------------------------------------------------+
  !PURPOSE  : Create bath mask
  !+-------------------------------------------------------------------+
  function Hreplica_mask(wdiag,uplo) result(Hmask)
    logical,optional                            :: wdiag,uplo
    logical                                     :: wdiag_,uplo_
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hloc
    logical,dimension(Nspin,Nspin,Norb,Norb)    :: Hmask
    integer                                     :: iorb,jorb,ispin,jspin,io,jo
    !
    wdiag_=.false.;if(present(wdiag))wdiag_=wdiag
    uplo_ =.false.;if(present(uplo))  uplo_=uplo
    !
    Hloc = Hreplica_build(Hreplica_lambda)
    Hmask=.false.
    where(abs(Hloc)>1d-6)Hmask=.true.
    !
    !
    if(wdiag_)then
       do ispin=1,Nspin
          do iorb=1,Norb
             Hmask(ispin,ispin,iorb,iorb)=.true.
          enddo
       enddo
    endif
    !
    if(uplo_)then
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = index_stride_so(ispin,iorb)
                   jo = index_stride_so(jspin,jorb)
                   if(io>jo)Hmask(ispin,jspin,iorb,jorb)=.false.
                enddo
             enddo
          enddo
       enddo
    endif
    !
  end function Hreplica_mask




  !+-------------------------------------------------------------------+
  !PURPOSE  : Check if a matrix is the identity
  !+-------------------------------------------------------------------+
  function is_identity_nn(mnnn) result(flag)
    complex(8),dimension(nspin,nspin,norb,norb) :: mnnn
    real(8),dimension(nspin*norb,nspin*norb) :: mtmp
    integer                                  :: i,j
    logical                                  :: flag
    !
    flag=.true.
    !
    mtmp=nn2so_reshape(mnnn,nspin,norb)
    !
    do i=1,Nspin*Norb-1
       if((mtmp(i,i).ne.mtmp(i+1,i+1)).or.(mtmp(i,i).lt.1.d-6))flag=.false.
    enddo
    !
    do i=1,Nspin*Norb
       do j=1,Nspin*Norb
          if((i.ne.j).and.(mtmp(i,j).gt.1.d-6))flag=.false.
       enddo
    enddo
    !
  end function is_identity_nn

  function is_identity_so(mlso) result(flag)
    complex(8),dimension(nspin*norb,nspin*norb) :: mlso
    real(8),dimension(nspin*norb,nspin*norb) :: mtmp
    integer                                  :: i,j
    logical                                  :: flag
    !
    flag=.true.
    !
    mtmp=mlso
    !
    do i=1,Nspin*Norb-1
       if((mtmp(i,i).ne.mtmp(i+1,i+1)).or.(mtmp(i,i).lt.1.d-6))flag=.false.
    enddo
    !
    do i=1,Nspin*Norb
       do j=1,Nspin*Norb
          if((i.ne.j).and.(mtmp(i,j).gt.1.d-6))flag=.false.
       enddo
    enddo
    !
  end function is_identity_so



  !+-------------------------------------------------------------------+
  !PURPOSE  : Check if a matrix is diagonal
  !+-------------------------------------------------------------------+
  function is_diagonal_nn(mnnn) result(flag)
    complex(8),dimension(nspin,nspin,norb,norb) :: mnnn
    complex(8),dimension(nspin*norb,nspin*norb) :: mtmp
    integer                                     :: i,j
    logical                                     :: flag
    !
    flag=.true.
    !
    mtmp=abs((nn2so_reshape(mnnn,nspin,norb)))
    !
    do i=1,Nspin*Norb
       do j=1,Nspin*Norb
          if((i.ne.j).and.(abs(mtmp(i,j)).gt.1.d-6))flag=.false.
       enddo
    enddo
    !
  end function is_diagonal_nn

  function is_diagonal_so(mlso) result(flag)
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: mlso
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: mtmp
    integer                                     :: i,j
    logical                                     :: flag
    !
    flag=.true.
    !
    mtmp=abs((mlso))
    !
    do i=1,Nspin*Norb
       do j=1,Nspin*Norb
          if((i.ne.j).and.(abs(mtmp(i,j)).gt.1.d-6))flag=.false.
       enddo
    enddo
    !
  end function is_diagonal_so










END MODULE ED_BATH
