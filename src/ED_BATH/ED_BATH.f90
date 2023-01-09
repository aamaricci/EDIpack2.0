MODULE ED_BATH
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,str
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
     module procedure init_Hreplica_symmetries_site
     module procedure init_Hreplica_symmetries_legacy  ! (deprecation-cycle)
     module procedure init_Hreplica_symmetries_lattice
  end interface set_Hreplica

  interface set_Hgeneral
     module procedure init_Hgeneral_symmetries_site
     module procedure init_Hgeneral_symmetries_legacy  ! (deprecation-cycle)
     module procedure init_Hgeneral_symmetries_lattice
  end interface set_Hgeneral

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
     module procedure orb_symmetrize_bath_site_o1o2
     module procedure orb_symmetrize_bath_lattice_o1o2
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
  public :: set_Hgeneral

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
  !
  public :: hgeneral_build                   !INTERNAL (for effective_bath)
  public :: hgeneral_mask                    !INTERNAL (for effective_bath)
  public :: hgeneral_site                    !INTERNAL (for effective_bath)


  integer :: ibath,ilat,iorb



contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : Inquire the correct bath size to allocate the
  ! the bath array in the calling program.
  !+-------------------------------------------------------------------+
  function get_bath_dimension_direct(H_nn) result(bath_size)
    complex(8),optional,intent(in) :: H_nn(:,:,:,:)
    integer                        :: bath_size,ndx,ispin,iorb,jspin,jorb,io,jo,Maxspin
    complex(8),allocatable         :: H(:,:,:,:)
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
       allocate(H(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb))
       if(present(H_nn))then    !User defined H_nn
          H=H_nn
       elseif(Hreplica_status)then !User defined Hreplica_basis
          H=Hreplica_build(Hreplica_lambda(Nbath,:))
       else                        !Error:
          deallocate(H)
          stop "ERROR get_bath_dimension_direct: ed_mode=replica neither H_nn present nor Hreplica_basis defined"
       endif
       !
       !Check Hermiticity:
       ! do ispin=1,Nspin
       !    do iorb=1,Norb
       !       if(abs(dimag(H(ispin,ispin,iorb,iorb))).gt.1d-6)stop "H is not Hermitian"
       !    enddo
       ! enddo
       if( all(abs(nn2so_reshape(H,Nnambu*Nspin,Norb) - conjg(transpose(nn2so_reshape(H,Nnambu*Nspin,Norb))))<1d-6)  )stop "H is not Hermitian"
       !
       !Re/Im off-diagonal non-vanishing elements
       ndx=0
       do ispin=1,Nnambu*Nspin
          do jspin=1,Nnambu*Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   if(io > jo)cycle
                   if(dreal(H(ispin,jspin,iorb,jorb)) /= 0d0)ndx=ndx+1
                   if(dimag(H(ispin,jspin,iorb,jorb)) /= 0d0)ndx=ndx+1
                enddo
             enddo
          enddo
       enddo
       !
       ndx = ndx * Nbath !number of non vanishing elements for each replica
       ndx = ndx + Nbath !diagonal hybridizations: Vs (different per spin)
       ndx = ndx + 1     !we also print Nbasis
       bath_size = ndx
    case('general')
       allocate(H(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb))
       if(present(H_nn))then    !User defined H_nn
          H=H_nn
       elseif(Hgeneral_status)then !User defined Hgeneral_basis
          H=Hgeneral_build(Hgeneral_lambda(Nbath,:))
       else                        !Error:
          deallocate(H)
          stop "ERROR get_bath_dimension_direct: ed_mode=general neither H_nn present nor Hgeneral_basis defined"
       endif
       !
       !Check Hermiticity:
       ! do ispin=1,Nspin
       !    do iorb=1,Norb
       !       if(abs(dimag(H(ispin,ispin,iorb,iorb))).gt.1d-6)stop "H is not Hermitian"
       !    enddo
       ! enddo
       if( all(abs(nn2so_reshape(H,Nnambu*Nspin,Norb) - conjg(transpose(nn2so_reshape(H,Nnambu*Nspin,Norb))))<1d-6)  )stop "H is not Hermitian"
       !
       !Re/Im off-diagonal non-vanishing elements
       ndx=0
       do ispin=1,Nnambu*Nspin
          do jspin=1,Nnambu*Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   if(io > jo)cycle
                   if(dreal(H(ispin,jspin,iorb,jorb)) /= 0d0)ndx=ndx+1
                   if(dimag(H(ispin,jspin,iorb,jorb)) /= 0d0)ndx=ndx+1
                enddo
             enddo
          enddo
       enddo
       !
       ndx = ndx * Nbath            !number of non vanishing elements for each general
       ndx = ndx + Nbath*Norb*Nspin !diagonal hybridizations: Vs (different per spin and orbitals)
       ndx = ndx + 1                !we also print Nbasis
       bath_size = ndx
    end select
  end function get_bath_dimension_direct

  function get_bath_dimension_symmetries(Nsym) result(bath_size)
    integer :: bath_size,ndx,isym,Nsym
    !
    select case(bath_type)
    case("replica")
       if(.not.Hreplica_status)STOP "get_bath_dimension_symmetries: H(replica/general)_basis  not allocated"
       if(Nsym/=size(Hreplica_lambda,2))&
            stop "ERROR get_bath_dimension_symmetries:  size(Hreplica_basis) != size(Hreplica_lambda,2)"
    case("general")
       if(.not.Hgeneral_status)STOP "get_bath_dimension_symmetries: H(general/general)_basis  not allocated"
       if(Nsym/=size(Hgeneral_lambda,2))&
            stop "ERROR get_bath_dimension_symmetries:  size(Hgeneral_basis) != size(Hgeneral_lambda,2)"
    case default
       stop "ERROR get_bath_dimension_symmetris wiht bath_type!=replica/general"
    end select
    !
    ndx=Nsym
    !
    !number of replicas
    ndx = ndx * Nbath
    !diagonal hybridizations: Vs
    select case(bath_type)
    case("replica") ! Vk depends only on bath site
       ndx = ndx + Nbath
    case("general")
       ndx = ndx + Nbath*Norb*Nspin
    end select
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
    !complex(8),allocatable         :: Hreplica(:,:,:,:,:)![Nspin][:][Norb][:][Nsym]
    select case (bath_type)
    case default
       Ntrue = get_bath_dimension()
    case ('replica')
       Ntrue   = get_bath_dimension_symmetries(size(Hreplica_basis))
    case ('general')
       Ntrue   = get_bath_dimension_symmetries(size(Hgeneral_basis))
    end select
    bool  = ( size(bath_) == Ntrue )
  end function check_bath_dimension




















  !##################################################################
  !
  !     USER BATH  SYMMETRIES: PREDEFINED AND USER CONTROLLED
  !
  !##################################################################

  !+-------------------------------------------------------------------+
  !PURPOSE  : given a bath array apply a specific transformation or
  ! impose a given symmetry:
  ! - break spin symmetry by applying a symmetry breaking field
  ! - given a bath array set both spin components to have
  !    the same bath, i.e. impose non-magnetic solution
  ! - given a bath array enforces the particle-hole symmetry
  !    by setting the positive energies in modulo identical to the negative
  !    ones.
  ! - given a bath enforce normal (i.e. non superconducting) solution
  ! - given a dmft bath pull/push the components W^{ss'}_\a(l) of the Hybridization
  !    matrix
  ! - given a dmft bath pull/push the nonsu2 components
  !+-------------------------------------------------------------------+
  subroutine impose_equal_lambda(bath_,ibath,lambdaindex_vec)
    real(8),dimension(:) :: bath_
    type(effective_bath) :: dmft_bath_
    real(8)              :: val
    integer,dimension(:) :: lambdaindex_vec
    integer              :: i,N,ibath
    !
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    !
    N=size(lambdaindex_vec)
    val=0.d0
    do i=1,N
       val=val+dmft_bath_%item(ibath)%lambda(lambdaindex_vec(i))/N
    enddo
    !
    do i=1,N
       dmft_bath_%item(ibath)%lambda(lambdaindex_vec(i))=val
    enddo
    !
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine impose_equal_lambda


  subroutine impose_bath_offset(bath_,ibath,offset)
    real(8),dimension(:) :: bath_
    type(effective_bath) :: dmft_bath_
    real(8)              :: offset
    integer              :: isym,N,ibath
    !
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    !
    select case(bath_type)
    case default
       if(size(Hreplica_basis) .ne. dmft_bath_%Nbasis)then
          dmft_bath_%item(ibath)%lambda(dmft_bath_%Nbasis)=offset
       else
          do isym=1,size(Hreplica_basis)
             if(is_identity(Hreplica_basis(isym)%O)) dmft_bath_%item(ibath)%lambda(isym)=offset
             return
          enddo
       endif
    case("general")
       if(size(Hgeneral_basis) .ne. dmft_bath_%Nbasis)then
          dmft_bath_%item(ibath)%lambda(dmft_bath_%Nbasis)=offset
       else
          do isym=1,size(Hgeneral_basis)
             if(is_identity(Hgeneral_basis(isym)%O)) dmft_bath_%item(ibath)%lambda(isym)=offset
             return
          enddo
       endif
    end select
       
    !
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
    !
  end subroutine impose_bath_offset


  subroutine break_symmetry_bath_site(bath_,field,sign,save)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    real(8)                :: field
    real(8)                :: sign
    logical,optional       :: save
    logical                :: save_
    if(bath_type=="replica")stop "break_symmetry_bath_site ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "break_symmetry_bath_site ERROR: can not be used with bath_type=general"
    save_=.true.;if(present(save))save_=save
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    dmft_bath_%e(1,:,:)    =dmft_bath_%e(1,:,:)      + sign*field
    dmft_bath_%e(Nspin,:,:)=dmft_bath_%e(Nspin,:,:)  - sign*field
    if(save_)call save_dmft_bath(dmft_bath_)
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine break_symmetry_bath_site
  !
  subroutine break_symmetry_bath_lattice(bath_,field,sign,save)
    real(8),dimension(:,:) :: bath_
    real(8)                :: field
    real(8)                :: sign
    logical,optional       :: save
    logical                :: save_
    integer                :: Nsites,ilat
    save_=.true.;if(present(save))save_=save
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix=reg(ineq_site_suffix)//reg(str(ilat,site_indx_padding))
       call break_symmetry_bath_site(bath_(ilat,:),field,sign,save_)
    enddo
    ed_file_suffix=""
  end subroutine break_symmetry_bath_lattice

  !---------------------------------------------------------!


  subroutine spin_symmetrize_bath_site(bath_,save)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    logical,optional       :: save
    logical                :: save_
    integer :: ibath
    if(bath_type=="replica")stop "spin_symmetry_bath_site ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "spin_symmetry_bath_site ERROR: can not be used with bath_type=general"
    save_=.true.;if(present(save))save_=save
    if(Nspin==1)then
       write(LOGfile,"(A)")"spin_symmetrize_bath: Nspin=1 nothing to symmetrize"
       return
    endif
    !
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    select case(ed_mode)
    case default
       dmft_bath_%e(Nspin,:,:)=dmft_bath_%e(1,:,:)
       dmft_bath_%v(Nspin,:,:)=dmft_bath_%v(1,:,:)
    case ("superc")
       dmft_bath_%e(Nspin,:,:)=dmft_bath_%e(1,:,:)
       dmft_bath_%v(Nspin,:,:)=dmft_bath_%v(1,:,:)
       dmft_bath_%d(Nspin,:,:)=dmft_bath_%d(1,:,:)
    end select
    if(save_)call save_dmft_bath(dmft_bath_)
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine spin_symmetrize_bath_site
  !
  subroutine spin_symmetrize_bath_lattice(bath_,save)
    real(8),dimension(:,:) :: bath_
    logical,optional       :: save
    logical                :: save_
    integer                :: Nsites,ilat
    save_=.true.;if(present(save))save_=save
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix=reg(ineq_site_suffix)//reg(str(ilat,site_indx_padding))
       call spin_symmetrize_bath_site(bath_(ilat,:),save_)
    enddo
    ed_file_suffix=""
  end subroutine spin_symmetrize_bath_lattice

  !---------------------------------------------------------!

  subroutine orb_symmetrize_bath_site(bath_,save)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    logical,optional       :: save
    logical                :: save_
    integer                :: iorb
    real(8),allocatable    :: lvl(:,:),hyb(:,:)
    if(bath_type=="replica")stop "orb_symmetry_bath_site ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "orb_symmetry_bath_site ERROR: can not be used with bath_type=general"
    save_=.true.;if(present(save))save_=save
    if(Norb==1)then
       write(LOGfile,"(A)")"orb_symmetrize_bath: Norb=1 nothing to symmetrize"
       return
    endif
    !
    call allocate_dmft_bath(dmft_bath_)
    ! if (bath_type=="replica")call init_dmft_bath_mask(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    !
    if(allocated(lvl))deallocate(lvl);allocate(lvl(Nspin,Nbath));lvl=0d0;lvl=sum(dmft_bath_%e,dim=2)/Norb
    if(allocated(hyb))deallocate(hyb);allocate(hyb(Nspin,Nbath));hyb=0d0;hyb=sum(dmft_bath_%v,dim=2)/Norb
    do iorb=1,Norb
       dmft_bath_%e(:,iorb,:)=lvl
       dmft_bath_%v(:,iorb,:)=hyb
    enddo
    !
    if(save_)call save_dmft_bath(dmft_bath_)
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine orb_symmetrize_bath_site
  subroutine orb_symmetrize_bath_lattice(bath_,save)
    real(8),dimension(:,:) :: bath_
    logical,optional       :: save
    logical                :: save_
    integer                :: Nsites,ilat
    if(bath_type=="replica")stop "orb_symmetry_bath_site ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "orb_symmetry_bath_site ERROR: can not be used with bath_type=general"
    save_=.true.;if(present(save))save_=save
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix=reg(ineq_site_suffix)//reg(str(ilat,site_indx_padding))
       call orb_symmetrize_bath_site(bath_(ilat,:),save_)
    enddo
    ed_file_suffix=""
  end subroutine orb_symmetrize_bath_lattice


  subroutine orb_symmetrize_bath_site_o1o2(bath_,orb1,orb2,save)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    logical,optional       :: save
    logical                :: save_
    integer                :: iorb,orb1,orb2
    real(8),allocatable    :: lvl(:,:),hyb(:,:)
    if(bath_type=="replica")stop "orb_symmetry_bath_site ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "orb_symmetry_bath_site ERROR: can not be used with bath_type=general"
    save_=.true.;if(present(save))save_=save
    if(Norb==1)then
       write(LOGfile,"(A)")"orb_symmetrize_bath: Norb=1 nothing to symmetrize"
       return
    endif
    !
    call allocate_dmft_bath(dmft_bath_)
    ! if (bath_type=="replica")call init_dmft_bath_mask(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    !
    if(allocated(lvl))deallocate(lvl);allocate(lvl(Nspin,Nbath));lvl=0d0
    if(allocated(hyb))deallocate(hyb);allocate(hyb(Nspin,Nbath));hyb=0d0
    !
    lvl=(dmft_bath_%e(:,orb1,:)+dmft_bath_%e(:,orb2,:))/2d0
    hyb=(dmft_bath_%v(:,orb1,:)+dmft_bath_%v(:,orb2,:))/2d0
    !
    dmft_bath_%e(:,orb1,:)=lvl
    dmft_bath_%v(:,orb1,:)=hyb
    dmft_bath_%e(:,orb2,:)=lvl
    dmft_bath_%v(:,orb2,:)=hyb
    !
    if(save_)call save_dmft_bath(dmft_bath_)
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine orb_symmetrize_bath_site_o1o2
  subroutine orb_symmetrize_bath_lattice_o1o2(bath_,orb1,orb2,save)
    real(8),dimension(:,:) :: bath_
    logical,optional       :: save
    logical                :: save_
    integer                :: Nsites,ilat,orb1,orb2
    if(bath_type=="replica")stop "orb_symmetry_bath_site ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "orb_symmetry_bath_site ERROR: can not be used with bath_type=general"
    save_=.true.;if(present(save))save_=save
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix=reg(ineq_site_suffix)//reg(str(ilat,site_indx_padding))
       call orb_symmetrize_bath_site_o1o2(bath_(ilat,:),orb1,orb2,save_)
    enddo
    ed_file_suffix=""
  end subroutine orb_symmetrize_bath_lattice_o1o2

  !---------------------------------------------------------!


  subroutine orb_equality_bath_site(bath_,indx,save)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    integer,optional       :: indx
    logical,optional       :: save
    integer                :: indx_
    logical                :: save_
    integer                :: iorb
    real(8),allocatable    :: lvl(:,:),hyb(:,:)
    if(bath_type=="replica")stop "orb_equality_bath_site ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "orb_equality_bath_site ERROR: can not be used with bath_type=general"
    indx_=1     ;if(present(indx))indx_=indx
    save_=.true.;if(present(save))save_=save
    if(Norb==1)then
       write(LOGfile,"(A)")"orb_symmetrize_bath: Norb=1 nothing to symmetrize"
       return
    endif
    !
    call allocate_dmft_bath(dmft_bath_)
    ! if (bath_type=="replica")call init_dmft_bath_mask(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    !
    if(allocated(lvl))deallocate(lvl);allocate(lvl(Nspin,Nbath));lvl=0d0;lvl=dmft_bath_%e(:,indx_,:)
    if(allocated(hyb))deallocate(hyb);allocate(hyb(Nspin,Nbath));hyb=0d0;hyb=dmft_bath_%v(:,indx_,:)
    do iorb=1,Norb
       if(iorb==indx_)cycle
       dmft_bath_%e(:,iorb,:)=lvl
       dmft_bath_%v(:,iorb,:)=hyb
    enddo
    !
    if(save_)call save_dmft_bath(dmft_bath_)
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine orb_equality_bath_site
  subroutine orb_equality_bath_lattice(bath_,indx,save)
    real(8),dimension(:,:) :: bath_
    integer,optional       :: indx
    logical,optional       :: save
    integer                :: indx_
    logical                :: save_
    integer                :: iorb
    integer                :: Nsites,ilat
    if(bath_type=="replica")stop "orb_equality_bath_site ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "orb_equality_bath_site ERROR: can not be used with bath_type=general"
    indx_=1     ;if(present(indx))indx_=indx
    save_=.true.;if(present(save))save_=save
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix=reg(ineq_site_suffix)//reg(str(ilat,site_indx_padding))
       call orb_equality_bath_site(bath_(ilat,:),indx_,save_)
    enddo
    ed_file_suffix=""
  end subroutine orb_equality_bath_lattice



  !---------------------------------------------------------!

  subroutine ph_symmetrize_bath_site(bath_,save)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    integer                :: i
    logical,optional       :: save
    logical                :: save_
    if(bath_type=="replica")stop "ph_symmetry_bath_site ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "ph_symmetry_bath_site ERROR: can not be used with bath_type=general"
    save_=.true.;if(present(save))save_=save
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    if(Nbath==1)return
    if(mod(Nbath,2)==0)then
       do i=1,Nbath/2
          dmft_bath_%e(:,:,Nbath+1-i)=-dmft_bath_%e(:,:,i)
          dmft_bath_%v(:,:,Nbath+1-i)= dmft_bath_%v(:,:,i)
          if(ed_mode=="superc")dmft_bath_%d(:,:,Nbath+1-i)=dmft_bath_%d(:,:,i)
       enddo
    else
       do i=1,(Nbath-1)/2
          dmft_bath_%e(:,:,Nbath+1-i)=-dmft_bath_%e(:,:,i)
          dmft_bath_%v(:,:,Nbath+1-i)= dmft_bath_%v(:,:,i)
          if(ed_mode=="superc")dmft_bath_%d(:,:,Nbath+1-i)=dmft_bath_%d(:,:,i)
       enddo
       dmft_bath_%e(:,:,(Nbath-1)/2+1)=0.d0
    endif
    if(save_)call save_dmft_bath(dmft_bath_)
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine ph_symmetrize_bath_site
  subroutine ph_symmetrize_bath_lattice(bath_,save)
    real(8),dimension(:,:) :: bath_
    logical,optional       :: save
    logical                :: save_
    integer                :: Nsites,ilat
    if(bath_type=="replica")stop "ph_symmetry_bath_site ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "ph_symmetry_bath_site ERROR: can not be used with bath_type=general"
    save_=.true.;if(present(save))save_=save
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix=reg(ineq_site_suffix)//reg(str(ilat,site_indx_padding))
       call ph_symmetrize_bath_site(bath_(ilat,:),save_)
    enddo
    ed_file_suffix=""
  end subroutine ph_symmetrize_bath_lattice

  !---------------------------------------------------------!

  subroutine ph_trans_bath_site(bath_,save)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    type(effective_bath)   :: tmp_dmft_bath
    integer                :: i
    logical,optional       :: save
    logical                :: save_
    if(bath_type=="replica")stop "ph_trans_bath_site ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "ph_trans_bath_site ERROR: can not be used with bath_type=general"
    save_=.true.;if(present(save))save_=save
    call allocate_dmft_bath(dmft_bath_)
    call allocate_dmft_bath(tmp_dmft_bath)
    call set_dmft_bath(bath_,dmft_bath_)
    if(Nbath==1)return
    do i=1,Nbath
       select case(Norb)
       case default
          ! do nothing
          dmft_bath_%e(:,:,i)= dmft_bath_%e(:,:,i)
          dmft_bath_%v(:,:,i)= dmft_bath_%v(:,:,i)
       case(1)
          dmft_bath_%e(:,:,i)= -dmft_bath_%e(:,:,i)
          dmft_bath_%v(:,:,i)=  dmft_bath_%v(:,:,i)
       case(2)
          tmp_dmft_bath%e(:,1,i) = -dmft_bath_%e(:,2,i)
          tmp_dmft_bath%e(:,2,i) = -dmft_bath_%e(:,1,i)
          dmft_bath_%e(:,:,i)    = tmp_dmft_bath%e(:,:,i)
          tmp_dmft_bath%v(:,1,i) = dmft_bath_%v(:,2,i)
          tmp_dmft_bath%v(:,2,i) = dmft_bath_%v(:,1,i)
          dmft_bath_%v(:,:,i)    = tmp_dmft_bath%v(:,:,i)
       end select
    end do
    if(save_)call save_dmft_bath(dmft_bath_)
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine ph_trans_bath_site
  subroutine ph_trans_bath_lattice(bath_,save)
    real(8),dimension(:,:) :: bath_
    logical,optional       :: save
    logical                :: save_
    integer                :: Nsites,ilat
    if(bath_type=="replica")stop "ph_trans_bath_site ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "ph_trans_bath_site ERROR: can not be used with bath_type=general"
    save_=.true.;if(present(save))save_=save
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix=reg(ineq_site_suffix)//reg(str(ilat,site_indx_padding))
       call ph_trans_bath_site(bath_(ilat,:),save_)
    enddo
    ed_file_suffix=""
  end subroutine ph_trans_bath_lattice

  !---------------------------------------------------------!

  subroutine enforce_normal_bath_site(bath_,save)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    logical,optional       :: save
    logical                :: save_
    if(bath_type=="replica")stop "enforce_normal_bath_site ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "enforce_normal_bath_site ERROR: can not be used with bath_type=general"
    save_=.true.;if(present(save))save_=save
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    if(ed_mode=="superc")dmft_bath_%d(:,:,:)=0.d0
    if(save_)call save_dmft_bath(dmft_bath_)
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine enforce_normal_bath_site
  subroutine enforce_normal_bath_lattice(bath_,save)
    real(8),dimension(:,:) :: bath_
    logical,optional       :: save
    logical                :: save_
    integer                :: Nsites,ilat
    save_=.true.;if(present(save))save_=save
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix=reg(ineq_site_suffix)//reg(str(ilat,site_indx_padding))
       call enforce_normal_bath_site(bath_(ilat,:),save_)
    enddo
    ed_file_suffix=""
  end subroutine enforce_normal_bath_lattice










  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  check if the specified itype is consistent with the input parameters.
  !+-----------------------------------------------------------------------------+!
  subroutine check_bath_component(type)
    character(len=1) :: type
    select case(bath_type)
    case default
       select case(ed_mode)
       case default
          if(type/="e".AND.type/='v')stop "check_bath_component error: type!=e,v"
       case ("superc")
          if(type/="e".AND.type/='v'.AND.type/='d')stop "check_bath_component error: type!=e,v,d"
       case ("nonsu2")
          if(type/="e".AND.type/='v'.AND.type/='u')stop "check_bath_component error: type!=e,v,u"
       end select
    case ("replica","general")
       if(type/="v".AND.type/="l")stop "check_bath_component error: type!=v,l"
    end select
    return
  end subroutine check_bath_component































  !##################################################################
  !
  !     DMFT BATH ROUTINES:
  !
  !##################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE  : Allocate the ED bath
  !+-------------------------------------------------------------------+
  subroutine allocate_dmft_bath(dmft_bath_)
    type(effective_bath) :: dmft_bath_
    integer              :: Nsym,ibath
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG allocate_dmft_bath"
#endif
    if(dmft_bath_%status)call deallocate_dmft_bath(dmft_bath_)
    !
    select case(bath_type)
    case default
       !
       select case(ed_mode)
       case default                                 !normal [N,Sz]
          allocate(dmft_bath_%e(Nspin,Norb,Nbath))  !local energies of the bath
          allocate(dmft_bath_%v(Nspin,Norb,Nbath))  !same-spin hybridization
       case ("superc")                              !superc [Sz]
          allocate(dmft_bath_%e(Nspin,Norb,Nbath))  !local energies of the bath
          allocate(dmft_bath_%d(Nspin,Norb,Nbath))  !local SC order parameters the bath
          allocate(dmft_bath_%v(Nspin,Norb,Nbath))  !same-spin hybridization
       case ("nonsu2")                              !nonsu2 [N]
          allocate(dmft_bath_%e(Nspin,Norb,Nbath))  !local energies of the bath
          allocate(dmft_bath_%v(Nspin,Norb,Nbath))  !same-spin hybridization
          allocate(dmft_bath_%u(Nspin,Norb,Nbath))  !spin-flip hybridization
       end select
       !
    case('hybrid')
       !
       select case(ed_mode)
       case default                                 !normal  [N,Sz]
          allocate(dmft_bath_%e(Nspin,1,Nbath))     !local energies of the bath
          allocate(dmft_bath_%v(Nspin,Norb,Nbath))  !same-spin hybridization
       case ("superc")                              !superc  [Sz]
          allocate(dmft_bath_%e(Nspin,1,Nbath))     !local energies of the bath
          allocate(dmft_bath_%d(Nspin,1,Nbath))     !local SC order parameters the bath
          allocate(dmft_bath_%v(Nspin,Norb,Nbath))  !same-spin hybridization
       case ("nonsu2")                              !nonsu2 case [N] qn
          allocate(dmft_bath_%e(Nspin,1,Nbath))     !local energies of the bath
          allocate(dmft_bath_%v(Nspin,Norb,Nbath))  !same-spin hybridization
          allocate(dmft_bath_%u(Nspin,Norb,Nbath))  !spin-flip hybridization
       end select
       !
    case('replica')
       !
       if(.not.Hreplica_status)stop "ERROR allocate_dmft_bath: Hreplica_basis not allocated"
       call deallocate_dmft_bath(dmft_bath_)     !
       Nsym=size(Hreplica_basis)
       !
       allocate(dmft_bath_%item(Nbath))
       dmft_Bath_%Nbasis=Nsym
       do ibath=1,Nbath
          allocate(dmft_bath_%item(ibath)%lambda(Nsym))
       enddo
       !
    case('general')
       !
       if(.not.Hgeneral_status)stop "ERROR allocate_dmft_bath: Hgeneral_basis not allocated"
       call deallocate_dmft_bath(dmft_bath_)     !
       Nsym=size(Hgeneral_basis)
       !
       allocate(dmft_bath_%item(Nbath))
       dmft_Bath_%Nbasis=Nsym
       do ibath=1,Nbath
          allocate(dmft_bath_%item(ibath)%lambda(Nsym))
          allocate(dmft_bath_%item(ibath)%vg(Norb*Nspin))
       enddo
       !
    end select
    !
    dmft_bath_%status=.true.
    !
  end subroutine allocate_dmft_bath


  !+-------------------------------------------------------------------+
  !PURPOSE  : Deallocate the ED bath
  !+-------------------------------------------------------------------+
  subroutine deallocate_dmft_bath(dmft_bath_)
    type(effective_bath) :: dmft_bath_
    integer              :: ibath,isym
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG deallocate_dmft_bath"
#endif
    if(.not.dmft_bath_%status)return
    if(allocated(dmft_bath_%e))   deallocate(dmft_bath_%e)
    if(allocated(dmft_bath_%d))   deallocate(dmft_bath_%d)
    if(allocated(dmft_bath_%v))   deallocate(dmft_bath_%v)
    if(allocated(dmft_bath_%u))   deallocate(dmft_bath_%u)
    if(bath_type=="replica")then
       dmft_bath_%Nbasis= 0
       do ibath=1,Nbath
          deallocate(dmft_bath_%item(ibath)%lambda)
       enddo
       deallocate(dmft_bath_%item)
    endif
    if(bath_type=="general")then
       dmft_bath_%Nbasis= 0
       do ibath=1,Nbath
          deallocate(dmft_bath_%item(ibath)%lambda)
          deallocate(dmft_bath_%item(ibath)%vg)
       enddo
       deallocate(dmft_bath_%item)
    endif
    dmft_bath_%status=.false.
  end subroutine deallocate_dmft_bath






  !+------------------------------------------------------------------+
  !PURPOSE  : Initialize the DMFT loop, builindg H parameters and/or
  !reading previous (converged) solution
  !+------------------------------------------------------------------+
  subroutine init_dmft_bath(dmft_bath_,used)
    type(effective_bath) :: dmft_bath_
    logical,optional     :: used
    integer              :: Nbasis
    integer              :: i,ibath,isym,unit,flen,Nh,Nsym
    integer              :: io,jo,iorb,ispin,jorb,jspin
    logical              :: IOfile,used_,diagonal_hsym,all_lambdas_are_equal
    real(8)              :: de
    real(8)              :: offset(Nbath)
    real(8)              :: one_lambdaval
    character(len=20)    :: hsuffix
    !
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG init_dmft_bath"
#endif
    used_   = .false.   ;if(present(used))used_=used
    hsuffix = ".restart";if(used_)hsuffix=reg(".used")
    if(.not.dmft_bath_%status)stop "ERROR init_dmft_bath error: bath not allocated"
    !
    select case(bath_type)
    case default
       !Get energies:
       dmft_bath_%e(:,:,1)    =-ed_hw_bath
       dmft_bath_%e(:,:,Nbath)= ed_hw_bath
       Nh=Nbath/2
       if(mod(Nbath,2)==0.and.Nbath>=4)then
          de=ed_hw_bath/max(Nh-1,1)
          dmft_bath_%e(:,:,Nh)  = -1.d-1
          dmft_bath_%e(:,:,Nh+1)=  1.d-1
          do i=2,Nh-1
             dmft_bath_%e(:,:,i)   =-ed_hw_bath + (i-1)*de
             dmft_bath_%e(:,:,Nbath-i+1)= ed_hw_bath - (i-1)*de
          enddo
       elseif(mod(Nbath,2)/=0.and.Nbath>=3)then
          de=ed_hw_bath/Nh
          dmft_bath_%e(:,:,Nh+1)= 0d0
          do i=2,Nh
             dmft_bath_%e(:,:,i)        =-ed_hw_bath + (i-1)*de
             dmft_bath_%e(:,:,Nbath-i+1)= ed_hw_bath - (i-1)*de
          enddo
       endif
       !Get spin-keep yhbridizations
       do i=1,Nbath
          dmft_bath_%v(:,:,i)=max(0.1d0,1d0/sqrt(dble(Nbath)))
       enddo
       !Get SC amplitudes
       if(ed_mode=="superc")dmft_bath_%d(:,:,:)=deltasc
       !Get spin-flip hybridizations
       if(ed_mode=="nonsu2")then
          do i=1,Nbath
             dmft_bath_%u(:,:,i) = dmft_bath_%v(:,:,i)!*ed_vsf_ratio
          enddo
       endif
       !
    case('replica')
       offset=0.d0
       if(Nbath>1) offset=linspace(-ed_offset_bath,ed_offset_bath,Nbath)
       !
       !BATH V INITIALIZATION
       do ibath=1,Nbath
          dmft_bath%item(ibath)%v=max(0.1d0,1d0/sqrt(dble(Nbath)))
       enddo
       !
       !BATH LAMBDAS INITIALIZATION
       !Do not need to check for Hreplica_basis: this is done at allocation time of the dmft_bath.
       Nsym = dmft_bath%Nbasis
       do isym=1,Nsym
          do ibath=1,Nbath
             dmft_bath%item(ibath)%lambda(isym) =  Hreplica_lambda(ibath,isym)
          enddo
          diagonal_hsym = is_diagonal(Hreplica_basis(isym)%O)
          one_lambdaval = Hreplica_lambda(Nbath,isym)
          all_lambdas_are_equal = all(Hreplica_lambda(:,isym)==one_lambdaval)
          if(diagonal_hsym.AND.all_lambdas_are_equal)then
             offset=linspace(-ed_offset_bath,ed_offset_bath,Nbath)
             if(is_identity(Hreplica_basis(isym)%O).AND.mod(Nbath,2)==0)then
                offset(Nbath/2) = max(-1.d-1,offset(Nbath/2))
                offset(Nbath/2 + 1) = min(1.d-1,offset(Nbath/2 + 1))
             endif
             do ibath=1,Nbath
                dmft_bath%item(ibath)%lambda(isym) =  Hreplica_lambda(ibath,isym) + offset(ibath)
             enddo
             write(*,*) "                                                                    "
             write(*,*) "WARNING: some of your lambdasym values have been internally changed "
             write(*,*) "         while calling ed_init_solver. This happens whenever the    "
             write(*,*) "         corresponding Hsym is diagonal and all the replicas receive"
             write(*,*) "         the same initial lambda value, due to the deprecated legacy"
             write(*,*) "         practice of defining a unique lambda vector forall replicas"
             write(*,*) "         and let the solver decide how to handle these degeneracies."
             write(*,*) "         >>> If you really intend to have a degenerate diagonal term"
             write(*,*) "             in the bath you can define a suitable restart file.    "
             write(*,*) "         >>> If instead this is what you expected please consider to"
             write(*,*) "             move the desired rescaling in your driver, since this  "
             write(*,*) "             funcionality might be removed in a future update.      "
             write(*,*) "                                                                    "
          endif
       enddo
       !
    case('general')
       offset=0.d0
       if(Nbath>1) offset=linspace(-ed_offset_bath,ed_offset_bath,Nbath)
       !
       !BATH V INITIALIZATION
       do ibath=1,Nbath
          dmft_bath%item(ibath)%vg(:)=max(0.1d0,1d0/sqrt(dble(Nbath)))
       enddo
       !
       !BATH LAMBDAS INITIALIZATION
       !Do not need to check for Hgeneral_basis: this is done at allocation time of the dmft_bath.
       Nsym = dmft_bath%Nbasis
       do isym=1,Nsym
          do ibath=1,Nbath
             dmft_bath%item(ibath)%lambda(isym) =  Hgeneral_lambda(ibath,isym)
          enddo
          diagonal_hsym = is_diagonal(Hgeneral_basis(isym)%O)
          one_lambdaval = Hgeneral_lambda(Nbath,isym)
          all_lambdas_are_equal = all(Hgeneral_lambda(:,isym)==one_lambdaval)
          if(diagonal_hsym.AND.all_lambdas_are_equal)then
             offset=linspace(-ed_offset_bath,ed_offset_bath,Nbath)
             if(is_identity(Hgeneral_basis(isym)%O).AND.mod(Nbath,2)==0)then
                offset(Nbath/2) = max(-1.d-1,offset(Nbath/2))
                offset(Nbath/2 + 1) = min(1.d-1,offset(Nbath/2 + 1))
             endif
             do ibath=1,Nbath
                dmft_bath%item(ibath)%lambda(isym) =  Hgeneral_lambda(ibath,isym) + offset(ibath)
             enddo
             write(*,*) "                                                                    "
             write(*,*) "WARNING: some of your lambdasym values have been internally changed "
             write(*,*) "         while calling ed_init_solver. This happens whenever the    "
             write(*,*) "         corresponding Hsym is diagonal and all the generals receive"
             write(*,*) "         the same initial lambda value, due to the deprecated legacy"
             write(*,*) "         practice of defining a unique lambda vector forall generals"
             write(*,*) "         and let the solver decide how to handle these degeneracies."
             write(*,*) "         >>> If you really intend to have a degenerate diagonal term"
             write(*,*) "             in the bath you can define a suitable restart file.    "
             write(*,*) "         >>> If instead this is what you expected please consider to"
             write(*,*) "             move the desired rescaling in your driver, since this  "
             write(*,*) "             funcionality might be removed in a future update.      "
             write(*,*) "                                                                    "
          endif
       enddo
       !
    end select
    !
    !
    !
    !Read from file if exist:
    inquire(file=trim(Hfile)//trim(ed_file_suffix)//trim(hsuffix),exist=IOfile)
    if(IOfile)then
       write(LOGfile,"(A)")'Reading bath from file'//trim(Hfile)//trim(ed_file_suffix)//trim(hsuffix)
       unit = free_unit()
       flen = file_length(trim(Hfile)//trim(ed_file_suffix)//trim(hsuffix))
       !
       open(unit,file=trim(Hfile)//trim(ed_file_suffix)//trim(hsuffix))
       !
       select case(bath_type)
       case default
          !
          read(unit,*)
          select case(ed_mode)
          case default
             do i=1,min(flen,Nbath)
                read(unit,*)((&
                     dmft_bath_%e(ispin,iorb,i),&
                     dmft_bath_%v(ispin,iorb,i),&
                     iorb=1,Norb),ispin=1,Nspin)
             enddo
          case ("superc")
             do i=1,min(flen,Nbath)
                read(unit,*)((&
                     dmft_bath_%e(ispin,iorb,i),&
                     dmft_bath_%d(ispin,iorb,i),&
                     dmft_bath_%v(ispin,iorb,i),&
                     iorb=1,Norb),ispin=1,Nspin)
             enddo
          case("nonsu2")
             do i=1,min(flen,Nbath)
                read(unit,*)((&
                     dmft_bath_%e(ispin,iorb,i),&
                     dmft_bath_%v(ispin,iorb,i),&
                     dmft_bath_%u(ispin,iorb,i),&
                     iorb=1,Norb),ispin=1,Nspin)
             enddo
          end select
          !
       case ('hybrid')
          read(unit,*)
          !
          select case(ed_mode)
          case default
             do i=1,min(flen,Nbath)
                read(unit,*)(&
                     dmft_bath_%e(ispin,1,i),&
                     (&
                     dmft_bath_%v(ispin,iorb,i),&
                     iorb=1,Norb),&
                     ispin=1,Nspin)
             enddo
          case ("superc")
             do i=1,min(flen,Nbath)
                read(unit,*)(&
                     dmft_bath_%e(ispin,1,i),&
                     dmft_bath_%d(ispin,1,i),&
                     (&
                     dmft_bath_%v(ispin,iorb,i),&
                     iorb=1,Norb),&
                     ispin=1,Nspin)
             enddo
          case ("nonsu2")
             do i=1,min(flen,Nbath)
                read(unit,*)(&
                     dmft_bath_%e(ispin,1,i),&
                     (&
                     dmft_bath_%v(ispin,iorb,i),&
                     dmft_bath_%u(ispin,iorb,i),&
                     iorb=1,Norb),&
                     ispin=1,Nspin)
             enddo
          end select
          !
       case ('replica')
          read(unit,*)
          !
          read(unit,*)dmft_bath%Nbasis
          do i=1,Nbath
             read(unit,*)dmft_bath_%item(i)%v,&
                  (dmft_bath_%item(i)%lambda(io),io=1,dmft_bath_%Nbasis)
          enddo
          !
          !
       case ('general')
          read(unit,*)
          !
          read(unit,*)dmft_bath%Nbasis
          do i=1,Nbath
             read(unit,*)dmft_bath_%item(i)%vg(:),&
                  (dmft_bath_%item(i)%lambda(io),io=1,dmft_bath_%Nbasis)
          enddo
          !
       end select
       close(unit)
       !
    endif
  end subroutine init_dmft_bath



  !+-------------------------------------------------------------------+
  !PURPOSE  : write out the bath to a given unit with
  ! the following column formatting:
  ! [(Ek_iorb,Vk_iorb)_iorb=1,Norb]_ispin=1,Nspin
  !+-------------------------------------------------------------------+
  subroutine write_dmft_bath(dmft_bath_,unit)
    type(effective_bath) :: dmft_bath_
    integer,optional     :: unit
    integer              :: unit_
    integer              :: i,Nsym
    integer              :: io,jo,iorb,ispin,isym
    real(8)              :: hybr_aux
    complex(8)           :: ho(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb)
    character(len=64)    :: string_fmt,string_fmt_first
    !
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG write_dmft_bath"
#endif
    unit_=LOGfile;if(present(unit))unit_=unit
    if(.not.dmft_bath_%status)stop "write_dmft_bath error: bath not allocated"
    select case(bath_type)
    case default
       !
       select case(ed_mode)
       case default
          write(unit_,"(90(A21,1X))")&
               ((&
               "#Ek_l"//reg(str(iorb))//"_s"//reg(str(ispin)),&
               "Vk_l"//reg(str(iorb))//"_s"//reg(str(ispin)),&
               iorb=1,Norb),ispin=1,Nspin)
          do i=1,Nbath
             write(unit_,"(90(ES21.12,1X))")((&
                  dmft_bath_%e(ispin,iorb,i),&
                  dmft_bath_%v(ispin,iorb,i),&
                  iorb=1,Norb),ispin=1,Nspin)
          enddo
       case ("superc")
          write(unit_,"(90(A21,1X))")&
               ((&
               "#Ek_l"//reg(str(iorb))//"_s"//reg(str(ispin)),&
               "Dk_l"//reg(str(iorb))//"_s"//reg(str(ispin)) ,&
               "Vk_l"//reg(str(iorb))//"_s"//reg(str(ispin)),&
               iorb=1,Norb),ispin=1,Nspin)
          do i=1,Nbath
             write(unit_,"(90(ES21.12,1X))")((&
                  dmft_bath_%e(ispin,iorb,i),&
                  dmft_bath_%d(ispin,iorb,i),&
                  dmft_bath_%v(ispin,iorb,i),&
                  iorb=1,Norb),ispin=1,Nspin)
          enddo
       case ("nonsu2")
          write(unit_,"(90(A21,1X))")&
               ((&
               "#Ek_l"//reg(str(iorb))//"_s"//reg(str(ispin)),&
               "Vak_l"//reg(str(iorb))//"_s"//reg(str(ispin)),&
               "Vbk_l"//reg(str(iorb))//"_s"//reg(str(ispin)),&
               iorb=1,Norb), ispin=1,Nspin)
          do i=1,Nbath
             write(unit,"(90(ES21.12,1X))")((&
                  dmft_bath_%e(ispin,iorb,i),&
                  dmft_bath_%v(ispin,iorb,i),&
                  dmft_bath_%u(ispin,iorb,i),&
                  iorb=1,Norb),ispin=1,Nspin)
          enddo
       end select
       !
    case('hybrid')
       !
       select case(ed_mode)
       case default
          write(unit_,"(90(A21,1X))")(&
               "#Ek_s"//reg(str(ispin)),&
               ("Vk_l"//reg(str(iorb))//"_s"//reg(str(ispin)),iorb=1,Norb),&
               ispin=1,Nspin)
          do i=1,Nbath
             write(unit_,"(90(ES21.12,1X))")(&
                  dmft_bath_%e(ispin,1,i),&
                  (dmft_bath_%v(ispin,iorb,i),iorb=1,Norb),&
                  ispin=1,Nspin)
          enddo
       case ("superc")
          write(unit_,"(90(A21,1X))")(&
               "#Ek_s"//reg(str(ispin)),&
               "Dk_s"//reg(str(ispin)) ,&
               ("Vk_l"//reg(str(iorb))//"_s"//reg(str(ispin)),iorb=1,Norb),&
               ispin=1,Nspin)
          do i=1,Nbath
             write(unit_,"(90(ES21.12,1X))")(&
                  dmft_bath_%e(ispin,1,i),&
                  dmft_bath_%d(ispin,1,i),&
                  (dmft_bath_%v(ispin,iorb,i),iorb=1,Norb),&
                  ispin=1,Nspin)
          enddo
       case ("nonsu2")
          write(unit_,"(90(A21,1X))")(&
               "#Ek_s"//reg(str(ispin)),&
               (&
               "Vak_l"//reg(str(iorb))//"_s"//reg(str(ispin)),&
               "Vbk_l"//reg(str(iorb))//"_s"//reg(str(ispin)),&
               iorb=1,Norb),&
               ispin=1,Nspin)
          do i=1,Nbath
             write(unit_,"(90(ES21.12,1X))")(&
                  dmft_bath_%e(ispin,1,i),    &
                  (dmft_bath_%v(ispin,iorb,i),dmft_bath_%u(ispin,iorb,i),iorb=1,Norb),&
                  ispin=1,Nspin)
          enddo
       end select
       !
    case ('replica')
       !
       string_fmt      ="("//str(Nnambu*Nspin*Norb)//"(A1,F5.2,A1,F5.2,A1,2x))"
       !
       write(unit_,"(90(A21,1X))")"#V_i",("Lambda_i"//reg(str(io)),io=1,dmft_bath_%Nbasis)
       write(unit_,"(I3)")dmft_bath_%Nbasis
       do i=1,Nbath
          write(unit_,"(90(ES21.12,1X))")dmft_bath_%item(i)%v,&
               (dmft_bath_%item(i)%lambda(io),io=1,dmft_bath_%Nbasis)
       enddo
       !
       if(unit_/=LOGfile)then
          write(unit_,*)""
          do isym=1,size(Hreplica_basis)
             Ho = nn2so_reshape(Hreplica_basis(isym)%O,Nnambu*Nspin,Norb)
             do io=1,Nnambu*Nspin*Norb
                write(unit_,string_fmt)&
                     ('(',dreal(Ho(io,jo)),',',dimag(Ho(io,jo)),')',jo =1,Nnambu*Nspin*Norb)
             enddo
             write(unit_,*)""
          enddo
       endif
    case ('general')
       !
       string_fmt      ="("//str(Nnambu*Nspin*Norb)//"(A1,F5.2,A1,F5.2,A1,2x))"
       !
       write(unit_,"(A1,90(A21,1X))")"#",("V_i"//reg(str(io)),io=1,Nspin*Norb),("Lambda_i"//reg(str(io)),io=1,dmft_bath_%Nbasis)
       write(unit_,"(I3)")dmft_bath_%Nbasis
       do i=1,Nbath
          write(unit_,"(90(ES21.12,1X))")(dmft_bath_%item(i)%vg(io),io=1,Nspin*Norb),&
               (dmft_bath_%item(i)%lambda(io),io=1,dmft_bath_%Nbasis) 
       enddo
       !
       if(unit_/=LOGfile)then
          write(unit_,*)""
          do isym=1,size(Hgeneral_basis)
             Ho = nn2so_reshape(Hgeneral_basis(isym)%O,Nnambu*Nspin,Norb)
             do io=1,Nnambu*Nspin*Norb
                write(unit_,string_fmt)&
                     ('(',dreal(Ho(io,jo)),',',dimag(Ho(io,jo)),')',jo =1,Nnambu*Nspin*Norb)
             enddo
             write(unit_,*)""
          enddo
       endif

    end select
  end subroutine write_dmft_bath






  !+-------------------------------------------------------------------+
  !PURPOSE  : save the bath to a given file using the write bath
  ! procedure and formatting:
  !+-------------------------------------------------------------------+
  subroutine save_dmft_bath(dmft_bath_,file,used)
    type(effective_bath)      :: dmft_bath_
    character(len=*),optional :: file
    character(len=256)        :: file_
    logical,optional          :: used
    logical                   :: used_
    character(len=16)         :: extension
    integer                   :: unit_
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG save_dmft_bath"
#endif
    if(.not.dmft_bath_%status)stop "save_dmft_bath error: bath is not allocated"
    used_=.false.;if(present(used))used_=used
    extension=".restart";if(used_)extension=".used"
    file_=str(str(Hfile)//str(ed_file_suffix)//str(extension))
    if(present(file))file_=str(file)    
    unit_=free_unit()
    if(MpiMaster)then
       open(unit_,file=str(file_))
       call write_dmft_bath(dmft_bath_,unit_)
       close(unit_)
    endif
  end subroutine save_dmft_bath




  !+-------------------------------------------------------------------+
  !PURPOSE  : set the bath components from a given user provided
  ! bath-array
  !+-------------------------------------------------------------------+
  subroutine set_dmft_bath(bath_,dmft_bath_)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    integer                :: stride,io,jo,i
    integer                :: iorb,ispin,jorb,jspin,ibath
    logical                :: check
    !
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG set_dmft_bath: dmft_bath <- user_bath"
#endif
    if(.not.dmft_bath_%status)stop "set_dmft_bath error: bath not allocated"
    check = check_bath_dimension(bath_)
    if(.not.check)stop "set_dmft_bath error: wrong bath dimensions"
    !
    select case(bath_type)
    case default
       !
       select case(ed_mode)
          !
       case default
          stride = 0
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   dmft_bath_%e(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
          stride = Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   dmft_bath_%v(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
          !
       case ("superc")
          stride = 0
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   dmft_bath_%e(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
          stride = Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   dmft_bath_%d(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
          stride = 2*Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   dmft_bath_%v(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
          !
       case("nonsu2")
          stride = 0
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   dmft_bath_%e(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
          stride = Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   dmft_bath_%v(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
          stride = 2*Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   dmft_bath_%u(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
       end select
       !
       !
    case ('hybrid')
       !
       select case(ed_mode)
       case default
          stride = 0
          do ispin=1,Nspin
             do i=1,Nbath
                io = stride + i + (ispin-1)*Nbath
                dmft_bath_%e(ispin,1,i) = bath_(io)
             enddo
          enddo
          stride = Nspin*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                   dmft_bath_%v(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
          !
       case ("superc")
          stride = 0
          do ispin=1,Nspin
             do i=1,Nbath
                io = stride + i + (ispin-1)*Nbath
                dmft_bath_%e(ispin,1,i) = bath_(io)
             enddo
          enddo
          stride = Nspin*Nbath
          do ispin=1,Nspin
             do i=1,Nbath
                io = stride + i + (ispin-1)*Nbath
                dmft_bath_%d(ispin,1,i) = bath_(io)
             enddo
          enddo
          stride = 2*Nspin*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                   dmft_bath_%v(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
          !
       case("nonsu2")
          stride = 0
          do ispin=1,Nspin
             do i=1,Nbath
                io = stride + i + (ispin-1)*Nbath
                dmft_bath_%e(ispin,1,i) = bath_(io)
             enddo
          enddo
          stride = Nspin*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                   dmft_bath_%v(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
          stride = Nspin*Nbath + Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                   dmft_bath_%u(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
       end select
       !
       !
    case ('replica')
       !
       stride = 1
       !Get Nbasis
       dmft_bath_%Nbasis = NINT(bath_(stride))
       !get Lambdas
       do ibath=1,Nbath
          stride = stride + 1
          dmft_bath_%item(ibath)%v = bath_(stride)
          dmft_bath_%item(ibath)%lambda=bath_(stride+1 :stride+dmft_bath_%Nbasis)
          stride=stride+dmft_bath_%Nbasis
       enddo
    case ('general')
       !
       stride = 1
       !Get Nbasis
       dmft_bath_%Nbasis = NINT(bath_(stride))
       !get Lambdas
       do ibath=1,Nbath
          dmft_bath_%item(ibath)%vg(:) = bath_(stride+1:stride+Nspin*Norb)
          stride = stride + Nspin*Norb
          dmft_bath_%item(ibath)%lambda=bath_(stride+1 :stride+dmft_bath_%Nbasis)
          stride=stride+dmft_bath_%Nbasis
       enddo
    end select
  end subroutine set_dmft_bath



  !+-------------------------------------------------------------------+
  !PURPOSE  : copy the bath components back to a 1-dim array
  !+-------------------------------------------------------------------+
  subroutine get_dmft_bath(dmft_bath_,bath_)
    type(effective_bath)   :: dmft_bath_
    real(8),dimension(:)   :: bath_
    real(8)                :: hrep_aux(Nspin*Norb,Nspin*Norb)
    integer                :: stride,io,jo,i
    integer                :: iorb,ispin,jorb,jspin,ibath,maxspin
    logical                :: check
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG get_dmft_bath: dmft_bath -> user_bath"
#endif
    if(.not.dmft_bath_%status)stop "get_dmft_bath error: bath not allocated"
    check=check_bath_dimension(bath_)
    if(.not.check)stop "get_dmft_bath error: wrong bath dimensions"
    !
    bath_ = 0d0
    !
    select case(bath_type)
    case default
       !
       select case(ed_mode)
       case default
          stride = 0
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   bath_(io) = dmft_bath_%e(ispin,iorb,i)
                enddo
             enddo
          enddo
          stride = Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   bath_(io) = dmft_bath_%v(ispin,iorb,i)
                enddo
             enddo
          enddo
          !
       case ("superc")
          stride = 0
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   bath_(io) = dmft_bath_%e(ispin,iorb,i)
                enddo
             enddo
          enddo
          stride = Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   bath_(io) =  dmft_bath_%d(ispin,iorb,i)
                enddo
             enddo
          enddo
          stride = 2*Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   bath_(io) =  dmft_bath_%v(ispin,iorb,i)
                enddo
             enddo
          enddo
          !
       case("nonsu2")
          stride = 0
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   bath_(io) = dmft_bath_%e(ispin,iorb,i)
                enddo
             enddo
          enddo
          stride = Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   bath_(io) = dmft_bath_%v(ispin,iorb,i)
                enddo
             enddo
          enddo
          stride = 2*Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   bath_(io) = dmft_bath_%u(ispin,iorb,i)
                enddo
             enddo
          enddo
       end select
       !
    case ('hybrid')
       !
       select case(ed_mode)
       case default
          stride = 0
          do ispin=1,Nspin
             do i=1,Nbath
                io = stride + i + (ispin-1)*Nbath
                bath_(io) =  dmft_bath_%e(ispin,1,i)
             enddo
          enddo
          stride = Nspin*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                   bath_(io) =  dmft_bath_%v(ispin,iorb,i)
                enddo
             enddo
          enddo
          !
       case ("superc")
          stride = 0
          do ispin=1,Nspin
             do i=1,Nbath
                io = stride + i + (ispin-1)*Nbath
                bath_(io) =  dmft_bath_%e(ispin,1,i)
             enddo
          enddo
          stride = Nspin*Nbath
          do ispin=1,Nspin
             do i=1,Nbath
                io = stride + i + (ispin-1)*Nbath
                bath_(io) =  dmft_bath_%d(ispin,1,i)
             enddo
          enddo
          stride = 2*Nspin*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                   bath_(io) =  dmft_bath_%v(ispin,iorb,i)
                enddo
             enddo
          enddo
          !
       case("nonsu2")
          stride = 0
          do ispin=1,Nspin
             do i=1,Nbath
                io = stride + i + (ispin-1)*Nbath
                bath_(io) = dmft_bath_%e(ispin,1,i)
             enddo
          enddo
          stride = Nspin*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                   bath_(io) = dmft_bath_%v(ispin,iorb,i)
                enddo
             enddo
          enddo
          stride = Nspin*Nbath + Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                   bath_(io) = dmft_bath_%u(ispin,iorb,i)
                enddo
             enddo
          enddo
       end select
       !
       !
    case ('replica')
       !
       stride = 1
       bath_(stride)=dmft_bath_%Nbasis
       do ibath=1,Nbath
          stride = stride + 1
          bath_(stride)=dmft_bath_%item(ibath)%v
          bath_(stride+1 : stride+dmft_bath_%Nbasis)=dmft_bath_%item(ibath)%lambda
          stride=stride+dmft_bath_%Nbasis
       enddo
    case ('general')
       !
       stride = 1
       bath_(stride)=dmft_bath_%Nbasis
       do ibath=1,Nbath
          bath_(stride+1:stride+Nspin*Norb)=dmft_bath_%item(ibath)%vg(:)
          stride = stride + Nspin*Norb
          bath_(stride+1 : stride+dmft_bath_%Nbasis)=dmft_bath_%item(ibath)%lambda
          stride=stride+dmft_bath_%Nbasis
       enddo
    end select    
  end subroutine get_dmft_bath





  !##################################################################
  !
  !     W_hyb PROCEDURES
  !
  !##################################################################
  !+-----------------------------------------------------------------------------+!
  !PURPOSE: build up the all-spin hybridization matrix W_{ss`}
  !+-----------------------------------------------------------------------------+!
  function get_Whyb_matrix_1orb(v,u) result(w)
    real(8),dimension(Nspin,Nbath)       :: v,u
    real(8),dimension(Nspin,Nspin,Nbath) :: w
    integer                              :: ispin
    !
    ! if(ed_para)then
    !    do ispin=1,Nspin
    !       w(ispin,ispin,:) = v(1,:)
    !    enddo
    !    w(1,Nspin,:) = u(1,:)
    !    w(Nspin,1,:) = u(1,:)
    ! else
    do ispin=1,Nspin
       w(ispin,ispin,:) = v(ispin,:)
    enddo
    w(1,Nspin,:) = u(1,:)
    w(Nspin,1,:) = u(2,:)
    ! endif
    !
  end function get_Whyb_matrix_1orb

  function get_Whyb_matrix_Aorb(v,u) result(w)
    real(8),dimension(Nspin,Norb,Nbath)       :: v,u
    real(8),dimension(Nspin,Nspin,Norb,Nbath) :: w
    integer                                   :: ispin
    !
    ! if(ed_para)then
    !    do ispin=1,Nspin
    !       w(ispin,ispin,:,:) = v(1,:,:)
    !    enddo
    !    w(1,Nspin,:,:) = u(1,:,:)
    !    w(Nspin,1,:,:) = u(1,:,:)
    ! else
    do ispin=1,Nspin
       w(ispin,ispin,:,:) = v(ispin,:,:)
    enddo
    w(1,Nspin,:,:) = u(1,:,:)
    w(Nspin,1,:,:) = u(2,:,:)
    ! endif
    !
  end function get_Whyb_matrix_Aorb

  function get_Whyb_matrix_dmft_bath(dmft_bath_) result(w)
    type(effective_bath)                      :: dmft_bath_
    real(8),dimension(Nspin,Nspin,Norb,Nbath) :: w
    integer                                   :: ispin
    !
    ! if(ed_para)then
    !    do ispin=1,Nspin
    !       w(ispin,ispin,:,:) = dmft_bath_%v(1,:,:)
    !    enddo
    !    w(1,Nspin,:,:) = dmft_bath_%u(1,:,:)
    !    w(Nspin,1,:,:) = dmft_bath_%u(1,:,:)
    ! else
    do ispin=1,Nspin
       w(ispin,ispin,:,:) = dmft_bath_%v(ispin,:,:)
    enddo
    w(1,Nspin,:,:) = dmft_bath_%u(1,:,:)
    w(Nspin,1,:,:) = dmft_bath_%u(2,:,:)
    ! endif
    !
  end function get_Whyb_matrix_dmft_bath




























  !##################################################################
  !
  !     H_REPLICA ROUTINES:
  !
  !##################################################################
  !-------------------------------------------------------------------!
  ! PURPOSE: INITIALIZE INTERNAL Hreplica STRUCTURES
  !-------------------------------------------------------------------!

  !allocate GLOBAL basis for H (used for Hbath) and vectors coefficient
  subroutine allocate_hreplica(Nsym)
    integer          :: Nsym
    integer          :: isym
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG allocate_Hreplica"
#endif
    if(allocated(Hreplica_basis))deallocate(Hreplica_basis)
    if(allocated(Hreplica_lambda))deallocate(Hreplica_lambda)
    !
    allocate(Hreplica_basis(Nsym))
    allocate(Hreplica_lambda(Nbath,Nsym))
    do isym=1,Nsym
       allocate(Hreplica_basis(isym)%O(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb))
       Hreplica_basis(isym)%O=zero
       Hreplica_lambda(:,isym)=0d0
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


  subroutine init_Hreplica_symmetries_site(Hvec,lambdavec)
    complex(8),dimension(:,:,:,:,:) :: Hvec      ![Nnambu*Nspin,Nnambu*Nspin,Norb,Norb,Nsym]
    real(8),dimension(:,:)          :: lambdavec ![Nbath,Nsym]
    integer                         :: isym,Nsym
    logical                         :: bool
    !
    if(ed_mode=="superc")Nnambu=2
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG init_Hreplica_symmetries_site: from {[Hs,Lam]}_b"
#endif
    !
    if(size(lambdavec,1)/=Nbath)then
       write(*,*) "                                                                               "
       write(*,*) "ERROR: if you are trying to init Hreplica for inequivalent sites please note   "
       write(*,*) "       that the lambdasym array /MUST/ have [Nineq]x[Nbath]x[Nsym] shape.      "
       write(*,*) "       The legacy [Nineq]x[Nsym] is not supported anymore, for it would shadow "
       write(*,*) "       the new recommended [Nbath]x[Nsym] shape for the single impurity case.  "
       write(*,*) "                                                                               "
       stop ! This unfortunately still leaves room for nasty problems if Nbath==Nineq, but that's it...
    else
       Nsym=size(lambdavec,2)
    endif
    !
    call assert_shape(Hvec,[Nnambu*Nspin,Nnambu*Nspin,Norb,Norb,Nsym],"init_Hreplica_symmetries","Hvec")
    !
    !CHECK NAMBU and HERMITICTY of each Hvec
    do isym=1,Nsym
       select case(ed_mode)
       case default
          bool = check_herm(nn2so_reshape(Hvec(:,:,:,:,isym),Nspin,Norb),Nspin*Norb)
       case("superc")
          bool = check_nambu(nn2so_reshape(Hvec(:,:,:,:,isym),Nnambu*Nspin,Norb),Nspin*Norb)
       end select
       if(.not.bool)then
          write(LOGfile,"(A)")"init_Hreplica_symmetries_site ERROR: not Hermitian/Nambu of replica basis O_"//str(isym)
          stop
       endif
    enddo
    !
    call allocate_hreplica(Nsym)
    !
    !
    do isym=1,Nsym
       Hreplica_lambda(:,isym)  = lambdavec(:,isym)
       Hreplica_basis(isym)%O = Hvec(:,:,:,:,isym)
    enddo
    !
    if(ed_verbose>2)then
       do ibath=1,Nbath
          write(*,*) "Hreplica #"//str(ibath)//":"
          call print_hloc(Hreplica_build(Hreplica_lambda(ibath,:)))
       enddo
    endif
    !
  end subroutine init_Hreplica_symmetries_site

  subroutine init_Hreplica_symmetries_legacy(Hvec,lambdavec)
    complex(8),dimension(:,:,:,:,:) :: Hvec      ![size(H),Nsym]
    real(8),dimension(:)            :: lambdavec ![Nsym]
    integer                         :: isym,Nsym
    logical                         :: bool
    !
    if(ed_mode=="superc")Nnambu=2
    Nsym=size(lambdavec)
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG init_Hreplica_symmetries_legacy: from {[Hs,Lam]}_b"
#endif
    call assert_shape(Hvec,[Nnambu*Nspin,Nnambu*Nspin,Norb,Norb,Nsym],"init_Hreplica_symmetries","Hvec")
    !
    !CHECK NAMBU and HERMITICTY of each Hvec
    do isym=1,Nsym
       select case(ed_mode)
       case default
          bool = check_herm(nn2so_reshape(Hvec(:,:,:,:,isym),Nspin,Norb),Nspin*Norb)
       case("superc")
          bool = check_nambu(nn2so_reshape(Hvec(:,:,:,:,isym),Nnambu*Nspin,Norb),Nspin*Norb)
       end select
       if(.not.bool)then
          write(LOGfile,"(A)")"init_Hreplica_symmetries_site ERROR: not Hermitian/Nambu of replica basis O_"//str(isym)
          stop
       endif
    enddo
    !
    call allocate_hreplica(Nsym)
    !
    do isym=1,Nsym
       do ibath=1,Nbath
          !> BACK-COMPATIBILITY PATCH (cfr init_dmft_bath)
          Hreplica_lambda(ibath,isym) = lambdavec(isym)
       enddo
       Hreplica_basis(isym)%O = Hvec(:,:,:,:,isym)
    enddo
    !
    ! PRINT DEPRECATION MESSAGE TO LOG
    write(*,*) "                                                                               "
    write(*,*) "WARNING: Passing a single lambdasym vector to ed_set_Hreplica is /deprecated/. "
    write(*,*) "         You should instead define a different lambda for each bath component, "
    write(*,*) "         namely passing a [Nbath]x[Nsym] array instead of a [Nsym] vector.     "
    write(*,*) "         Your single lambda vector has been internally copied into the required"
    write(*,*) "         higher-rank array, so giving each replica the same set of lambdas.    "
    write(*,*) "         >>> This back-compatibility patch might be removed in a future update."
    write(*,*) "                                                                               "
    !
    if(ed_verbose>2)then
       do ibath=1,Nbath
          write(*,*) "Hreplica #"//str(ibath)//":"
          call print_hloc(Hreplica_build(Hreplica_lambda(ibath,:)))
       enddo
    endif
    !
  end subroutine init_Hreplica_symmetries_legacy

  subroutine init_Hreplica_symmetries_lattice(Hvec,lambdavec)
    complex(8),dimension(:,:,:,:,:) :: Hvec      ![size(H),Nsym]
    real(8),dimension(:,:,:)        :: lambdavec ![Nlat,Nbath,Nsym]
    integer                         :: ilat,Nlat
    integer                         :: isym,Nsym
    logical                         :: bool
    !
    if(ed_mode=="superc")Nnambu=2
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG init_Hreplica_symmetries_lattice: from ({[Hs,Lam]}_b)_site"
#endif
    !
    Nlat=size(lambdavec,1)
    Nsym=size(lambdavec,3)
    call assert_shape(Hvec,[Nnambu*Nspin,Nnambu*Nspin,Norb,Norb,Nsym],"init_Hreplica_symmetries","Hvec")
    !
    !CHECK NAMBU and HERMITICTY of each Hvec
    do isym=1,Nsym
       select case(ed_mode)
       case default
          bool = check_herm(nn2so_reshape(Hvec(:,:,:,:,isym),Nspin,Norb),Nspin*Norb)
       case("superc")
          bool = check_nambu(nn2so_reshape(Hvec(:,:,:,:,isym),Nnambu*Nspin,Norb),Nspin*Norb)
       end select
       if(.not.bool)then
          write(LOGfile,"(A)")"init_Hreplica_symmetries_site ERROR: not Hermitian/Nambu of replica basis O_"//str(isym)
          stop
       endif
    enddo
    !
    if(allocated(Hreplica_lambda_ineq))deallocate(Hreplica_lambda_ineq)
    allocate(Hreplica_lambda_ineq(Nlat,Nbath,Nsym))
    call allocate_hreplica(Nsym)
    !
    do isym=1,Nsym
       Hreplica_lambda_ineq(:,:,isym)  = lambdavec(:,:,isym)
       Hreplica_basis(isym)%O = Hvec(:,:,:,:,isym)
    enddo
    !
    if(ed_verbose>2)then
       do ilat=1,Nlat
          write(*,*) "Inequivalent #"//str(ilat)//":"
          do ibath=1,Nbath
             write(*,*) "> Hreplica #"//str(ibath)//":"
             call print_hloc(Hreplica_build(Hreplica_lambda_ineq(ilat,ibath,:)))
          enddo
       enddo
    endif
    !
  end subroutine init_Hreplica_symmetries_lattice

  
  !##################################################################
  !
  !     AUX FUNCTIONS REPLICA:
  !
  !##################################################################
  subroutine Hreplica_site(site)
    integer :: site
    if(site<1.OR.site>size(Hreplica_lambda_ineq,1))stop "ERROR Hreplica_site: site not in [1,Nlat]"
    if(.not.allocated(Hreplica_lambda_ineq))stop "ERROR Hreplica_site: Hreplica_lambda_ineq not allocated"
    Hreplica_lambda(:,:)  = Hreplica_lambda_ineq(site,:,:)
  end subroutine Hreplica_site


  !reconstruct [Nspin,Nspin,Norb,Norb] hamiltonian from basis expansion given [lambda]
  function Hreplica_build(lambdavec) result(H)
    real(8),dimension(:)                                      :: lambdavec ![Nsym]
    integer                                                   :: isym
    complex(8),dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb) :: H
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
    logical,optional                                          :: wdiag,uplo
    logical                                                   :: wdiag_,uplo_
    complex(8),dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb) :: H
    logical,dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb)    :: Hmask
    integer                                                   :: iorb,jorb,ispin,jspin,io,jo
    !
    wdiag_=.false.;if(present(wdiag))wdiag_=wdiag
    uplo_ =.false.;if(present(uplo))  uplo_=uplo
    !
    H = Hreplica_build(Hreplica_lambda(Nbath,:)) !The mask should be replica-independent
    Hmask=.false.
    where(abs(H)>1d-6)Hmask=.true.
    !
    !
    if(wdiag_)then
       do ispin=1,Nnambu*Nspin
          do iorb=1,Norb
             Hmask(ispin,ispin,iorb,iorb)=.true.
          enddo
       enddo
    endif
    !
    if(uplo_)then
       do ispin=1,Nnambu*Nspin
          do jspin=1,Nnambu*Nspin
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



  
  !##################################################################
  !
  !     H_GENERAL ROUTINES:
  !
  !##################################################################
  !-------------------------------------------------------------------!
  ! PURPOSE: INITIALIZE INTERNAL Hgeneral STRUCTURES
  !-------------------------------------------------------------------!

  !allocate GLOBAL basis for H (used for Hbath) and vectors coefficient
  subroutine allocate_hgeneral(Nsym)
    integer          :: Nsym
    integer          :: isym
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG allocate_Hgeneral"
#endif
    if(allocated(Hgeneral_basis))deallocate(Hgeneral_basis)
    if(allocated(Hgeneral_lambda))deallocate(Hgeneral_lambda)
    !
    allocate(Hgeneral_basis(Nsym))
    allocate(Hgeneral_lambda(Nbath,Nsym))
    do isym=1,Nsym
       allocate(Hgeneral_basis(isym)%O(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb))
       Hgeneral_basis(isym)%O=zero
       Hgeneral_lambda(:,isym)=0d0
    enddo
    Hgeneral_status=.true.
  end subroutine allocate_hgeneral


  !deallocate GLOBAL basis for H (used for impHloc and bath) and vectors coefficient
  subroutine deallocate_hgeneral()
    integer              :: isym
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG deallocate_Hgeneral"
#endif
    do isym=1,size(Hgeneral_basis)
       deallocate(Hgeneral_basis(isym)%O)
    enddo
    deallocate(Hgeneral_basis)
    deallocate(Hgeneral_lambda)
    Hgeneral_status=.false.
  end subroutine deallocate_hgeneral


  subroutine init_Hgeneral_symmetries_site(Hvec,lambdavec)
    complex(8),dimension(:,:,:,:,:) :: Hvec      ![Nnambu*Nspin,Nnambu*Nspin,Norb,Norb,Nsym]
    real(8),dimension(:,:)          :: lambdavec ![Nbath,Nsym]
    integer                         :: isym,Nsym
    logical                         :: bool
    !
    if(ed_mode=="superc")Nnambu=2
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG init_Hgeneral_symmetries_site: from {[Hs,Lam]}_b"
#endif
    !
    if(size(lambdavec,1)/=Nbath)then
       write(*,*) "                                                                               "
       write(*,*) "ERROR: if you are trying to init Hgeneral for inequivalent sites please note   "
       write(*,*) "       that the lambdasym array /MUST/ have [Nineq]x[Nbath]x[Nsym] shape.      "
       write(*,*) "       The legacy [Nineq]x[Nsym] is not supported anymore, for it would shadow "
       write(*,*) "       the new recommended [Nbath]x[Nsym] shape for the single impurity case.  "
       write(*,*) "                                                                               "
       stop ! This unfortunately still leaves room for nasty problems if Nbath==Nineq, but that's it...
    else
       Nsym=size(lambdavec,2)
    endif
    !
    call assert_shape(Hvec,[Nnambu*Nspin,Nnambu*Nspin,Norb,Norb,Nsym],"init_Hgeneral_symmetries","Hvec")
    !
    !CHECK NAMBU and HERMITICTY of each Hvec
    do isym=1,Nsym
       select case(ed_mode)
       case default
          bool = check_herm(nn2so_reshape(Hvec(:,:,:,:,isym),Nspin,Norb),Nspin*Norb)
       case("superc")
          bool = check_nambu(nn2so_reshape(Hvec(:,:,:,:,isym),Nnambu*Nspin,Norb),Nspin*Norb)
       end select
       if(.not.bool)then
          write(LOGfile,"(A)")"init_Hgeneral_symmetries_site ERROR: not Hermitian/Nambu of general basis O_"//str(isym)
          stop
       endif
    enddo
    !
    call allocate_hgeneral(Nsym)
    !
    !
    do isym=1,Nsym
       Hgeneral_lambda(:,isym)  = lambdavec(:,isym)
       Hgeneral_basis(isym)%O = Hvec(:,:,:,:,isym)
    enddo
    !
    if(ed_verbose>2)then
       do ibath=1,Nbath
          write(*,*) "Hgeneral #"//str(ibath)//":"
          call print_hloc(Hgeneral_build(Hgeneral_lambda(ibath,:)))
       enddo
    endif
    !
  end subroutine init_Hgeneral_symmetries_site

  subroutine init_Hgeneral_symmetries_legacy(Hvec,lambdavec)
    complex(8),dimension(:,:,:,:,:) :: Hvec      ![size(H),Nsym]
    real(8),dimension(:)            :: lambdavec ![Nsym]
    integer                         :: isym,Nsym
    logical                         :: bool
    !
    if(ed_mode=="superc")Nnambu=2
    Nsym=size(lambdavec)
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG init_Hgeneral_symmetries_legacy: from {[Hs,Lam]}_b"
#endif
    call assert_shape(Hvec,[Nnambu*Nspin,Nnambu*Nspin,Norb,Norb,Nsym],"init_Hgeneral_symmetries","Hvec")
    !
    !CHECK NAMBU and HERMITICTY of each Hvec
    do isym=1,Nsym
       select case(ed_mode)
       case default
          bool = check_herm(nn2so_reshape(Hvec(:,:,:,:,isym),Nspin,Norb),Nspin*Norb)
       case("superc")
          bool = check_nambu(nn2so_reshape(Hvec(:,:,:,:,isym),Nnambu*Nspin,Norb),Nspin*Norb)
       end select
       if(.not.bool)then
          write(LOGfile,"(A)")"init_Hgeneral_symmetries_site ERROR: not Hermitian/Nambu of general basis O_"//str(isym)
          stop
       endif
    enddo
    !
    call allocate_hgeneral(Nsym)
    !
    do isym=1,Nsym
       do ibath=1,Nbath
          !> BACK-COMPATIBILITY PATCH (cfr init_dmft_bath)
          Hgeneral_lambda(ibath,isym) = lambdavec(isym)
       enddo
       Hgeneral_basis(isym)%O = Hvec(:,:,:,:,isym)
    enddo
    !
    ! PRINT DEPRECATION MESSAGE TO LOG
    write(*,*) "                                                                               "
    write(*,*) "WARNING: Passing a single lambdasym vector to ed_set_Hgeneral is /deprecated/. "
    write(*,*) "         You should instead define a different lambda for each bath component, "
    write(*,*) "         namely passing a [Nbath]x[Nsym] array instead of a [Nsym] vector.     "
    write(*,*) "         Your single lambda vector has been internally copied into the required"
    write(*,*) "         higher-rank array, so giving each general the same set of lambdas.    "
    write(*,*) "         >>> This back-compatibility patch might be removed in a future update."
    write(*,*) "                                                                               "
    !
    if(ed_verbose>2)then
       do ibath=1,Nbath
          write(*,*) "Hgeneral #"//str(ibath)//":"
          call print_hloc(Hgeneral_build(Hgeneral_lambda(ibath,:)))
       enddo
    endif
    !
  end subroutine init_Hgeneral_symmetries_legacy

  subroutine init_Hgeneral_symmetries_lattice(Hvec,lambdavec)
    complex(8),dimension(:,:,:,:,:) :: Hvec      ![size(H),Nsym]
    real(8),dimension(:,:,:)        :: lambdavec ![Nlat,Nbath,Nsym]
    integer                         :: ilat,Nlat
    integer                         :: isym,Nsym
    logical                         :: bool
    !
    if(ed_mode=="superc")Nnambu=2
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG init_Hgeneral_symmetries_lattice: from ({[Hs,Lam]}_b)_site"
#endif
    !
    Nlat=size(lambdavec,1)
    Nsym=size(lambdavec,3)
    call assert_shape(Hvec,[Nnambu*Nspin,Nnambu*Nspin,Norb,Norb,Nsym],"init_Hgeneral_symmetries","Hvec")
    !
    !CHECK NAMBU and HERMITICTY of each Hvec
    do isym=1,Nsym
       select case(ed_mode)
       case default
          bool = check_herm(nn2so_reshape(Hvec(:,:,:,:,isym),Nspin,Norb),Nspin*Norb)
       case("superc")
          bool = check_nambu(nn2so_reshape(Hvec(:,:,:,:,isym),Nnambu*Nspin,Norb),Nspin*Norb)
       end select
       if(.not.bool)then
          write(LOGfile,"(A)")"init_Hgeneral_symmetries_site ERROR: not Hermitian/Nambu of general basis O_"//str(isym)
          stop
       endif
    enddo
    !
    if(allocated(Hgeneral_lambda_ineq))deallocate(Hgeneral_lambda_ineq)
    allocate(Hgeneral_lambda_ineq(Nlat,Nbath,Nsym))
    call allocate_hgeneral(Nsym)
    !
    do isym=1,Nsym
       Hgeneral_lambda_ineq(:,:,isym)  = lambdavec(:,:,isym)
       Hgeneral_basis(isym)%O = Hvec(:,:,:,:,isym)
    enddo
    !
    if(ed_verbose>2)then
       do ilat=1,Nlat
          write(*,*) "Inequivalent #"//str(ilat)//":"
          do ibath=1,Nbath
             write(*,*) "> Hgeneral #"//str(ibath)//":"
             call print_hloc(Hgeneral_build(Hgeneral_lambda_ineq(ilat,ibath,:)))
          enddo
       enddo
    endif
    !
  end subroutine init_Hgeneral_symmetries_lattice

  
  !##################################################################
  !
  !     AUX FUNCTIONS GENERAL:
  !
  !##################################################################
  subroutine Hgeneral_site(site)
    integer :: site
    if(site<1.OR.site>size(Hgeneral_lambda_ineq,1))stop "ERROR Hgeneral_site: site not in [1,Nlat]"
    if(.not.allocated(Hgeneral_lambda_ineq))stop "ERROR Hgeneral_site: Hgeneral_lambda_ineq not allocated"
    Hgeneral_lambda(:,:)  = Hgeneral_lambda_ineq(site,:,:)
  end subroutine Hgeneral_site


  !reconstruct [Nspin,Nspin,Norb,Norb] hamiltonian from basis expansion given [lambda]
  function Hgeneral_build(lambdavec) result(H)
    real(8),dimension(:)                                      :: lambdavec ![Nsym]
    integer                                                   :: isym
    complex(8),dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb) :: H
    !
    if(.not.Hgeneral_status)STOP "ERROR Hgeneral_build: Hgeneral_basis is not setup"
    if(size(lambdavec)/=size(Hgeneral_basis)) STOP "ERROR Hgeneral_build: Wrong coefficient vector size"
    H=zero
    do isym=1,size(lambdavec)
       H=H+lambdavec(isym)*Hgeneral_basis(isym)%O
    enddo
  end function Hgeneral_build



  !+-------------------------------------------------------------------+
  !PURPOSE  : Create bath mask
  !+-------------------------------------------------------------------+
  function Hgeneral_mask(wdiag,uplo) result(Hmask)
    logical,optional                                          :: wdiag,uplo
    logical                                                   :: wdiag_,uplo_
    complex(8),dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb) :: H
    logical,dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb)    :: Hmask
    integer                                                   :: iorb,jorb,ispin,jspin,io,jo
    !
    wdiag_=.false.;if(present(wdiag))wdiag_=wdiag
    uplo_ =.false.;if(present(uplo))  uplo_=uplo
    !
    H = Hgeneral_build(Hgeneral_lambda(Nbath,:)) !The mask should be general-independent
    Hmask=.false.
    where(abs(H)>1d-6)Hmask=.true.
    !
    !
    if(wdiag_)then
       do ispin=1,Nnambu*Nspin
          do iorb=1,Norb
             Hmask(ispin,ispin,iorb,iorb)=.true.
          enddo
       enddo
    endif
    !
    if(uplo_)then
       do ispin=1,Nnambu*Nspin
          do jspin=1,Nnambu*Nspin
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
  end function Hgeneral_mask

  !+-------------------------------------------------------------------+
  !AUXILIARY FUNCTIONS AND SUBROUTINES
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !PURPOSE  : Check if a matrix is the identity
  !+-------------------------------------------------------------------+
  function is_identity_nn(mnnn) result(flag)
    complex(8),dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb) :: mnnn
    real(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb) :: mtmp
    integer                                  :: i,j
    logical                                  :: flag
    !
    flag=.true.
    !
    mtmp=nn2so_reshape(mnnn,Nnambu*Nspin,Norb)
    !
    do i=1,Nnambu*Nspin*Norb-1
       if((mtmp(i,i).ne.mtmp(i+1,i+1)).or.(mtmp(i,i).lt.1.d-6))flag=.false.
    enddo
    !
    do i=1,Nnambu*Nspin*Norb
       do j=1,Nnambu*Nspin*Norb
          if((i.ne.j).and.(mtmp(i,j).gt.1.d-6))flag=.false.
       enddo
    enddo
    !
  end function is_identity_nn

  function is_identity_so(mlso) result(flag)
    complex(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb) :: mlso
    real(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb) :: mtmp
    integer                                  :: i,j
    logical                                  :: flag
    !
    flag=.true.
    !
    mtmp=mlso
    !
    do i=1,Nnambu*Nspin*Norb-1
       if((mtmp(i,i).ne.mtmp(i+1,i+1)).or.(mtmp(i,i).lt.1.d-6))flag=.false.
    enddo
    !
    do i=1,Nnambu*Nspin*Norb
       do j=1,Nnambu*Nspin*Norb
          if((i.ne.j).and.(mtmp(i,j).gt.1.d-6))flag=.false.
       enddo
    enddo
    !
  end function is_identity_so



  !+-------------------------------------------------------------------+
  !PURPOSE  : Check if a matrix is diagonal
  !+-------------------------------------------------------------------+
  function is_diagonal_nn(mnnn) result(flag)
    complex(8),dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb) :: mnnn
    complex(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb) :: mtmp
    integer                                     :: i,j
    logical                                     :: flag
    !
    flag=.true.
    !
    mtmp=abs((nn2so_reshape(mnnn,Nnambu*Nspin,Norb)))
    !
    do i=1,Nnambu*Nspin*Norb
       do j=1,Nnambu*Nspin*Norb
          if((i.ne.j).and.(abs(mtmp(i,j)).gt.1.d-6))flag=.false.
       enddo
    enddo
    !
  end function is_diagonal_nn

  function is_diagonal_so(mlso) result(flag)
    complex(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb) :: mlso
    complex(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb) :: mtmp
    integer                                     :: i,j
    logical                                     :: flag
    !
    flag=.true.
    !
    mtmp=abs((mlso))
    !
    do i=1,Nnambu*Nspin*Norb
       do j=1,Nnambu*Nspin*Norb
          if((i.ne.j).and.(abs(mtmp(i,j)).gt.1.d-6))flag=.false.
       enddo
    enddo
    !
  end function is_diagonal_so






  function check_herm(A,N,error) result(bool)
    integer,intent(in)                   :: N
    complex(8),dimension(N,N),intent(in) :: A
    logical                              :: bool
    real(8),optional                     :: error
    real(8)                              :: error_
    error_ = 1d-6 ; if(present(error))error_=error
    bool   = all(abs(A - conjg(transpose(A)))<error_)
  end function check_herm


  function check_nambu(A,N,error) result(bool)
    integer,intent(in)                       :: N
    complex(8),dimension(2*N,2*N),intent(in) :: A
    complex(8),dimension(N,N)                :: h11,h22
    logical                                  :: bool
    real(8),optional                         :: error
    real(8)                                  :: error_
    error_ = 1d-6 ; if(present(error))error_=error
    h11    = A(1:N    ,1:N)
    h22    = A(N+1:2*N,N+1:2*N)
    bool   = check_herm(A,2*N,error_) !this checks also for F = A_12, s.t. A_21=herm(A_12)
    bool   = bool.AND.( all(abs(h22 + conjg(h11))<error_) )
  end function check_nambu




END MODULE ED_BATH










! public :: get_bath_component_dimension
! public :: get_bath_component
! public :: set_bath_component
! public :: copy_bath_component
! !



! !+-------------------------------------------------------------------+
! !PURPOSE  : Inquire the correct bath size to allocate the
! ! the bath array in the calling program.
! !
! ! Get size of each dimension of the component array.
! ! The Result is an rank 1 integer array Ndim with dimension:
! ! 3 for get_component_size_bath
! ! 2 for get_spin_component_size_bath & get_orb_component_size_bath
! ! 1 for get_spin_orb_component_size_bath
! !+-------------------------------------------------------------------+
! function get_bath_component_dimension(type) result(Ndim)
!   character(len=1) :: type
!   integer          :: Ndim(3)
!   call  check_bath_component(type)
!   select case(bath_type)
!   case default
!      Ndim=[Nspin,Norb,Nbath]
!   case('hybrid')
!      select case(ed_mode)
!      case default
!         select case(type)
!         case('e')
!            Ndim=[Nspin,1,Nbath]
!         case('v')
!            Ndim=[Nspin,Norb,Nbath]
!         end select
!      case ("superc")
!         select case(type)
!         case('e','d')
!            Ndim=[Nspin,1,Nbath]
!         case('v')
!            Ndim=[Nspin,Norb,Nbath]
!         end select
!      case ("nonsu2")
!         select case(type)
!         case('e')
!            Ndim=[Nspin,1,Nbath]
!         case('v','u')
!            Ndim=[Nspin,Norb,Nbath]
!         end select
!      end select
!   end select
! end function get_bath_component_dimension


! !+-----------------------------------------------------------------------------+!
! !PURPOSE: check that the input array hsa the correct dimensions specified
! ! for the choice of itype and possiblty ispin and/or iorb.
! !+-----------------------------------------------------------------------------+!
! subroutine assert_bath_component_size(array,type,string1,string2)
!   real(8),dimension(:,:,:) :: array
!   character(len=1)         :: type
!   character(len=*)         :: string1,string2
!   integer                  :: Ndim(3)
!   Ndim = get_bath_component_dimension(type)
!   call assert_shape(Array,Ndim,reg(string1),reg(string2))
! end subroutine assert_bath_component_size








! !+-----------------------------------------------------------------------------+!
! !PURPOSE: Get a specified itype,ispin,iorb component of the user bath.
! ! The component is returned into an Array of rank D
! ! get_full_component_bath    : return the entire itype component (D=3)
! ! get_spin_component_bath    : return the itype component for the select ispin (D=2)
! ! get_spin_orb_component_bath: return the itype component for the select ispin & iorb (D=1)
! !+-----------------------------------------------------------------------------+!
! subroutine get_bath_component(array,bath_,type)
!   real(8),dimension(:,:,:) :: array
!   real(8),dimension(:)     :: bath_
!   character(len=1)         :: type
!   logical                  :: check
!   type(effective_bath)     :: dmft_bath_
!   !
!   check= check_bath_dimension(bath_)
!   if(.not.check)stop "get_component_bath error: wrong bath dimensions"
!   call allocate_dmft_bath(dmft_bath_)
!   call set_dmft_bath(bath_,dmft_bath_)
!   call assert_bath_component_size(array,type,"get_bath_component","Array")
!   call check_bath_component(type)
!   select case(ed_mode)
!   case default
!      select case(type)
!      case('e')
!         Array = dmft_bath_%e(:,:,:)
!      case('v')
!         Array = dmft_bath_%v(:,:,:)
!      end select
!   case ("superc")
!      select case(type)
!      case('e')
!         Array = dmft_bath_%e(:,:,:)
!      case('d')
!         Array = dmft_bath_%d(:,:,:)
!      case('v')
!         Array = dmft_bath_%v(:,:,:)
!      end select
!   case ("nonsu2")
!      select case(type)
!      case('e')
!         Array = dmft_bath_%e(:,:,:)
!      case('v')
!         Array = dmft_bath_%v(:,:,:)
!      case('u')
!         Array = dmft_bath_%u(:,:,:)
!      end select
!   end select
!   call deallocate_dmft_bath(dmft_bath_)
! end subroutine get_bath_component


! !+-----------------------------------------------------------------------------+!
! !PURPOSE: Set a specified itype,ispin,iorb component of the user bath.
! !+-----------------------------------------------------------------------------+!
! subroutine set_bath_component(array,bath_,type)
!   real(8),dimension(:,:,:) :: array
!   real(8),dimension(:)     :: bath_
!   character(len=1)         :: type
!   logical                  :: check
!   type(effective_bath)     :: dmft_bath_
!   !
!   check= check_bath_dimension(bath_)
!   if(.not.check)stop "set_component_bath error: wrong bath dimensions"
!   call allocate_dmft_bath(dmft_bath_)
!   call set_dmft_bath(bath_,dmft_bath_)
!   call assert_bath_component_size(array,type,"set_bath_component","Array")
!   call check_bath_component(type)
!   select case(ed_mode)
!   case default
!      select case(type)
!      case('e')
!         dmft_bath_%e(:,:,:) = Array
!      case('v')
!         dmft_bath_%v(:,:,:) = Array
!      end select
!   case ("superc")
!      select case(type)
!      case('e')
!         dmft_bath_%e(:,:,:) = Array
!      case('d')
!         dmft_bath_%d(:,:,:) = Array
!      case('v')
!         dmft_bath_%v(:,:,:) = Array
!      end select
!   case ("nonsu2")
!      select case(type)
!      case('e')
!         dmft_bath_%e(:,:,:) = Array
!      case('v')
!         dmft_bath_%v(:,:,:) = Array
!      case('u')
!         dmft_bath_%u(:,:,:) = Array
!      end select
!   end select
!   call get_dmft_bath(dmft_bath_,bath_)
!   call deallocate_dmft_bath(dmft_bath_)
! end subroutine set_bath_component



! !+-----------------------------------------------------------------------------+!
! !PURPOSE: Copy a specified component of IN bath to the OUT bath.
! !+-----------------------------------------------------------------------------+!
! subroutine copy_bath_component(bathIN,bathOUT,type)
!   real(8),dimension(:)     :: bathIN,bathOUT
!   character(len=1)         :: type
!   logical                  :: check
!   type(effective_bath)     :: dIN,dOUT
!   !
!   check= check_bath_dimension(bathIN)
!   if(.not.check)stop "copy_component_bath error: wrong bath dimensions IN"
!   check= check_bath_dimension(bathOUT)
!   if(.not.check)stop "copy_component_bath error: wrong bath dimensions OUT"
!   call allocate_dmft_bath(dIN)
!   call allocate_dmft_bath(dOUT)
!   call set_dmft_bath(bathIN,dIN)
!   call set_dmft_bath(bathOUT,dOUT)
!   call check_bath_component(type)
!   select case(ed_mode)
!   case default
!      select case(type)
!      case('e')
!         dOUT%e(:,:,:)  = dIN%e(:,:,:)
!      case('v')
!         dOUT%v(:,:,:)  = dIN%v(:,:,:)
!      end select
!   case ("superc")
!      select case(type)
!      case('e')
!         dOUT%e(:,:,:)  = dIN%e(:,:,:)
!      case('d')
!         dOUT%d(:,:,:)  = dIN%d(:,:,:)
!      case('v')
!         dOUT%v(:,:,:)  = dIN%v(:,:,:)
!      end select
!   case ("nonsu2")
!      select case(type)
!      case('e')
!         dOUT%e(:,:,:)  = dIN%e(:,:,:)
!      case('v')
!         dOUT%v(:,:,:)  = dIN%v(:,:,:)
!      case('u')
!         dOUT%u(:,:,:)  = dIN%u(:,:,:)
!      end select
!   end select
!   call get_dmft_bath(dOUT,bathOUT)
!   call deallocate_dmft_bath(dIN)
!   call deallocate_dmft_bath(dOUT)
! end subroutine copy_bath_component
