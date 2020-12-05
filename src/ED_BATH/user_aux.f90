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
  if(size(lambda_impHloc) .ne. dmft_bath_%Nbasis)then
     dmft_bath_%item(ibath)%lambda(dmft_bath_%Nbasis)=offset
  else
     do isym=1,size(lambda_impHloc)
        if(is_identity(H_basis(isym)%O)) dmft_bath_%item(ibath)%lambda(isym)=offset
        return
     enddo
  endif
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
     ed_file_suffix=reg(ineq_site_suffix)//reg(txtfy(ilat,site_indx_padding))
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
subroutine spin_symmetrize_bath_lattice(bath_,save)
  real(8),dimension(:,:) :: bath_
  logical,optional       :: save
  logical                :: save_
  integer                :: Nsites,ilat
  save_=.true.;if(present(save))save_=save
  Nsites=size(bath_,1)
  do ilat=1,Nsites
     ed_file_suffix=reg(ineq_site_suffix)//reg(txtfy(ilat,site_indx_padding))
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
  save_=.true.;if(present(save))save_=save
  Nsites=size(bath_,1)
  do ilat=1,Nsites
     ed_file_suffix=reg(ineq_site_suffix)//reg(txtfy(ilat,site_indx_padding))
     call orb_symmetrize_bath_site(bath_(ilat,:),save_)
  enddo
  ed_file_suffix=""
end subroutine orb_symmetrize_bath_lattice

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
  indx_=1     ;if(present(indx))indx_=indx
  save_=.true.;if(present(save))save_=save
  Nsites=size(bath_,1)
  do ilat=1,Nsites
     ed_file_suffix=reg(ineq_site_suffix)//reg(txtfy(ilat,site_indx_padding))
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
  save_=.true.;if(present(save))save_=save
  Nsites=size(bath_,1)
  do ilat=1,Nsites
     ed_file_suffix=reg(ineq_site_suffix)//reg(txtfy(ilat,site_indx_padding))
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
  save_=.true.;if(present(save))save_=save
  Nsites=size(bath_,1)
  do ilat=1,Nsites
     ed_file_suffix=reg(ineq_site_suffix)//reg(txtfy(ilat,site_indx_padding))
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
     ed_file_suffix=reg(ineq_site_suffix)//reg(txtfy(ilat,site_indx_padding))
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
        if(type/="e".OR.type/='v')stop "check_bath_component error: type!=e,v"
     case ("superc")
        if(type/="e".OR.type/='v'.OR.type/='d')stop "check_bath_component error: type!=e,v,d"
     case ("nonsu2")
        if(type/="e".OR.type/='v'.OR.type/='u')stop "check_bath_component error: type!=e,v,u"
     end select
  case ("replica")
     if(type/="v".OR.type/="l")stop "check_bath_component error: type!=v,l"
  end select
  return
end subroutine check_bath_component


!+-------------------------------------------------------------------+
!PURPOSE  : Inquire the correct bath size to allocate the 
! the bath array in the calling program.
!
! Get size of each dimension of the component array. 
! The Result is an rank 1 integer array Ndim with dimension:
! 3 for get_component_size_bath
! 2 for get_spin_component_size_bath & get_orb_component_size_bath
! 1 for get_spin_orb_component_size_bath
!+-------------------------------------------------------------------+
function get_bath_component_dimension(type) result(Ndim)
  character(len=1) :: type
  integer          :: Ndim(3)
  call  check_bath_component(type)
  select case(bath_type)
  case default
     Ndim=[Nspin,Norb,Nbath]
  case('hybrid')
     select case(ed_mode)
     case default
        select case(type)
        case('e')
           Ndim=[Nspin,1,Nbath]
        case('v')
           Ndim=[Nspin,Norb,Nbath]
        end select
     case ("superc")
        select case(type)
        case('e','d')
           Ndim=[Nspin,1,Nbath]
        case('v')
           Ndim=[Nspin,Norb,Nbath]
        end select
     case ("nonsu2")
        select case(type)
        case('e')
           Ndim=[Nspin,1,Nbath]
        case('v','u')
           Ndim=[Nspin,Norb,Nbath]
        end select
     end select
  end select
end function get_bath_component_dimension


!+-----------------------------------------------------------------------------+!
!PURPOSE: check that the input array hsa the correct dimensions specified 
! for the choice of itype and possiblty ispin and/or iorb.
!+-----------------------------------------------------------------------------+!
subroutine assert_bath_component_size(array,type,string1,string2)
  real(8),dimension(:,:,:) :: array
  character(len=1)         :: type
  character(len=*)         :: string1,string2
  integer                  :: Ndim(3)
  Ndim = get_bath_component_dimension(type)
  call assert_shape(Array,Ndim,reg(string1),reg(string2))
end subroutine assert_bath_component_size








!+-----------------------------------------------------------------------------+!
!PURPOSE: Get a specified itype,ispin,iorb component of the user bath.
! The component is returned into an Array of rank D
! get_full_component_bath    : return the entire itype component (D=3)
! get_spin_component_bath    : return the itype component for the select ispin (D=2)
! get_spin_orb_component_bath: return the itype component for the select ispin & iorb (D=1)
!+-----------------------------------------------------------------------------+!
subroutine get_bath_component(array,bath_,type)
  real(8),dimension(:,:,:) :: array
  real(8),dimension(:)     :: bath_
  character(len=1)         :: type
  logical                  :: check
  type(effective_bath)     :: dmft_bath_
  !
  check= check_bath_dimension(bath_)
  if(.not.check)stop "get_component_bath error: wrong bath dimensions"
  call allocate_dmft_bath(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  call assert_bath_component_size(array,type,"get_bath_component","Array")
  call check_bath_component(type)
  select case(ed_mode)
  case default
     select case(type)
     case('e')
        Array = dmft_bath_%e(:,:,:)
     case('v')
        Array = dmft_bath_%v(:,:,:)
     end select
  case ("superc")
     select case(type)
     case('e')
        Array = dmft_bath_%e(:,:,:)           
     case('d')
        Array = dmft_bath_%d(:,:,:)
     case('v')
        Array = dmft_bath_%v(:,:,:)
     end select
  case ("nonsu2")
     select case(type)
     case('e')
        Array = dmft_bath_%e(:,:,:)           
     case('v')
        Array = dmft_bath_%v(:,:,:)
     case('u')
        Array = dmft_bath_%u(:,:,:)
     end select
  end select
  call deallocate_dmft_bath(dmft_bath_)
end subroutine get_bath_component


!+-----------------------------------------------------------------------------+!
!PURPOSE: Set a specified itype,ispin,iorb component of the user bath.
!+-----------------------------------------------------------------------------+!
subroutine set_bath_component(array,bath_,type)
  real(8),dimension(:,:,:) :: array
  real(8),dimension(:)     :: bath_
  character(len=1)         :: type
  logical                  :: check
  type(effective_bath)     :: dmft_bath_
  !
  check= check_bath_dimension(bath_)
  if(.not.check)stop "set_component_bath error: wrong bath dimensions"
  call allocate_dmft_bath(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  call assert_bath_component_size(array,type,"set_bath_component","Array")
  call check_bath_component(type)
  select case(ed_mode)
  case default
     select case(type)
     case('e')
        dmft_bath_%e(:,:,:) = Array 
     case('v')
        dmft_bath_%v(:,:,:) = Array
     end select
  case ("superc")
     select case(type)
     case('e')
        dmft_bath_%e(:,:,:) = Array   
     case('d')
        dmft_bath_%d(:,:,:) = Array
     case('v')
        dmft_bath_%v(:,:,:) = Array
     end select
  case ("nonsu2")
     select case(type)
     case('e')
        dmft_bath_%e(:,:,:) = Array   
     case('v')
        dmft_bath_%v(:,:,:) = Array
     case('u')
        dmft_bath_%u(:,:,:) = Array
     end select
  end select
  call get_dmft_bath(dmft_bath_,bath_)
  call deallocate_dmft_bath(dmft_bath_)
end subroutine set_bath_component



!+-----------------------------------------------------------------------------+!
!PURPOSE: Copy a specified component of IN bath to the OUT bath.
!+-----------------------------------------------------------------------------+!
subroutine copy_bath_component(bathIN,bathOUT,type)
  real(8),dimension(:)     :: bathIN,bathOUT
  character(len=1)         :: type
  logical                  :: check
  type(effective_bath)     :: dIN,dOUT
  !
  check= check_bath_dimension(bathIN)
  if(.not.check)stop "copy_component_bath error: wrong bath dimensions IN"
  check= check_bath_dimension(bathOUT)
  if(.not.check)stop "copy_component_bath error: wrong bath dimensions OUT"
  call allocate_dmft_bath(dIN)
  call allocate_dmft_bath(dOUT)
  call set_dmft_bath(bathIN,dIN)
  call set_dmft_bath(bathOUT,dOUT)
  call check_bath_component(type)
  select case(ed_mode)
  case default
     select case(type)
     case('e')
        dOUT%e(:,:,:)  = dIN%e(:,:,:)
     case('v')
        dOUT%v(:,:,:)  = dIN%v(:,:,:)
     end select
  case ("superc")
     select case(type)
     case('e')
        dOUT%e(:,:,:)  = dIN%e(:,:,:)
     case('d')
        dOUT%d(:,:,:)  = dIN%d(:,:,:)
     case('v')
        dOUT%v(:,:,:)  = dIN%v(:,:,:)
     end select
  case ("nonsu2")
     select case(type)
     case('e')
        dOUT%e(:,:,:)  = dIN%e(:,:,:)
     case('v')
        dOUT%v(:,:,:)  = dIN%v(:,:,:)
     case('u')
        dOUT%u(:,:,:)  = dIN%u(:,:,:)
     end select
  end select
  call get_dmft_bath(dOUT,bathOUT)
  call deallocate_dmft_bath(dIN)
  call deallocate_dmft_bath(dOUT)
end subroutine copy_bath_component








