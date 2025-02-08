!
!Rank _nX refers here to the rank of Self WITHOUT the frequency dimension
!


subroutine ed_get_gimp_site_n2(self,axis,type,z)
  complex(8),dimension(:,:,:),intent(inout) :: self ! Green's function matrix
  character(len=*),optional                 :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  character(len=*),optional                 :: type ! Can be :f:var:`"n"` for Normal (default), :f:var:`"a"` for anomalous
  complex(8),dimension(:),optional          :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                          :: axis_
  character(len=1)                          :: type_
  complex(8),dimension(:),allocatable       :: z_
  !
  axis_='m';if(present(axis))axis_=trim(axis)
  type_='n';if(present(type))type_=trim(type)
  !
  !
  call allocate_grids
  !
  if(present(z))then
     allocate(z_, source=z)
  else
     select case(axis_)
     case default;stop "ed_get_gimp ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        z_ = dcmplx(0d0,wm)
     case ('r','R')
        z_ = dcmplx(wr,eps)
     end select
  endif
  !
  L = size(z_)
  !
  call assert_shape(self,[Nspin*Norb,Nspin*Norb,L],'ed_get_gimp','self')
  !
  select case(type_)
  case default; stop "ed_get_gimp ERROR: type is neither Normal, nor Anomalous"
  case ('n','N')
     self = nn2so_reshape( get_G_impurity(z_,axis_), Nspin,Norb,L)
  case('a','A')
     self = nn2so_reshape( get_F_impurity(z_,axis_), Nspin,Norb,L)
  end select
  !
  call deallocate_grids
  !
end subroutine ed_get_gimp_site_n2


subroutine ed_get_gimp_site_n4(self,axis,type,z)
  complex(8),dimension(:,:,:,:,:),intent(inout) :: self
  character(len=*),optional                     :: axis
  character(len=*),optional                     :: type
  complex(8),dimension(:),optional              :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                              :: axis_
  character(len=1)                              :: type_
  complex(8),dimension(:),allocatable           :: z_
  !
  !
  axis_='m';if(present(axis))axis_=trim(axis)
  type_='n';if(present(type))type_=trim(type)
  !
  !
  call allocate_grids
  !
  if(present(z))then
     allocate(z_, source=z)
  else
     select case(axis_)
     case default;stop "ed_get_gimp ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        z_ = dcmplx(0d0,wm)
     case ('r','R')
        z_ = dcmplx(wr,eps)
     end select
  endif
  !
  L = size(z_)
  !
  call assert_shape(self,[Nspin,Nspin,Norb,Norb,L],'ed_get_gimp','self')
  !  
  select case(type_)
  case default; stop "ed_get_gimp ERROR: type is neither Normal, nor Anomalous"
  case ('n','N')
     self = get_G_impurity(z_,axis_)
  case('a','A')
     self = get_F_impurity(z_,axis_)
  end select
  !
  call deallocate_grids
  !
end subroutine ed_get_gimp_site_n4


!##################################################################
!##################################################################
!##################################################################


subroutine ed_get_gimp_lattice_n3(self,nlat,axis,type,z)
  complex(8),dimension(:,:,:),intent(inout) :: self
  integer,intent(in)                        :: nlat  ! Number of inequivalent impurity sites for real-space DMFT
  character(len=*),optional                 :: axis
  character(len=*),optional                 :: type
  complex(8),dimension(:),optional          :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                          :: axis_
  character(len=1)                          :: type_
  complex(8),dimension(:),allocatable       :: z_
  integer                                   :: ilat

  !
  axis_='m';if(present(axis))axis_=trim(axis)
  type_='n';if(present(type))type_=trim(type)
  !
  !
  call allocate_grids
  !
  if(present(z))then
     allocate(z_, source=z)
  else
     select case(axis_)
     case default;stop "ed_get_gimp ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        z_ = dcmplx(0d0,wm)
     case ('r','R')
        z_ = dcmplx(wr,eps)
     end select
  endif
  !
  L = size(z_)
  !


  if(Nlat/=size(Gmats_ineq,1))stop "ERROR ed_get_gimp: wrong Nlat"
  select case(type_)
  case default; stop "ed_get_gimp ERROR: type is neither Normal, nor Anomalous"
  case ('n','N')
     call assert_shape(self,[Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lreal],'ed_get_gimp','self')
     self = nnn2lso_reshape(Freal_ineq,Nlat,Nspin,Norb,Lreal)
  case('a','A')
     call assert_shape(self,[Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lreal],'ed_get_gimp','self')
     self = nnn2lso_reshape(Freal_ineq,Nlat,Nspin,Norb,Lreal)
  end select
end subroutine ed_get_gimp_lattice_n3

subroutine ed_get_gimp_lattice_n4(self,nlat,axis,type)
  complex(8),dimension(:,:,:,:),intent(inout) :: self
  integer,intent(in)                     :: nlat
  character(len=*),optional              :: axis
  character(len=*),optional              :: type
  character(len=1)                       :: axis_
  character(len=1)                       :: type_
  integer                                :: ilat
  axis_='m';if(present(axis))axis_=trim(axis)
  type_='n';if(present(type))type_=trim(type)
  if(Nlat/=size(Gmats_ineq,1))stop "ERROR ed_get_gimp: wrong Nlat"
  select case(type_)
  case default; stop "ed_get_gimp ERROR: type is neither Normal, nor Anomalous"
  case ('n','N')
     select case(axis_)
     case default;stop "ed_get_gimp ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        call assert_shape(self,[Nlat,Nspin*Norb,Nspin*Norb,Lmats],'ed_get_gimp','self')
        do ilat=1,Nlat
           self(ilat,:,:,:) = nn2so_reshape(Gmats_ineq(ilat,:,:,:,:,:),Nspin,Norb,Lmats)
        enddo
     case('r','R')
        call assert_shape(self,[Nlat,Nspin*Norb,Nspin*Norb,Lreal],'ed_get_gimp','self')
        do ilat=1,Nlat
           self(ilat,:,:,:) = nn2so_reshape(Freal_ineq(ilat,:,:,:,:,:),Nspin,Norb,Lreal)
        enddo
     end select
  case('a','A')
     select case(axis_)
     case default;stop "ed_get_gimp ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        call assert_shape(self,[Nlat,Nspin*Norb,Nspin*Norb,Lmats],'ed_get_gimp','self')
        do ilat=1,Nlat
           self(ilat,:,:,:) = nn2so_reshape(Fmats_ineq(ilat,:,:,:,:,:),Nspin,Norb,Lmats)
        enddo
     case('r','R')
        call assert_shape(self,[Nlat,Nspin*Norb,Nspin*Norb,Lreal],'ed_get_gimp','self')
        do ilat=1,Nlat
           self(ilat,:,:,:) = nn2so_reshape(Freal_ineq(ilat,:,:,:,:,:),Nspin,Norb,Lreal)
        enddo
     end select
  end select
end subroutine ed_get_gimp_lattice_n4

subroutine ed_get_gimp_lattice_n6(self,nlat,axis,type)
  complex(8),dimension(:,:,:,:,:,:),intent(inout) :: self
  integer,intent(in)                              :: nlat
  character(len=*),optional                       :: axis
  character(len=*),optional                       :: type
  character(len=1)                                :: axis_
  character(len=1)                                :: type_
  integer                                         :: ilat
  axis_='m';if(present(axis))axis_=trim(axis)
  type_='n';if(present(type))type_=trim(type)
  if(Nlat/=size(Gmats_ineq,1))stop "ERROR ed_get_gimp: wrong Nlat"
  select case(type_)
  case default; stop "ed_get_gimp ERROR: type is neither Normal, nor Anomalous"
  case ('n','N')
     select case(axis_)
     case default;stop "ed_get_gimp ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        call assert_shape(self,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],'ed_get_gimp','self')
        self = Gmats_ineq
     case('r','R')
        call assert_shape(self,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],'ed_get_gimp','self')
        self = Freal_ineq
     end select
  case('a','A')
     select case(axis_)
     case default;stop "ed_get_gimp ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        call assert_shape(self,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],'ed_get_gimp','self')
        self = Fmats_ineq
     case('r','R')
        call assert_shape(self,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],'ed_get_gimp','self')
        self = Freal_ineq
     end select
  end select
end subroutine ed_get_gimp_lattice_n6



