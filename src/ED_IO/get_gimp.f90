!
!Rank _nX refers here to the rank of Self WITHOUT the frequency dimension
!


subroutine ed_get_gimp_site_n2(self,axis,type,z)
  complex(8),dimension(:,:,:),intent(inout)   :: self ! Green's function matrix
  character(len=*),optional                   :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  character(len=*),optional                   :: type ! Can be :f:var:`"n"` for Normal (default), :f:var:`"a"` for anomalous
  complex(8),dimension(:),optional            :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                            :: axis_
  character(len=1)                            :: type_
  complex(8),dimension(:),allocatable         :: z_
  complex(8),dimension(:,:,:,:,:),allocatable :: gf
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
  allocate(gf(Nspin,Nspin,Norb,Norb,L))
  !
  select case(type_)
  case default; stop "ed_get_gimp ERROR: type is neither Normal, nor Anomalous"
  case ('n','N');gf = get_impG(z_,axis_)
  case ('a','A');gf = get_impF(z_,axis_)
  end select
  !
  self = nn2so_reshape( gf, Nspin,Norb,L)
  !
  call deallocate_grids
  !
end subroutine ed_get_gimp_site_n2


subroutine ed_get_gimp_site_n4(self,axis,type,z)
  complex(8),dimension(:,:,:,:,:),intent(inout) :: self
  character(len=*),optional                   :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  character(len=*),optional                   :: type ! Can be :f:var:`"n"` for Normal (default), :f:var:`"a"` for anomalous
  complex(8),dimension(:),optional            :: z    ! User provided array of complex frequency where to evaluate Self
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
     self = get_impG(z_,axis_)
  case('a','A')
     self = get_impF(z_,axis_)
  end select
  !
  call deallocate_grids
  !
end subroutine ed_get_gimp_site_n4


!##################################################################
!LATTICE EXTENSION: need to read to the 
!##################################################################


subroutine ed_get_gimp_lattice_n2(self,nlat,axis,type,z)
  complex(8),dimension(:,:,:),intent(inout)     :: self !! [Nlso,Nlso,:]
  integer,intent(in)                            :: nlat  ! Number of inequivalent impurity sites for real-space DMFT
  character(len=*),optional                     :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  character(len=*),optional                     :: type ! Can be :f:var:`"n"` for Normal (default), :f:var:`"a"` for anomalous
  complex(8),dimension(:),optional              :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                              :: axis_
  character(len=1)                              :: type_
  complex(8),dimension(:),allocatable           :: z_
  integer                                       :: ilat
  complex(8),dimension(:,:,:,:,:,:),allocatable :: gf
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
  call assert_shape(self,[Nlat*Nspin*Norb,Nlat*Nspin*Norb,L],'ed_get_gimp','self')
  !
  allocate(gf(Nlat,Nspin,Nspin,Norb,Norb,L))
  gf = zero
  !
  do ilat=1,Nlat
     call ed_set_suffix(ilat)
     call read_impGmatrix()
     select case(type_)
     case default; stop "ed_get_gimp ERROR: type is neither Normal, nor Anomalous"
     case ('n','N');gf(ilat,:,:,:,:,:) = get_impG(z_,axis_)
     case ('a','A');gf(ilat,:,:,:,:,:) = get_impF(z_,axis_)
     end select
  enddo
  !
  self = nnn2lso_reshape(gf,Nlat,Nspin,Norb,Lreal)
  !
  call ed_reset_suffix()
  call deallocate_grids()
  deallocate(gf)
  if(allocated(impGmatrix))call deallocate_GFmatrix(impGmatrix)
  if(allocated(impGmatrix))deallocate(impGmatrix)
  !
end subroutine ed_get_gimp_lattice_n2

subroutine ed_get_gimp_lattice_n4(self,nlat,axis,type,z)
  complex(8),dimension(:,:,:,:),intent(inout) :: self !! [Nlat,Nso,Nso,:]
  integer,intent(in)                          :: nlat
  character(len=*),optional                   :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  character(len=*),optional                   :: type ! Can be :f:var:`"n"` for Normal (default), :f:var:`"a"` for anomalous
  complex(8),dimension(:),optional            :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                            :: axis_
  character(len=1)                            :: type_
  complex(8),dimension(:),allocatable         :: z_
  integer                                     :: ilat
  complex(8),dimension(:,:,:,:,:),allocatable :: gf
  !
  axis_='m';if(present(axis))axis_=trim(axis)
  type_='n';if(present(type))type_=trim(type)
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
  call assert_shape(self,[Nlat,Nspin*Norb,Nspin*Norb,L],'ed_get_gimp','self')
  !
  allocate(gf(Nspin,Nspin,Norb,Norb,L))
  gf = zero

  do ilat=1,Nlat
     call ed_set_suffix(ilat)
     call read_impGmatrix()
     !
     select case(type_)
     case default; stop "ed_get_gimp ERROR: type is neither Normal, nor Anomalous"
     case ('n','N');gf = get_impG(z_,axis_)
     case ('a','A');gf = get_impF(z_,axis_)
     end select
     !
     self(ilat,:,:,:) = nn2so_reshape(gf,Nspin,Norb,L)
     !
  enddo
  !
  call ed_reset_suffix()
  call deallocate_grids()
  deallocate(gf)
  if(allocated(impGmatrix))call deallocate_GFmatrix(impGmatrix)
  if(allocated(impGmatrix))deallocate(impGmatrix)
  !
end subroutine ed_get_gimp_lattice_n4

subroutine ed_get_gimp_lattice_n6(self,nlat,axis,type,z)
  complex(8),dimension(:,:,:,:,:,:),intent(inout) :: self
  integer,intent(in)                              :: nlat
  character(len=*),optional                       :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  character(len=*),optional                       :: type ! Can be :f:var:`"n"` for Normal (default), :f:var:`"a"` for anomalous
  complex(8),dimension(:),optional                :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                                :: axis_
  character(len=1)                                :: type_
  complex(8),dimension(:),allocatable             :: z_
  integer                                         :: ilat
  !
  axis_='m';if(present(axis))axis_=trim(axis)
  type_='n';if(present(type))type_=trim(type)
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
  call assert_shape(self,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],'ed_get_gimp','self')
  !
  do ilat=1,Nlat
     call ed_set_suffix(ilat)
     call read_impGmatrix()
     select case(type_)
     case default; stop "ed_get_gimp ERROR: type is neither Normal, nor Anomalous"
     case ('n','N');self(ilat,:,:,:,:,:) = get_impG(z_,axis_)
     case ('a','A');self(ilat,:,:,:,:,:) = get_impF(z_,axis_)
     end select
  enddo
  !
  call ed_reset_suffix()
  call deallocate_grids()
  if(allocated(impGmatrix))call deallocate_GFmatrix(impGmatrix)
  if(allocated(impGmatrix))deallocate(impGmatrix)
  !
end subroutine ed_get_gimp_lattice_n6



