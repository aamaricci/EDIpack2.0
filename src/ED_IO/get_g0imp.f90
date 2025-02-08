!
!Rank _nX refers here to the rank of Self WITHOUT the frequency dimension
!


subroutine ed_get_g0imp_site_n2(self,bath,axis,type,z)
  complex(8),dimension(:,:,:),intent(inout)   :: self ! Green's function matrix
  real(8),dimension(:)                        :: bath ! The bath vector
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
     case default;stop "ed_get_g0imp ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        z_ = dcmplx(0d0,wm)
     case ('r','R')
        z_ = dcmplx(wr,eps)
     end select
  endif
  !
  L = size(z_)
  !
  call assert_shape(self,[Nspin*Norb,Nspin*Norb,L],'ed_get_g0imp','self')
  !
  allocate(gf(Nspin,Nspin,Norb,Norb,L))
  !
  call ed_get_g0and(z_,bath,gf,axis=axis_,type=type_)
  !
  self = nn2so_reshape( gf, Nspin,Norb,L)
  !
  call deallocate_grids
  !
end subroutine ed_get_g0imp_site_n2


subroutine ed_get_g0imp_site_n4(self,bath,axis,type,z)
  complex(8),dimension(:,:,:,:,:),intent(inout) :: self
  real(8),dimension(:)                          :: bath ! The bath vector
  character(len=*),optional                     :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  character(len=*),optional                     :: type ! Can be :f:var:`"n"` for Normal (default), :f:var:`"a"` for anomalous
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
     case default;stop "ed_get_g0imp ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        z_ = dcmplx(0d0,wm)
     case ('r','R')
        z_ = dcmplx(wr,eps)
     end select
  endif
  !
  L = size(z_)
  !
  call assert_shape(self,[Nspin,Nspin,Norb,Norb,L],'ed_get_g0imp','self')
  !
  call ed_get_g0and(z_,bath,self,axis=axis_,type=type_)
  !
  call deallocate_grids
  !
end subroutine ed_get_g0imp_site_n4


!##################################################################
!LATTICE EXTENSION: need to read to the 
!##################################################################


subroutine ed_get_g0imp_lattice_n2(self,bath,axis,type,z)
  complex(8),dimension(:,:,:),intent(inout)     :: self !! [Nlso,Nlso,:]
  real(8),dimension(:,:)                        :: bath
  character(len=*),optional                     :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  character(len=*),optional                     :: type ! Can be :f:var:`"n"` for Normal (default), :f:var:`"a"` for anomalous
  complex(8),dimension(:),optional              :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                              :: axis_
  character(len=1)                              :: type_
  complex(8),dimension(:),allocatable           :: z_
  integer                                       :: ilat,Nlat
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
     case default;stop "ed_get_g0imp ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        z_ = dcmplx(0d0,wm)
     case ('r','R')
        z_ = dcmplx(wr,eps)
     end select
  endif
  !
  L = size(z_)
  Nlat = size(bath,1)
  call assert_shape(self,[Nlat*Nspin*Norb,Nlat*Nspin*Norb,L],'ed_get_g0imp','self')
  !
  allocate(gf(Nlat,Nspin,Nspin,Norb,Norb,L))
  gf = zero
  !
  do ilat=1,Nlat
     call ed_set_suffix(ilat)
     call ed_get_g0and(z_,bath(ilat,:),gf(ilat,:,:,:,:,:),axis=axis_,type=type_)
  enddo
  !
  self = nnn2lso_reshape(gf,Nlat,Nspin,Norb,Lreal)
  !
  call ed_reset_suffix()
  call deallocate_grids()
  deallocate(gf)
  !
end subroutine ed_get_g0imp_lattice_n2

subroutine ed_get_g0imp_lattice_n4(self,bath,axis,type,z)
  complex(8),dimension(:,:,:,:),intent(inout) :: self !! [Nlat,Nso,Nso,:]
  real(8),dimension(:,:)                      :: bath
  character(len=*),optional                   :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  character(len=*),optional                   :: type ! Can be :f:var:`"n"` for Normal (default), :f:var:`"a"` for anomalous
  complex(8),dimension(:),optional            :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                            :: axis_
  character(len=1)                            :: type_
  complex(8),dimension(:),allocatable         :: z_
  integer                                     :: ilat,Nlat
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
     case default;stop "ed_get_g0imp ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        z_ = dcmplx(0d0,wm)
     case ('r','R')
        z_ = dcmplx(wr,eps)
     end select
  endif
  !
  L    = size(z_)
  Nlat = size(bath,1)
  !
  call assert_shape(self,[Nlat,Nspin*Norb,Nspin*Norb,L],'ed_get_g0imp','self')
  !
  allocate(gf(Nspin,Nspin,Norb,Norb,L))
  gf = zero

  do ilat=1,Nlat
     call ed_set_suffix(ilat)
     call ed_get_g0and(z_,bath(ilat,:),gf,axis=axis_,type=type_)
     !
     self(ilat,:,:,:) = nn2so_reshape(gf,Nspin,Norb,L)
     !
  enddo
  !
  call ed_reset_suffix()
  call deallocate_grids()
  deallocate(gf)
  !
end subroutine ed_get_g0imp_lattice_n4

subroutine ed_get_g0imp_lattice_n6(self,bath,axis,type,z)
  complex(8),dimension(:,:,:,:,:,:),intent(inout) :: self
  real(8),dimension(:,:)                          :: bath
  character(len=*),optional                       :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  character(len=*),optional                       :: type ! Can be :f:var:`"n"` for Normal (default), :f:var:`"a"` for anomalous
  complex(8),dimension(:),optional                :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                                :: axis_
  character(len=1)                                :: type_
  complex(8),dimension(:),allocatable             :: z_
  integer                                         :: ilat,Nlat
  !
  axis_='m';if(present(axis))axis_=trim(axis)
  type_='n';if(present(type))type_=trim(type)
  call allocate_grids
  !
  if(present(z))then
     allocate(z_, source=z)
  else
     select case(axis_)
     case default;stop "ed_get_g0imp ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        z_ = dcmplx(0d0,wm)
     case ('r','R')
        z_ = dcmplx(wr,eps)
     end select
  endif
  !
  L = size(z_)
  Nlat=size(bath,1)
  !
  call assert_shape(self,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],'ed_get_g0imp','self')
  !
  do ilat=1,Nlat
     call ed_set_suffix(ilat)
     call ed_get_g0and(z_,bath(ilat,:),self(ilat,:,:,:,:,:),axis=axis_,type=type_)
  enddo
  !
  call ed_reset_suffix()
  call deallocate_grids()
  !
end subroutine ed_get_g0imp_lattice_n6



