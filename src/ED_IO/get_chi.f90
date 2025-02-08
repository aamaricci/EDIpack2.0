subroutine ed_get_spinChi_site_n2(self,axis,z)
  complex(8),dimension(:,:,:),intent(inout) :: self ! spin susceptibility 
  character(len=*),optional                 :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  complex(8),dimension(:),optional          :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                          :: axis_
  complex(8),dimension(:),allocatable       :: z_
  !
  axis_='m';if(present(axis))axis_=trim(axis)
  !
  !
  call allocate_grids
  !
  if(present(z))then
     allocate(z_, source=z)
  else
     select case(axis_)
     case default;stop "ed_get_dimp ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        z_ = dcmplx(0d0,wm)
     case ('r','R')
        z_ = dcmplx(wr,eps)
     end select
  endif
  !
  L = size(z_)
  call assert_shape(self,[Norb,Norb,L],'ed_get_spinChi','self')
  !
  self = get_spinChi(z_,axis_)
  call deallocate_grids
  !
end subroutine ed_get_spinChi_site_n2


subroutine ed_get_densChi_site_n2(self,axis,z)
  complex(8),dimension(:,:,:),intent(inout) :: self ! spin susceptibility 
  character(len=*),optional                 :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  complex(8),dimension(:),optional          :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                          :: axis_
  complex(8),dimension(:),allocatable       :: z_
  !
  axis_='m';if(present(axis))axis_=trim(axis)
  !
  !
  call allocate_grids
  !
  if(present(z))then
     allocate(z_, source=z)
  else
     select case(axis_)
     case default;stop "ed_get_dimp ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        z_ = dcmplx(0d0,wm)
     case ('r','R')
        z_ = dcmplx(wr,eps)
     end select
  endif
  !
  L = size(z_)
  call assert_shape(self,[Norb,Norb,L],'ed_get_spinChi','self')
  !
  self = get_densChi(z_,axis_)
  call deallocate_grids
  !
end subroutine ed_get_densChi_site_n2


subroutine ed_get_pairChi_site_n2(self,axis,z)
  complex(8),dimension(:,:,:),intent(inout) :: self ! spin susceptibility 
  character(len=*),optional                 :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  complex(8),dimension(:),optional          :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                          :: axis_
  complex(8),dimension(:),allocatable       :: z_
  !
  axis_='m';if(present(axis))axis_=trim(axis)
  !
  !
  call allocate_grids
  !
  if(present(z))then
     allocate(z_, source=z)
  else
     select case(axis_)
     case default;stop "ed_get_dimp ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        z_ = dcmplx(0d0,wm)
     case ('r','R')
        z_ = dcmplx(wr,eps)
     end select
  endif
  !
  L = size(z_)
  call assert_shape(self,[Norb,Norb,L],'ed_get_spinChi','self')
  !
  self = get_pairChi(z_,axis_)
  call deallocate_grids
  !
end subroutine ed_get_pairChi_site_n2

subroutine ed_get_exctChi_site_n2(self,axis,z)
  complex(8),dimension(:,:,:,:),intent(inout) :: self ! spin susceptibility 
  character(len=*),optional                   :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  complex(8),dimension(:),optional            :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                            :: axis_
  complex(8),dimension(:),allocatable         :: z_
  !
  axis_='m';if(present(axis))axis_=trim(axis)
  !
  !
  call allocate_grids
  !
  if(present(z))then
     allocate(z_, source=z)
  else
     select case(axis_)
     case default;stop "ed_get_dimp ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        z_ = dcmplx(0d0,wm)
     case ('r','R')
        z_ = dcmplx(wr,eps)
     end select
  endif
  !
  L = size(z_)
  call assert_shape(self,[3,Norb,Norb,L],'ed_get_spinChi','self')
  !
  self = get_exctChi(z_,axis_)
  call deallocate_grids
  !
end subroutine ed_get_exctChi_site_n2




!##################################################################
!LATTICE EXTENSION: need to read to the 
!##################################################################


subroutine ed_get_spinChi_lattice_n2(self,nlat,axis,z)
  complex(8),dimension(:,:,:,:),intent(inout)       :: self !! [Nlat,Norb,Norb,:]
  integer,intent(in)                            :: nlat  ! Number of inequivalent impurity sites for real-space DMFT
  character(len=*),optional                     :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  complex(8),dimension(:),optional              :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                              :: axis_
  complex(8),dimension(:),allocatable           :: z_
  integer                                       :: ilat
  !
  axis_='m';if(present(axis))axis_=trim(axis)
  !
  call allocate_grids
  !
  if(present(z))then
     allocate(z_, source=z)
  else
     select case(axis_)
     case default;stop "ed_get_dimp ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        z_ = dcmplx(0d0,wm)
     case ('r','R')
        z_ = dcmplx(wr,eps)
     end select
  endif
  !
  L = size(z_)
  !
  call assert_shape(self,[Nlat,Norb,Norb,L],'ed_get_dimp','self')
  !
  do ilat=1,Nlat
     call ed_set_suffix(ilat)
     call read_spinChimatrix()
     self(ilat,:,:,:) =  get_spinChi(z_,axis_)
  enddo
  !
  call ed_reset_suffix()
  call deallocate_grids()
  if(allocated(spinChimatrix))call deallocate_GFmatrix(spinChimatrix)
  if(allocated(spinChimatrix))deallocate(spinChimatrix)
  !
end subroutine ed_get_spinChi_lattice_n2





subroutine ed_get_densChi_lattice_n2(self,nlat,axis,z)
  complex(8),dimension(:,:,:,:),intent(inout)       :: self !! [Nlat,Norb,Norb,:]
  integer,intent(in)                            :: nlat  ! Number of inequivalent impurity sites for real-space DMFT
  character(len=*),optional                     :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  complex(8),dimension(:),optional              :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                              :: axis_
  complex(8),dimension(:),allocatable           :: z_
  integer                                       :: ilat
  !
  axis_='m';if(present(axis))axis_=trim(axis)
  !
  call allocate_grids
  !
  if(present(z))then
     allocate(z_, source=z)
  else
     select case(axis_)
     case default;stop "ed_get_dimp ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        z_ = dcmplx(0d0,wm)
     case ('r','R')
        z_ = dcmplx(wr,eps)
     end select
  endif
  !
  L = size(z_)
  !
  call assert_shape(self,[Nlat,Norb,Norb,L],'ed_get_dimp','self')
  !
  do ilat=1,Nlat
     call ed_set_suffix(ilat)
     call read_densChimatrix()
     self(ilat,:,:,:) =  get_densChi(z_,axis_)
  enddo
  !
  call ed_reset_suffix()
  call deallocate_grids()
  if(allocated(densChimatrix))call deallocate_GFmatrix(densChimatrix)
  if(allocated(densChimatrix))deallocate(densChimatrix)
  !
end subroutine ed_get_densChi_lattice_n2



subroutine ed_get_pairChi_lattice_n2(self,nlat,axis,z)
  complex(8),dimension(:,:,:,:),intent(inout)       :: self !! [Nlat,Norb,Norb,:]
  integer,intent(in)                            :: nlat  ! Number of inequivalent impurity sites for real-space DMFT
  character(len=*),optional                     :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  complex(8),dimension(:),optional              :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                              :: axis_
  complex(8),dimension(:),allocatable           :: z_
  integer                                       :: ilat
  !
  axis_='m';if(present(axis))axis_=trim(axis)
  !
  call allocate_grids
  !
  if(present(z))then
     allocate(z_, source=z)
  else
     select case(axis_)
     case default;stop "ed_get_dimp ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        z_ = dcmplx(0d0,wm)
     case ('r','R')
        z_ = dcmplx(wr,eps)
     end select
  endif
  !
  L = size(z_)
  !
  call assert_shape(self,[Nlat,Norb,Norb,L],'ed_get_dimp','self')
  !
  do ilat=1,Nlat
     call ed_set_suffix(ilat)
     call read_pairChimatrix()
     self(ilat,:,:,:) =  get_pairChi(z_,axis_)
  enddo
  !
  call ed_reset_suffix()
  call deallocate_grids()
  if(allocated(pairChimatrix))call deallocate_GFmatrix(pairChimatrix)
  if(allocated(pairChimatrix))deallocate(pairChimatrix)
  !
end subroutine ed_get_pairChi_lattice_n2




subroutine ed_get_exctChi_lattice_n2(self,nlat,axis,z)
  complex(8),dimension(:,:,:,:,:),intent(inout)       :: self !! [Nlat,3,Norb,Norb,:]
  integer,intent(in)                            :: nlat  ! Number of inequivalent impurity sites for real-space DMFT
  character(len=*),optional                     :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  complex(8),dimension(:),optional              :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                              :: axis_
  complex(8),dimension(:),allocatable           :: z_
  integer                                       :: ilat
  !
  axis_='m';if(present(axis))axis_=trim(axis)
  !
  call allocate_grids
  !
  if(present(z))then
     allocate(z_, source=z)
  else
     select case(axis_)
     case default;stop "ed_get_dimp ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        z_ = dcmplx(0d0,wm)
     case ('r','R')
        z_ = dcmplx(wr,eps)
     end select
  endif
  !
  L = size(z_)
  !
  call assert_shape(self,[Nlat,3,Norb,Norb,L],'ed_get_dimp','self')
  !
  do ilat=1,Nlat
     call ed_set_suffix(ilat)
     call read_exctChimatrix()
     self(ilat,:,:,:,:) =  get_exctChi(z_,axis_)
  enddo
  !
  call ed_reset_suffix()
  call deallocate_grids()
  if(allocated(exctChimatrix))call deallocate_GFmatrix(exctChimatrix)
  if(allocated(exctChimatrix))deallocate(exctChimatrix)
  !
end subroutine ed_get_exctChi_lattice_n2
