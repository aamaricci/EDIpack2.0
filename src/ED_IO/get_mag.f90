subroutine ed_get_mag_1(mag,component) 
  real(8),dimension(Norb)   :: mag
  character(len=1),optional :: component
  character(len=1)          :: char
  integer                   :: id
  char='Z';if(present(component))char=component
  select case(char)
  case default;id=3
  case('X','x');id=1
  case('Y','y');id=2
  end select
  mag = ed_mag(id,:)
end subroutine ed_get_mag_1

subroutine ed_get_mag_2(mag,iorb,component) 
  real(8)   :: mag
  integer   :: iorb
  character(len=1),optional :: component
  character(len=1)          :: char
  integer                   :: id
  char='Z';if(present(component))char=component
  select case(char)
  case default;id=3
  case('X','x');id=1
  case('Y','y');id=2
  end select
  if(iorb>Norb)stop "ed_get_dens_up error: orbital index > N_orbital"
  mag = ed_mag(id,iorb)
end subroutine ed_get_mag_2

subroutine ed_get_mag_lattice_1(yii,Nlat,component)
  integer                      :: Nlat
  real(8),dimension(Nlat,Norb) :: yii
  character(len=1),optional :: component
  character(len=1)          :: char
  integer                   :: id
  char='Z';if(present(component))char=component
  select case(char)
  case default;id=3
  case('X','x');id=1
  case('Y','y');id=2
  end select
  yii=0d0
  if(allocated(mag_ineq))then
     if(Nlat>size(mag_ineq,1)) stop "ed_get_mag error: required N_sites > evaluated N_sites"
     yii=mag_ineq(:Nlat,id,:)
  endif
end subroutine ed_get_mag_lattice_1

subroutine ed_get_mag_lattice_2(yii,Nlat,iorb,component) 
  integer                 :: Nlat
  real(8),dimension(Nlat) :: yii
  integer                 :: iorb
  character(len=1),optional :: component
  character(len=1)          :: char
  integer                   :: id
  char='Z';if(present(component))char=component
  select case(char)
  case default;id=3
  case('X','x');id=1
  case('Y','y');id=2
  end select
  if(iorb>Norb)stop "ed_get_mag error: orbital index > N_orbital"
  yii=0d0
  if(allocated(mag_ineq))then
     if(Nlat>size(mag_ineq,1)) stop "ed_get_mag error: required N_sites > evaluated N_sites"
     yii=mag_ineq(:Nlat,id,iorb)
  endif
end subroutine ed_get_mag_lattice_2
