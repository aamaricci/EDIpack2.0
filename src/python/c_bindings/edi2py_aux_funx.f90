!SET HLOC
subroutine ed_set_Hloc_single_N2_c(Hloc,d) bind(c, name='ed_set_Hloc_single_N2')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                        :: d(2)
  complex(c_double_complex),dimension(d(1),d(2)),intent(in) :: Hloc
  call ed_set_Hloc(Hloc)
end subroutine ed_set_Hloc_single_N2_c

subroutine ed_set_Hloc_single_N4_c(Hloc,d) bind(c, name='ed_set_Hloc_single_N4')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                                   :: d(4)
  complex(c_double_complex),dimension(d(1),d(2),d(3),d(4)),intent(in)  :: Hloc
  call ed_set_Hloc(Hloc)
end subroutine ed_set_Hloc_single_N4_c

subroutine ed_set_Hloc_lattice_N2_c(Hloc,d,Nlat) bind(c, name='ed_set_Hloc_lattice_N2')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                        :: d(2)
  complex(c_double_complex),dimension(d(1),d(2)),intent(in) :: Hloc
  integer(c_int),value                                      :: Nlat
  call ed_set_Hloc(Hloc,Nlat)
end subroutine ed_set_Hloc_lattice_N2_c


subroutine ed_set_Hloc_lattice_N3_c(Hloc,d,Nlat) bind(c, name='ed_set_Hloc_lattice_N3')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                             :: d(3)
  complex(c_double_complex),dimension(d(1),d(2),d(3)),intent(in) :: Hloc
  integer(c_int),value                                           :: Nlat
  call ed_set_Hloc(Hloc,Nlat)
end subroutine ed_set_Hloc_lattice_N3_c

subroutine ed_set_Hloc_lattice_N5_c(Hloc,d,Nlat) bind(c, name='ed_set_Hloc_lattice_N5')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                                       :: d(5)
  complex(c_double_complex),dimension(d(1),d(2),d(3),d(4),d(5)),intent(in) :: Hloc
  integer(c_int),value                                                     :: Nlat
  call ed_set_Hloc(Hloc,Nlat)
end subroutine ed_set_Hloc_lattice_N5_c



!SEARCH VARIABLE:
subroutine search_variable(var,ntmp,converged) bind(c, name='search_variable')
  use, intrinsic :: iso_c_binding
  real(c_double),dimension(1)         :: var(1)
  real(c_double),dimension(1)         :: ntmp(1)
  integer(c_int),dimension(1)         :: converged(1)
  logical                             :: bool
  converged(1)=0
  call ed_search_variable(var(1),ntmp(1),bool)
  if (bool) converged(1)=1
end subroutine search_variable






