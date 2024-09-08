!ED_IO:

!Sigma
subroutine ed_get_sigma_site_n3_c(self,d,axis,typ) bind(c, name='ed_get_sigma_site_n3')
  integer(c_int64_t)                                                :: d(3)
  complex(c_double_complex),dimension(d(1),d(2),d(3)),intent(inout) :: self
  character(kind=c_char), dimension(1),optional                     :: axis,typ
  character(len=1)                                                  :: axis_,typ_
  typ_(1:1)=typ(1)
  axis_(1:1)=axis(1)
  call ed_get_sigma(self,axis_,typ_)
end subroutine ed_get_sigma_site_n3_c

subroutine ed_get_sigma_site_n5_c(self,d,axis,typ) bind(c, name='ed_get_sigma_site_n5')
  integer(c_int64_t)                                                          :: d(5)
  complex(c_double_complex),dimension(d(1),d(2),d(3),d(4),d(5)),intent(inout) :: self
  character(kind=c_char), dimension(1),optional                               :: axis,typ
  character(len=1)                                                            :: axis_,typ_
  typ_(1:1)=typ(1)
  axis_(1:1)=axis(1)
  call ed_get_sigma(self,axis_,typ_)
end subroutine ed_get_sigma_site_n5_c

subroutine ed_get_sigma_lattice_n3_c(self,d,nlat,axis,typ) bind(c, name='ed_get_sigma_lattice_n3')
  integer(c_int64_t)                                                          :: d(3)
  complex(c_double_complex),dimension(d(1),d(2),d(3)),intent(inout)           :: self
  integer(c_int),value                                                        :: nlat
  character(kind=c_char), dimension(1),optional                               :: axis,typ
  character(len=1)                                                            :: axis_,typ_
  typ_(1:1)=typ(1)
  axis_(1:1)=axis(1)
  call ed_get_sigma(self,nlat,axis_,typ_)
end subroutine ed_get_sigma_lattice_n3_c

subroutine ed_get_sigma_lattice_n4_c(self,d,nlat,axis,typ) bind(c, name='ed_get_sigma_lattice_n4')
  integer(c_int64_t)                                                          :: d(4)
  complex(c_double_complex),dimension(d(1),d(2),d(3),d(4)),intent(inout)      :: self
  integer(c_int),value                                                        :: nlat
  character(kind=c_char), dimension(1),optional                               :: axis,typ
  character(len=1)                                                            :: axis_,typ_
  typ_(1:1)=typ(1)
  axis_(1:1)=axis(1)
  call ed_get_sigma(self,nlat,axis_,typ_)
end subroutine ed_get_sigma_lattice_n4_c

subroutine ed_get_sigma_lattice_n6_c(self,nlat,d,axis,typ) bind(c, name='ed_get_sigma_lattice_n6')
  integer(c_int64_t)                                                                    :: d(6)
  complex(c_double_complex),dimension(d(1),d(2),d(3),d(4),d(5),d(6)),intent(inout)      :: self
  integer(c_int),value                                                                  :: nlat
  character(kind=c_char), dimension(1),optional                                         :: axis,typ
  character(len=1)                                                                      :: axis_,typ_
  typ_(1:1)=typ(1)
  axis_(1:1)=axis(1)
  call ed_get_sigma(self,nlat,axis_,typ_)
end subroutine ed_get_sigma_lattice_n6_c

!Gimp
subroutine ed_get_gimp_site_n3_c(self,d,axis,typ) bind(c, name='ed_get_gimp_site_n3')
  integer(c_int64_t)                                                :: d(3)
  complex(c_double_complex),dimension(d(1),d(2),d(3)),intent(inout) :: self
  character(kind=c_char), dimension(1),optional                     :: axis,typ
  character(len=1)                                                  :: axis_,typ_
  typ_(1:1)=typ(1)
  axis_(1:1)=axis(1)
  call ed_get_gimp(self,axis_,typ_)
end subroutine ed_get_gimp_site_n3_c

subroutine ed_get_gimp_site_n5_c(self,d,axis,typ) bind(c, name='ed_get_gimp_site_n5')
  integer(c_int64_t)                                                          :: d(5)
  complex(c_double_complex),dimension(d(1),d(2),d(3),d(4),d(5)),intent(inout) :: self
  character(kind=c_char), dimension(1),optional                               :: axis,typ
  character(len=1)                                                            :: axis_,typ_
  typ_(1:1)=typ(1)
  axis_(1:1)=axis(1)
  call ed_get_gimp(self,axis_,typ_)
end subroutine ed_get_gimp_site_n5_c

subroutine ed_get_gimp_lattice_n3_c(self,d,nlat,axis,typ) bind(c, name='ed_get_gimp_lattice_n3')
  integer(c_int64_t)                                                          :: d(3)
  complex(c_double_complex),dimension(d(1),d(2),d(3)),intent(inout)           :: self
  integer(c_int),value                                                        :: nlat
  character(kind=c_char), dimension(1),optional                               :: axis,typ
  character(len=1)                                                            :: axis_,typ_
  typ_(1:1)=typ(1)
  axis_(1:1)=axis(1)
  call ed_get_gimp(self,nlat,axis_,typ_)
end subroutine ed_get_gimp_lattice_n3_c

subroutine ed_get_gimp_lattice_n4_c(self,d,nlat,axis,typ) bind(c, name='ed_get_gimp_lattice_n4')
  integer(c_int64_t)                                                          :: d(4)
  complex(c_double_complex),dimension(d(1),d(2),d(3),d(4)),intent(inout)      :: self
  integer(c_int),value                                                        :: nlat
  character(kind=c_char), dimension(1),optional                               :: axis,typ
  character(len=1)                                                            :: axis_,typ_
  typ_(1:1)=typ(1)
  axis_(1:1)=axis(1)
  call ed_get_gimp(self,nlat,axis_,typ_)
end subroutine ed_get_gimp_lattice_n4_c

subroutine ed_get_gimp_lattice_n6_c(self,nlat,d,axis,typ) bind(c, name='ed_get_gimp_lattice_n6')
  integer(c_int64_t)                                                                    :: d(6)
  complex(c_double_complex),dimension(d(1),d(2),d(3),d(4),d(5),d(6)),intent(inout)      :: self
  integer(c_int),value                                                                  :: nlat
  character(kind=c_char), dimension(1),optional                                         :: axis,typ
  character(len=1)                                                                      :: axis_,typ_
  typ_(1:1)=typ(1)
  axis_(1:1)=axis(1)
  call ed_get_gimp(self,nlat,axis_,typ_)
end subroutine ed_get_gimp_lattice_n6_c


!OBSERVABLES

!density
subroutine ed_get_dens_n1_c(self) bind(c,name="ed_get_dens_n1")
  real(c_double)     :: self(Norb)
  call ed_get_dens(self)
end subroutine ed_get_dens_n1_c


subroutine ed_get_dens_n2(self,d,Nlat) bind(c,name="ed_get_dens_n2")
  integer(c_int64_t) :: d(2)
  real(c_double)     :: self(d(1),d(2))
  integer(c_int)     :: Nlat
  call ed_get_dens(self,Nlat)
end subroutine ed_get_dens_n2