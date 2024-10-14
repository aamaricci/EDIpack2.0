!ED_IO:

!Sigma
subroutine ed_get_sigma_site_n3_c(self,d,axis,typ) bind(c, name='ed_get_sigma_site_n3')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                                :: d(3)
  complex(c_double_complex),dimension(d(1),d(2),d(3)),intent(inout) :: self
  character(kind=c_char), dimension(1),optional                     :: axis,typ
  character(len=1)                                                  :: axis_,typ_
  typ_(1:1)=typ(1)
  axis_(1:1)=axis(1)
  call ed_get_sigma(self,axis_,typ_)
end subroutine ed_get_sigma_site_n3_c

subroutine ed_get_sigma_site_n5_c(self,d,axis,typ) bind(c, name='ed_get_sigma_site_n5')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                                          :: d(5)
  complex(c_double_complex),dimension(d(1),d(2),d(3),d(4),d(5)),intent(inout) :: self
  character(kind=c_char), dimension(1),optional                               :: axis,typ
  character(len=1)                                                            :: axis_,typ_
  typ_(1:1)=typ(1)
  axis_(1:1)=axis(1)
  call ed_get_sigma(self,axis_,typ_)
end subroutine ed_get_sigma_site_n5_c

subroutine ed_get_sigma_lattice_n3_c(self,d,nlat,axis,typ) bind(c, name='ed_get_sigma_lattice_n3')
  use, intrinsic :: iso_c_binding
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
  use, intrinsic :: iso_c_binding
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
  use, intrinsic :: iso_c_binding
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
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                                :: d(3)
  complex(c_double_complex),dimension(d(1),d(2),d(3)),intent(inout) :: self
  character(kind=c_char), dimension(1),optional                     :: axis,typ
  character(len=1)                                                  :: axis_,typ_
  typ_(1:1)=typ(1)
  axis_(1:1)=axis(1)
  call ed_get_gimp(self,axis_,typ_)
end subroutine ed_get_gimp_site_n3_c

subroutine ed_get_gimp_site_n5_c(self,d,axis,typ) bind(c, name='ed_get_gimp_site_n5')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                                          :: d(5)
  complex(c_double_complex),dimension(d(1),d(2),d(3),d(4),d(5)),intent(inout) :: self
  character(kind=c_char), dimension(1),optional                               :: axis,typ
  character(len=1)                                                            :: axis_,typ_
  typ_(1:1)=typ(1)
  axis_(1:1)=axis(1)
  call ed_get_gimp(self,axis_,typ_)
end subroutine ed_get_gimp_site_n5_c

subroutine ed_get_gimp_lattice_n3_c(self,d,nlat,axis,typ) bind(c, name='ed_get_gimp_lattice_n3')
  use, intrinsic :: iso_c_binding
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
  use, intrinsic :: iso_c_binding
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
  use, intrinsic :: iso_c_binding
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
  use, intrinsic :: iso_c_binding
  real(c_double)     :: self(Norb)
  call ed_get_dens(self)
end subroutine ed_get_dens_n1_c

subroutine ed_get_dens_n2_c(self,Nlat) bind(c,name="ed_get_dens_n2")
  use, intrinsic :: iso_c_binding
  real(c_double)           :: self(Nlat,Norb)
  integer(c_int),value     :: Nlat
  call ed_get_dens(self,Nlat)
end subroutine ed_get_dens_n2_c

!magnetization
subroutine ed_get_mag_n2_c(self) bind(c,name="ed_get_mag_n2")
  use, intrinsic :: iso_c_binding
  real(c_double)           :: self(3,Norb)
  integer(c_int)           :: icomp,iorb
  do iorb = 1,Norb
    call ed_get_mag(self(1,iorb),component="x",iorb=iorb)
    call ed_get_mag(self(2,iorb),component="y",iorb=iorb)
    call ed_get_mag(self(3,iorb),component="z",iorb=iorb)
  enddo
end subroutine ed_get_mag_n2_c

subroutine ed_get_mag_n3_c(self,Nlat) bind(c,name="ed_get_mag_n3")
  use, intrinsic :: iso_c_binding
  real(c_double)           :: self(Nlat,3,Norb)
  integer(c_int),value     :: Nlat
  call ed_get_mag(self(:,1,:),"x",Nlat)
  call ed_get_mag(self(:,2,:),"y",Nlat)
  call ed_get_mag(self(:,3,:),"z",Nlat)
end subroutine ed_get_mag_n3_c

!double occupation
subroutine ed_get_docc_n1_c(self) bind(c,name="ed_get_docc_n1")
  use, intrinsic :: iso_c_binding
  real(c_double)     :: self(Norb)
  call ed_get_docc(self)
end subroutine ed_get_docc_n1_c

subroutine ed_get_docc_n2_c(self,Nlat) bind(c,name="ed_get_docc_n2")
  use, intrinsic :: iso_c_binding
  real(c_double)           :: self(Nlat,Norb)
  integer(c_int),value     :: Nlat
  call ed_get_docc(self,Nlat)
end subroutine ed_get_docc_n2_c

!energy
subroutine ed_get_eimp_n1_c(self) bind(c,name="ed_get_eimp_n1")
  use, intrinsic :: iso_c_binding
  real(c_double) :: self(4)
  call ed_get_eimp(self)
end subroutine ed_get_eimp_n1_c

subroutine ed_get_eimp_n2_c(self,Nlat) bind(c,name="ed_get_eimp_n2")
  use, intrinsic :: iso_c_binding
  real(c_double)                   :: self(Nlat,4)
  integer(c_int),value             :: Nlat
  call ed_get_eimp(self,Nlat)
end subroutine ed_get_eimp_n2_c

!rebuild sigma
subroutine rebuild_sigma_single_n3_c(zeta,dz,sigma_normal,sigma_anomalous) bind(c,name="build_sigma_single_n3")
  use, intrinsic :: iso_c_binding
  complex(c_double_complex)                     :: zeta(dz)
  complex(c_double_complex)                     :: sigma_normal(Nspin*Norb,Nspin*Norb,dz)
  complex(c_double_complex)                     :: sigma_anomalous(Nspin*Norb,Nspin*Norb,dz)
  integer(c_int),value                          :: dz
  call ed_build_sigma(zeta,sigma_normal,sigma_anomalous)
end subroutine rebuild_sigma_single_n3_c

subroutine rebuild_sigma_single_n5_c(zeta,dz,sigma_normal,sigma_anomalous) bind(c,name="build_sigma_single_n5")
  use, intrinsic :: iso_c_binding
  complex(c_double_complex)                     :: zeta(dz)
  complex(c_double_complex)                     :: sigma_normal(Nspin,Nspin,Norb,Norb,dz)
  complex(c_double_complex)                     :: sigma_anomalous(Nspin,Nspin,Norb,Norb,dz)
  integer(c_int),value                          :: dz
  call ed_build_sigma(zeta,sigma_normal,sigma_anomalous)
end subroutine rebuild_sigma_single_n5_c

subroutine rebuild_sigma_ineq_n3_c(zeta,dz,Nineq,sigma_normal,sigma_anomalous) bind(c,name="build_sigma_ineq_n3")
  use, intrinsic :: iso_c_binding
  complex(c_double_complex)                     :: zeta(dz)
  complex(c_double_complex)                     :: sigma_normal(Nineq*Nspin*Norb,Nineq*Nspin*Norb,dz)
  complex(c_double_complex)                     :: sigma_anomalous(Nineq*Nspin*Norb,Nineq*Nspin*Norb,dz)
  integer(c_int),value                          :: dz,Nineq
  call ed_build_sigma(zeta,Nineq,sigma_normal,sigma_anomalous)
end subroutine rebuild_sigma_ineq_n3_c

subroutine rebuild_sigma_ineq_n4_c(zeta,dz,Nineq,sigma_normal,sigma_anomalous) bind(c,name="build_sigma_ineq_n4")
  use, intrinsic :: iso_c_binding
  complex(c_double_complex)                     :: zeta(dz)
  complex(c_double_complex)                     :: sigma_normal(Nineq,Nspin*Norb,Nspin*Norb,dz)
  complex(c_double_complex)                     :: sigma_anomalous(Nineq,Nspin*Norb,Nspin*Norb,dz)
  integer(c_int),value                          :: dz,Nineq
  call ed_build_sigma(zeta,Nineq,sigma_normal,sigma_anomalous)
end subroutine rebuild_sigma_ineq_n4_c

subroutine rebuild_sigma_ineq_n6_c(zeta,dz,Nineq,sigma_normal,sigma_anomalous) bind(c,name="build_sigma_ineq_n6")
  use, intrinsic :: iso_c_binding
  complex(c_double_complex)                     :: zeta(dz)
  complex(c_double_complex)                     :: sigma_normal(Nineq,Nspin,Nspin,Norb,Norb,dz)
  complex(c_double_complex)                     :: sigma_anomalous(Nineq,Nspin,Nspin,Norb,Norb,dz)
  integer(c_int),value                          :: dz,Nineq
  call ed_build_sigma(zeta,Nineq,sigma_normal,sigma_anomalous)
end subroutine rebuild_sigma_ineq_n6_c

!rebuild gimp
subroutine rebuild_gimp_single_n3_c(zeta,dz,gimp_normal,gimp_anomalous) bind(c,name="build_gimp_single_n3")
  use, intrinsic :: iso_c_binding
  complex(c_double_complex)                     :: zeta(dz)
  complex(c_double_complex)                     :: gimp_normal(Nspin*Norb,Nspin*Norb,dz)
  complex(c_double_complex)                     :: gimp_anomalous(Nspin*Norb,Nspin*Norb,dz)
  integer(c_int),value                          :: dz
  call ed_build_gimp(zeta,gimp_normal,gimp_anomalous)
end subroutine rebuild_gimp_single_n3_c

subroutine rebuild_gimp_single_n5_c(zeta,dz,gimp_normal,gimp_anomalous) bind(c,name="build_gimp_single_n5")
  use, intrinsic :: iso_c_binding
  complex(c_double_complex)                     :: zeta(dz)
  complex(c_double_complex)                     :: gimp_normal(Nspin,Nspin,Norb,Norb,dz)
  complex(c_double_complex)                     :: gimp_anomalous(Nspin,Nspin,Norb,Norb,dz)
  integer(c_int),value                          :: dz
  call ed_build_gimp(zeta,gimp_normal,gimp_anomalous)
end subroutine rebuild_gimp_single_n5_c

subroutine rebuild_gimp_ineq_n3_c(zeta,dz,Nineq,gimp_normal,gimp_anomalous) bind(c,name="build_gimp_ineq_n3")
  use, intrinsic :: iso_c_binding
  complex(c_double_complex)                     :: zeta(dz)
  complex(c_double_complex)                     :: gimp_normal(Nineq*Nspin*Norb,Nineq*Nspin*Norb,dz)
  complex(c_double_complex)                     :: gimp_anomalous(Nineq*Nspin*Norb,Nineq*Nspin*Norb,dz)
  integer(c_int),value                          :: dz,Nineq
  call ed_build_gimp(zeta,Nineq,gimp_normal,gimp_anomalous)
end subroutine rebuild_gimp_ineq_n3_c

subroutine rebuild_gimp_ineq_n4_c(zeta,dz,Nineq,gimp_normal,gimp_anomalous) bind(c,name="build_gimp_ineq_n4")
  use, intrinsic :: iso_c_binding
  complex(c_double_complex)                     :: zeta(dz)
  complex(c_double_complex)                     :: gimp_normal(Nineq,Nspin*Norb,Nspin*Norb,dz)
  complex(c_double_complex)                     :: gimp_anomalous(Nineq,Nspin*Norb,Nspin*Norb,dz)
  integer(c_int),value                          :: dz,Nineq
  call ed_build_gimp(zeta,Nineq,gimp_normal,gimp_anomalous)
end subroutine rebuild_gimp_ineq_n4_c

subroutine rebuild_gimp_ineq_n6_c(zeta,dz,Nineq,gimp_normal,gimp_anomalous) bind(c,name="build_gimp_ineq_n6")
  use, intrinsic :: iso_c_binding
  complex(c_double_complex)                     :: zeta(dz)
  complex(c_double_complex)                     :: gimp_normal(Nineq,Nspin,Nspin,Norb,Norb,dz)
  complex(c_double_complex)                     :: gimp_anomalous(Nineq,Nspin,Nspin,Norb,Norb,dz)
  integer(c_int),value                          :: dz,Nineq
  call ed_build_gimp(zeta,Nineq,gimp_normal,gimp_anomalous)
end subroutine rebuild_gimp_ineq_n6_c
