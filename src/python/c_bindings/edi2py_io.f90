!ED_IO:

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
  integer(c_int)           :: iorb
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

!---------!
!GET SIGMA!
!---------!

!SITE
subroutine get_sigma_site_n3_c(self,axis,typ,zeta,dz,zflag) bind(c,name="get_sigma_site_n3")
  use, intrinsic :: iso_c_binding
  integer(c_int),value                          :: dz,axis,typ,zflag
  character(len=1)                              :: axis_
  character(len=1)                              :: typ_
  complex(c_double_complex)                     :: zeta(dz)
  complex(c_double_complex)                     :: self(Nspin*Norb,Nspin*Norb,dz)
  !
  axis_="m"
  if(axis==1)axis_="r"
  typ_="n"
  if(typ==1)typ_="a"
  if(zflag==1)then
    call ed_get_sigma(self,axis_,typ_,zeta)
  else
    call ed_get_sigma(self,axis_,typ_)
  endif
end subroutine get_sigma_site_n3_c

subroutine get_sigma_site_n5_c(self,axis,typ,zeta,dz,zflag) bind(c,name="get_sigma_site_n5")
  use, intrinsic :: iso_c_binding
  integer(c_int),value                          :: dz,axis,typ,zflag
  character(len=1)                              :: axis_
  character(len=1)                              :: typ_
  complex(c_double_complex)                     :: zeta(dz)
  complex(c_double_complex)                     :: self(Nspin,Nspin,Norb,Norb,dz)
  !
  axis_="m"
  if(axis==1)axis_="r"
  typ_="n"
  if(typ==1)typ_="a"
  if(zflag==1)then
    call ed_get_sigma(self,axis_,typ_,zeta)
  else
    call ed_get_sigma(self,axis_,typ_)
  endif
end subroutine get_sigma_site_n5_c

!LATTICE
subroutine get_sigma_lattice_n3_c(self,Nineq,axis,typ,zeta,dz,zflag) bind(c,name="get_sigma_lattice_n3")
  use, intrinsic :: iso_c_binding
  integer(c_int),value                          :: dz,axis,Nineq,typ,zflag
  character(len=1)                              :: axis_
  character(len=1)                              :: typ_
  complex(c_double_complex)                     :: zeta(dz)
  complex(c_double_complex)                     :: self(Nineq*Nspin*Norb,Nineq*Nspin*Norb,dz)
  !
  axis_="m"
  if(axis==1)axis_="r"
  typ_="n"
  if(typ==1)typ_="a"
  if(zflag==1)then
    call ed_get_sigma(self,Nineq,axis_,typ_,zeta)
  else
    call ed_get_sigma(self,Nineq,axis_,typ_)
  endif
end subroutine get_sigma_lattice_n3_c

subroutine get_sigma_lattice_n4_c(self,Nineq,axis,typ,zeta,dz,zflag) bind(c,name="get_sigma_lattice_n4")
  use, intrinsic :: iso_c_binding
  integer(c_int),value                          :: dz,axis,Nineq,typ,zflag
  character(len=1)                              :: axis_
  character(len=1)                              :: typ_
  complex(c_double_complex)                     :: zeta(dz)
  complex(c_double_complex)                     :: self(Nineq,Nspin*Norb,Nspin*Norb,dz)
  !
  axis_="m"
  if(axis==1)axis_="r"
  typ_="n"
  if(typ==1)typ_="a"
  if(zflag==1)then
    call ed_get_sigma(self,Nineq,axis_,typ_,zeta)
  else
    call ed_get_sigma(self,Nineq,axis_,typ_)
  endif
end subroutine get_sigma_lattice_n4_c

subroutine get_sigma_lattice_n6_c(self,Nineq,axis,typ,zeta,dz,zflag) bind(c,name="get_sigma_lattice_n6")
  use, intrinsic :: iso_c_binding
  integer(c_int),value                          :: dz,axis,Nineq,typ,zflag
  character(len=1)                              :: axis_
  character(len=1)                              :: typ_
  complex(c_double_complex)                     :: zeta(dz)
  complex(c_double_complex)                     :: self(Nineq,Nspin,Nspin,Norb,Norb,dz)
  !
  axis_="m"
  if(axis==1)axis_="r"
  typ_="n"
  if(typ==1)typ_="a"
  if(zflag==1)then
    call ed_get_sigma(self,Nineq,axis_,typ_,zeta)
  else
    call ed_get_sigma(self,Nineq,axis_,typ_)
  endif
end subroutine get_sigma_lattice_n6_c

!---------!
!GET GIMP !
!---------!

!SITE
subroutine get_gimp_site_n3_c(gimp,axis,typ,zeta,dz,zflag) bind(c,name="get_gimp_site_n3")
  use, intrinsic :: iso_c_binding
  integer(c_int),value                          :: dz,axis,typ,zflag
  character(len=1)                              :: axis_
  character(len=1)                              :: typ_
  complex(c_double_complex)                     :: zeta(dz)
  complex(c_double_complex)                     :: gimp(Nspin*Norb,Nspin*Norb,dz)
  !
  axis_="m"
  if(axis==1)axis_="r"
  typ_="n"
  if(typ==1)typ_="a"
  if(zflag==1)then
    call ed_get_gimp(gimp,axis_,typ_,zeta)
  else
    call ed_get_gimp(gimp,axis_,typ_)
  endif
end subroutine get_gimp_site_n3_c

subroutine get_gimp_site_n5_c(gimp,axis,typ,zeta,dz,zflag) bind(c,name="get_gimp_site_n5")
  use, intrinsic :: iso_c_binding
  integer(c_int),value                          :: dz,axis,typ,zflag
  character(len=1)                              :: axis_
  character(len=1)                              :: typ_
  complex(c_double_complex)                     :: zeta(dz)
  complex(c_double_complex)                     :: gimp(Nspin,Nspin,Norb,Norb,dz)
  !
  axis_="m"
  if(axis==1)axis_="r"
  typ_="n"
  if(typ==1)typ_="a"
  if(zflag==1)then
    call ed_get_gimp(gimp,axis_,typ_,zeta)
  else
    call ed_get_gimp(gimp,axis_,typ_)
  endif
end subroutine get_gimp_site_n5_c

!LATTICE
subroutine get_gimp_lattice_n3_c(gimp,Nineq,axis,typ,zeta,dz,zflag) bind(c,name="get_gimp_lattice_n3")
  use, intrinsic :: iso_c_binding
  integer(c_int),value                          :: dz,axis,Nineq,typ,zflag
  character(len=1)                              :: axis_
  character(len=1)                              :: typ_
  complex(c_double_complex)                     :: zeta(dz)
  complex(c_double_complex)                     :: gimp(Nineq*Nspin*Norb,Nineq*Nspin*Norb,dz)
  !
  axis_="m"
  if(axis==1)axis_="r"
  typ_="n"
  if(typ==1)typ_="a"
  if(zflag==1)then
    call ed_get_gimp(gimp,Nineq,axis_,typ_,zeta)
  else
    call ed_get_gimp(gimp,Nineq,axis_,typ_)
  endif
end subroutine get_gimp_lattice_n3_c

subroutine get_gimp_lattice_n4_c(gimp,Nineq,axis,typ,zeta,dz,zflag) bind(c,name="get_gimp_lattice_n4")
  use, intrinsic :: iso_c_binding
  integer(c_int),value                          :: dz,axis,Nineq,typ,zflag
  character(len=1)                              :: axis_
  character(len=1)                              :: typ_
  complex(c_double_complex)                     :: zeta(dz)
  complex(c_double_complex)                     :: gimp(Nineq,Nspin*Norb,Nspin*Norb,dz)
  !
  axis_="m"
  if(axis==1)axis_="r"
  typ_="n"
  if(typ==1)typ_="a"
  if(zflag==1)then
    call ed_get_gimp(gimp,Nineq,axis_,typ_,zeta)
  else
    call ed_get_gimp(gimp,Nineq,axis_,typ_)
  endif
end subroutine get_gimp_lattice_n4_c

subroutine get_gimp_lattice_n6_c(gimp,Nineq,axis,typ,zeta,dz,zflag) bind(c,name="get_gimp_lattice_n6")
  use, intrinsic :: iso_c_binding
  integer(c_int),value                          :: dz,axis,Nineq,typ,zflag
  character(len=1)                              :: axis_
  character(len=1)                              :: typ_
  complex(c_double_complex)                     :: zeta(dz)
  complex(c_double_complex)                     :: gimp(Nineq,Nspin,Nspin,Norb,Norb,dz)
  !
  axis_="m"
  if(axis==1)axis_="r"
  typ_="n"
  if(typ==1)typ_="a"
  if(zflag==1)then
    call ed_get_gimp(gimp,Nineq,axis_,typ_,zeta)
  else
    call ed_get_gimp(gimp,Nineq,axis_,typ_)
  endif
end subroutine get_gimp_lattice_n6_c
