!NORMAL, MATSUBARA GREEN'S FUNCTIONS
subroutine ed_get_gimp_matsubara_main(Gmats)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Gmats
  Gmats(:,:,:,:,:) = impGmats(:,:,:,:,:)
end subroutine ed_get_gimp_matsubara_main

subroutine ed_get_gimp_matsubara_site(Gmats,ilat)
  integer                                                         :: ilat
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Gmats
  Gmats(:,:,:,:,:) = Gmats_ineq(ilat,:,:,:,:,:)
end subroutine ed_get_gimp_matsubara_site

subroutine ed_get_gimp_matsubara_lattice(Gmats,nlat)
  integer                                                              :: nlat
  complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Gmats
  if(Nlat/=size(Gmats_ineq,1))stop "ERROR ed_get_gimp_matsubara: wrong Nlat"
  Gmats = Gmats_ineq
end subroutine ed_get_gimp_matsubara_lattice



!NORMAL, REAL GREEN'S FUNCTION
subroutine ed_get_gimp_realaxis_main(Greal)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
  Greal(:,:,:,:,:) = impGreal(:,:,:,:,:)
end subroutine ed_get_gimp_realaxis_main

subroutine ed_get_gimp_realaxis_site(Greal,ilat)
  integer                                                         :: ilat
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
  Greal(:,:,:,:,:) = Greal_ineq(ilat,:,:,:,:,:)
end subroutine ed_get_gimp_realaxis_site

subroutine ed_get_gimp_realaxis_lattice(Greal,nlat)
  integer                                                              :: nlat
  complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
  if(Nlat/=size(Greal_ineq,1))stop "ERROR ed_get_gimp_realaxis: wrong Nlat"
  Greal = Greal_ineq
end subroutine ed_get_gimp_realaxis_lattice




!ANOMALous, MATSUBARA GREEN'S FUNCTION
subroutine ed_get_fimp_matsubara_main(Fmats)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Fmats
  Fmats(:,:,:,:,:) = impFmats(:,:,:,:,:)
end subroutine ed_get_fimp_matsubara_main

subroutine ed_get_fimp_matsubara_site(Fmats,ilat)
  integer                                                         :: ilat
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Fmats
  Fmats(:,:,:,:,:) = Fmats_ineq(ilat,:,:,:,:,:)
end subroutine ed_get_fimp_matsubara_site

subroutine ed_get_fimp_matsubara_lattice(Fmats,nlat)
  integer                                                              :: nlat
  complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Fmats
  if(Nlat/=size(Fmats_ineq,1))stop "ERROR ed_get_fimp_matsubara: wrong Nlat"
  Fmats = Fmats_ineq
end subroutine ed_get_fimp_matsubara_lattice







!ANOMALous, REAL GREEN'S FUNCTION
subroutine ed_get_fimp_realaxis_main(Freal)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Freal
  Freal(:,:,:,:,:) = impFreal(:,:,:,:,:)
end subroutine ed_get_fimp_realaxis_main

subroutine ed_get_fimp_realaxis_site(Freal,ilat)
  integer                                                         :: ilat
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Freal
  Freal(:,:,:,:,:) = Freal_ineq(ilat,:,:,:,:,:)
end subroutine ed_get_fimp_realaxis_site

subroutine ed_get_fimp_realaxis_lattice(Freal,nlat)
  integer                                                              :: nlat
  complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Freal
  if(Nlat/=size(Freal_ineq,1))stop "ERROR ed_get_fimp_realaxis: wrong Nlat"
  Freal = Freal_ineq
end subroutine ed_get_fimp_realaxis_lattice
