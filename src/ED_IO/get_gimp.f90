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
