!NORMAL, MATSUBARA SELF-ENEGRGY
subroutine ed_get_sigma_matsubara_main(Smats)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Smats
  Smats(:,:,:,:,:) = impSmats(:,:,:,:,:)
end subroutine ed_get_sigma_matsubara_main

subroutine ed_get_sigma_matsubara_site(Smats,ilat)
  integer                                                         :: ilat
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Smats
  Smats(:,:,:,:,:) = Smats_ineq(ilat,:,:,:,:,:)
end subroutine ed_get_sigma_matsubara_site


!NORMAL, REAL SELF-ENERGY
subroutine ed_get_sigma_realaxis_main(Sreal)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Sreal
  Sreal(:,:,:,:,:) = impSreal(:,:,:,:,:)
end subroutine ed_get_sigma_realaxis_main

subroutine ed_get_sigma_realaxis_site(Sreal,ilat)
  integer                                                         :: ilat
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Sreal
  Sreal(:,:,:,:,:) = Sreal_ineq(ilat,:,:,:,:,:)
end subroutine ed_get_sigma_realaxis_site



!ANOMALous, MATSUBARA SELF-ENERGY
subroutine ed_get_self_matsubara_main(SAmats)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: SAmats
  SAmats(:,:,:,:,:) = impSAmats(:,:,:,:,:)
end subroutine ed_get_self_matsubara_main

subroutine ed_get_self_matsubara_site(SAmats,ilat)
  integer                                                         :: ilat
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: SAmats
  SAmats(:,:,:,:,:) = SAmats_ineq(ilat,:,:,:,:,:)
end subroutine ed_get_self_matsubara_site


!ANOMALous, REAL SELF-ENERGY
subroutine ed_get_self_realaxis_main(SAreal)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: SAreal
  SAreal(:,:,:,:,:) = impSAreal(:,:,:,:,:)
end subroutine ed_get_self_realaxis_main

subroutine ed_get_self_realaxis_site(SAreal,ilat)
  integer                                                         :: ilat
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: SAreal
  SAreal(:,:,:,:,:) = SAreal_ineq(ilat,:,:,:,:,:)
end subroutine ed_get_self_realaxis_site

