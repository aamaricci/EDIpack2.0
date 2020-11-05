subroutine ed_get_phisc_1(phisc)
  real(8),dimension(Norb) :: phisc
  phisc = ed_phisc
end subroutine ed_get_phisc_1

subroutine ed_get_phisc_2(phisc,iorb)
  real(8)   :: phisc
  integer   :: iorb
  if(iorb>Norb)stop "ed_get_phisc error: orbital index > N_orbital"
  phisc = ed_phisc(iorb)
end subroutine ed_get_phisc_2

subroutine ed_get_phisc_lattice_1(yii,Nlat) 
  integer                      :: Nlat
  real(8),dimension(Nlat,Norb) :: yii
  yii=0d0
  if(allocated(phisc_ineq))then
     if(Nlat>size(phisc_ineq,1)) stop "ed_get_phisc error: required N_sites > evaluated N_sites"   
     yii=phisc_ineq
  endif
end subroutine ed_get_phisc_lattice_1


subroutine ed_get_phisc_lattice_2(yii,Nlat,iorb) 
  integer                 :: Nlat
  real(8),dimension(Nlat) :: yii
  integer                 :: iorb
  if(iorb>Norb)stop "ed_get_phisc error: orbital index > N_orbital"
  yii=0d0
  if(allocated(phisc_ineq))then
     if(Nlat>size(phisc_ineq,1)) stop "ed_get_phisc error: required N_sites > evaluated N_sites"
     yii=phisc_ineq(:,iorb)
  endif
end subroutine ed_get_phisc_lattice_2
