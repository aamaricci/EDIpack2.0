subroutine ed_get_phisc_n0(self,iorb,jorb)
  real(8)          :: self ! :math:`\phi` value or array of values
  integer,optional :: iorb ! first orbital index
  integer,optional :: jorb ! second orbital index
  integer          :: iorb_,jorb_
  iorb_=1;if(present(iorb))iorb_=iorb
  jorb_=1;if(present(jorb))jorb_=jorb
  if(iorb_>Norb.OR.jorb_>Norb)stop "ed_get_phisc error: orbital index > N_orbital"
  self = ed_phisc(iorb_,jorb_)
end subroutine ed_get_phisc_n0

!phi(a,a) or phi_i(a_,a_)
subroutine ed_get_phisc_n1(self,iorb,Nlat)
  real(8),dimension(:) :: self
  integer,optional     :: iorb
  integer,optional     :: Nlat ! number of inequivalent impurity sites for real-space DMFT
  integer              :: iorb_
  iorb_=1;if(present(iorb))iorb_=iorb
  if(iorb_>Norb)stop "ed_get_phisc error: orbital index > N_orbital"
  if(present(Nlat))then
     if(.not.allocated(phisc_ineq))stop "ed_get_phisc error: phisc_ineq not allocated"
     if(Nlat>size(phisc_ineq,1))stop "ed_get_phisc error: required N_sites > evaluated N_sites"
  endif
  if(present(Nlat))then
     call assert_shape(self,[Nlat],'ed_get_phisc','phisc')
     self = phisc_ineq(:,iorb_,iorb_)
  else
     call assert_shape(self,[Norb],'ed_get_phisc','phisc')
     self = diagonal(ed_phisc)
  endif
end subroutine ed_get_phisc_n1


!phi_ab
subroutine ed_get_phisc_n2(self)
  real(8),dimension(:,:) :: self
  call assert_shape(self,[Norb,Norb],'ed_get_phisc','phisc')
  self = ed_phisc
end subroutine ed_get_phisc_n2


!phi_i,ab
subroutine ed_get_phisc_n3(self,Nlat)
  real(8),dimension(:,:,:) :: self
  integer                  :: Nlat
  if(.not.allocated(phisc_ineq))stop "ed_get_phisc error: phisc_ineq not allocated"
  if(Nlat>size(phisc_ineq,1))stop "ed_get_phisc error: required N_sites > evaluated N_sites"
  call assert_shape(self,[Nlat,Norb,Norb],'ed_get_phisc','phisc')
  self = phisc_ineq
end subroutine ed_get_phisc_n3

