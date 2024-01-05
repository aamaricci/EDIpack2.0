subroutine ed_get_dens(self,iorb,Nlat)
  real(8),dimension(..) :: self
  integer,optional      :: iorb,Nlat
  integer               :: iorb_
  iorb_=1;if(present(iorb))iorb_=iorb
  if(iorb_>Norb)stop "ed_get_dens error: orbital index > N_orbital"
  if(present(Nlat))then
     if(.not.allocated(dens_ineq))stop "ed_get_dens error: dens_ineq not allocated"
     if(Nlat>size(dens_ineq,1))stop "ed_get_dens error: required N_sites > evaluated N_sites"
  endif  
  select rank(self)
  rank default;stop 'ed_get_dens error: dens has wrong rank'
  rank (0)                      !scalar
  self = ed_dens(iorb_)
  rank (1)                      !all orbitals or all sites 1orbital
  if(present(Nlat))then
     call assert_shape(self,[Nlat],'ed_get_dens','dens')
     self = dens_ineq(:,iorb_)
  else
     call assert_shape(self,[Norb],'ed_get_dens','dens')
     self = ed_dens
  endif
  rank (2)                      !all orbitals and sites
  if(.not.present(Nlat))stop "ed_get_dens error: Nlat require for rank-2 array"
  call assert_shape(self,[Nlat,Norb],'ed_get_dens','dens')
  self = dens_ineq
  end select
end subroutine ed_get_dens



subroutine ed_get_mag(self,component,iorb,Nlat)
  real(8),dimension(..)     :: self
  character(len=1),optional :: component
  integer,optional          :: iorb,Nlat
  !
  integer                   :: iorb_
  character(len=1)          :: char
  integer                   :: id
  !
  iorb_=1;if(present(iorb))iorb_=iorb
  char='Z';if(present(component))char=component
  select case(char)
  case default;id=3
  case('X','x');id=1
  case('Y','y');id=2
  end select
  if(iorb_>Norb)stop "ed_get_mag error: orbital index > N_orbital"
  if(present(Nlat))then
     if(.not.allocated(mag_ineq))stop "ed_get_mag error: mag_ineq not allocated"
     if(Nlat>size(mag_ineq,1))stop "ed_get_mag error: required N_sites > evaluated N_sites"
  endif
  select rank(self)
  rank default;stop 'ed_get_mag error: mag has wrong rank'
  rank (0)                      !scalar
  self = ed_mag(id,iorb_)
  rank (1)                      
  if(present(Nlat))then
     call assert_shape(self,[Nlat],'ed_get_mag','mag')
     self = mag_ineq(:,id,iorb_)
  else
     call assert_shape(self,[Norb],'ed_get_mag','mag')
     self = ed_mag(id,:)
  endif
  rank (2)                      
  if(.not.present(Nlat))stop "ed_get_mag error: Nlat require for rank-2 array"
  call assert_shape(self,[Nlat,Norb],'ed_get_mag','mag')
  self = mag_ineq(:,id,:)
  rank (3)                      
  if(.not.present(Nlat))stop "ed_get_mag error: Nlat require for rank-2 array"
  call assert_shape(self,[Nlat,3,Norb],'ed_get_mag','mag')
  self = mag_ineq
  end select
end subroutine ed_get_mag


subroutine ed_get_docc(self,iorb,Nlat)
  real(8),dimension(..) :: self
  integer,optional      :: iorb,Nlat
  integer               :: iorb_
  iorb_=1;if(present(iorb))iorb_=iorb
  if(iorb_>Norb)stop "ed_get_docc error: orbital index > N_orbital"
  if(present(Nlat))then
     if(.not.allocated(docc_ineq))stop "ed_get_docc error: docc_ineq not allocated"
     if(Nlat>size(docc_ineq,1))stop "ed_get_docc error: required N_sites > evaluated N_sites"
  endif  
  select rank(self)
  rank default;stop 'ed_get_docc error: docc has wrong rank'
  rank (0)                      !scalar
  self = ed_docc(iorb_)
  rank (1)                      !all orbitals or all sites 1orbital
  if(present(Nlat))then
     call assert_shape(self,[Nlat],'ed_get_docc','docc')
     self = docc_ineq(:,iorb_)
  else
     call assert_shape(self,[Norb],'ed_get_docc','docc')
     self = ed_docc
  endif
  rank (2)                      !all orbitals and sites
  if(.not.present(Nlat))stop "ed_get_docc error: Nlat require for rank-2 array"
  call assert_shape(self,[Nlat,Norb],'ed_get_docc','docc')
  self = docc_ineq
  end select
end subroutine ed_get_docc


subroutine ed_get_phisc(self,iorb,Nlat)
  real(8),dimension(..) :: self
  integer,optional      :: iorb,Nlat
  integer               :: iorb_
  iorb_=1;if(present(iorb))iorb_=iorb
  if(iorb_>Norb)stop "ed_get_phisc error: orbital index > N_orbital"
  if(present(Nlat))then
     if(.not.allocated(phisc_ineq))stop "ed_get_phisc error: phisc_ineq not allocated"
     if(Nlat>size(phisc_ineq,1))stop "ed_get_phisc error: required N_sites > evaluated N_sites"
  endif  
  select rank(self)
  rank default;stop 'ed_get_phisc error: phisc has wrong rank'
  rank (0)                      !scalar
  self = ed_phisc(iorb_)
  rank (1)                      !all orbitals or all sites 1orbital
  if(present(Nlat))then
     call assert_shape(self,[Nlat],'ed_get_phisc','phisc')
     self = phisc_ineq(:,iorb_)
  else
     call assert_shape(self,[Norb],'ed_get_phisc','phisc')
     self = ed_phisc
  endif
  rank (2)                      !all orbitals and sites
  if(.not.present(Nlat))stop "ed_get_phisc error: Nlat require for rank-2 array"
  call assert_shape(self,[Nlat,Norb],'ed_get_phisc','phisc')
  self = phisc_ineq
  end select
end subroutine ed_get_phisc


subroutine ed_get_eimp(self,Nlat)
  real(8),dimension(..) :: self
  integer,optional      :: Nlat
  !
  if(present(Nlat))then
     if(.not.allocated(e_ineq))stop "ed_get_eimp error: e_ineq not allocated"
     if(Nlat>size(e_ineq,1))stop "ed_get_eimp error: required N_sites > evaluated N_sites"
  endif
    select rank(self)
  rank default;stop 'ed_get_eimp error: eimp has wrong rank'
  rank (1)                      !all orbitals or all sites 1orbital
  call assert_shape(self,[4],'ed_get_eimp','eimp')
  self = [ed_Epot,ed_Eint,ed_Ehartree,ed_Eknot]
  rank (2)                      !all orbitals and sites
  if(.not.present(Nlat))stop "ed_get_eimp error: Nlat require for rank-2 array"
  call assert_shape(self,[2,Nlat],'ed_get_eimp','eimp')
  self = e_ineq
  end select
end subroutine ed_get_eimp




subroutine ed_get_epot(self,Nlat)
  real(8),dimension(..) :: self
  integer,optional      :: Nlat
  !
  if(present(Nlat))then
     if(.not.allocated(e_ineq))stop "ed_get_epot error: e_ineq not allocated"
     if(Nlat>size(e_ineq,1))stop "ed_get_epot error: required N_sites > evaluated N_sites"
  endif
    select rank(self)
  rank default;stop 'ed_get_epot error: epot has wrong rank'
  rank (0)
  self = ed_Epot
  rank (1)                      !all orbitals and sites
  if(.not.present(Nlat))stop "ed_get_epot error: Nlat require for rank-2 array"
  call assert_shape(self,[Nlat],'ed_get_epot','epot')
  self = e_ineq(:,1)
  end select
end subroutine ed_get_epot


subroutine ed_get_eint(self,Nlat)
  real(8),dimension(..) :: self
  integer,optional      :: Nlat
  !
  if(present(Nlat))then
     if(.not.allocated(e_ineq))stop "ed_get_eint error: e_ineq not allocated"
     if(Nlat>size(e_ineq,1))stop "ed_get_eint error: required N_sites > evaluated N_sites"
  endif
    select rank(self)
  rank default;stop 'ed_get_eint error: eint has wrong rank'
  rank (0)
  self = ed_Eint
  rank (1)                      !all orbitals and sites
  if(.not.present(Nlat))stop "ed_get_eint error: Nlat require for rank-2 array"
  call assert_shape(self,[Nlat],'ed_get_eint','eint')
  self = e_ineq(:,2)
  end select
end subroutine ed_get_eint

subroutine ed_get_ehartree(self,Nlat)
  real(8),dimension(..) :: self
  integer,optional      :: Nlat
  !
  if(present(Nlat))then
     if(.not.allocated(e_ineq))stop "ed_get_ehartree error: e_ineq not allocated"
     if(Nlat>size(e_ineq,1))stop "ed_get_ehartree error: required N_sites > evaluated N_sites"
  endif
    select rank(self)
  rank default;stop 'ed_get_ehartree error: ehartree has wrong rank'
  rank (0)
  self = ed_Ehartree
  rank (1)                      !all orbitals and sites
  if(.not.present(Nlat))stop "ed_get_ehartree error: Nlat require for rank-2 array"
  call assert_shape(self,[Nlat],'ed_get_ehartree','ehartree')
  self = e_ineq(:,3)
  end select
end subroutine ed_get_ehartree


subroutine ed_get_eknot(self,Nlat)
  real(8),dimension(..) :: self
  integer,optional      :: Nlat
  !
  if(present(Nlat))then
     if(.not.allocated(e_ineq))stop "ed_get_eknot error: e_ineq not allocated"
     if(Nlat>size(e_ineq,1))stop "ed_get_eknot error: required N_sites > evaluated N_sites"
  endif
    select rank(self)
  rank default;stop 'ed_get_eknot error: eknot has wrong rank'
  rank (0)
  self = ed_Eknot
  rank (1)                      !all orbitals and sites
  if(.not.present(Nlat))stop "ed_get_eknot error: Nlat require for rank-2 array"
  call assert_shape(self,[Nlat],'ed_get_eknot','eknot')
  self = e_ineq(:,4)
  end select
end subroutine ed_get_eknot










subroutine ed_get_doubles(self,Nlat)
  real(8),dimension(..) :: self
  integer,optional      :: Nlat
  !
  if(present(Nlat))then
     if(.not.allocated(dd_ineq))stop "ed_get_doubles error: dd_ineq not allocated"
     if(Nlat>size(dd_ineq,1))stop "ed_get_doubles error: required N_sites > evaluated N_sites"
  endif
    select rank(self)
  rank default;stop 'ed_get_doubles error: doubles has wrong rank'
  rank (1)                      !all orbitals or all sites 1orbital
  call assert_shape(self,[4],'ed_get_doubles','doubles')
  self = [ed_Dust,ed_Dund,ed_Dse,ed_Dph]
  rank (2)                      !all orbitals and sites
  if(.not.present(Nlat))stop "ed_get_doubles error: Nlat require for rank-2 array"
  call assert_shape(self,[2,Nlat],'ed_get_doubles','doubles')
  self = dd_ineq
  end select
end subroutine ed_get_doubles

subroutine ed_get_dust(self,Nlat)
  real(8),dimension(..) :: self
  integer,optional      :: Nlat
  !
  if(present(Nlat))then
     if(.not.allocated(dd_ineq))stop "ed_get_dust error: dd_ineq not allocated"
     if(Nlat>size(dd_ineq,1))stop "ed_get_dust error: required N_sites > evaluated N_sites"
  endif
    select rank(self)
  rank default;stop 'ed_get_dust error: dust has wrong rank'
  rank (0)
  self = ed_Dust
  rank (1)                      !all orbitals and sites
  if(.not.present(Nlat))stop "ed_get_dust error: Nlat require for rank-2 array"
  call assert_shape(self,[Nlat],'ed_get_dust','dust')
  self = dd_ineq(:,1)
  end select
end subroutine ed_get_dust


subroutine ed_get_dund(self,Nlat)
  real(8),dimension(..) :: self
  integer,optional      :: Nlat
  !
  if(present(Nlat))then
     if(.not.allocated(dd_ineq))stop "ed_get_dund error: dd_ineq not allocated"
     if(Nlat>size(dd_ineq,1))stop "ed_get_dund error: required N_sites > evaluated N_sites"
  endif
    select rank(self)
  rank default;stop 'ed_get_dund error: dund has wrong rank'
  rank (0)
  self = ed_Dund
  rank (1)                      !all orbitals and sites
  if(.not.present(Nlat))stop "ed_get_dund error: Nlat require for rank-2 array"
  call assert_shape(self,[Nlat],'ed_get_dund','dund')
  self = dd_ineq(:,2)
  end select
end subroutine ed_get_dund

subroutine ed_get_dse(self,Nlat)
  real(8),dimension(..) :: self
  integer,optional      :: Nlat
  !
  if(present(Nlat))then
     if(.not.allocated(dd_ineq))stop "ed_get_dse error: dd_ineq not allocated"
     if(Nlat>size(dd_ineq,1))stop "ed_get_dse error: required N_sites > evaluated N_sites"
  endif
    select rank(self)
  rank default;stop 'ed_get_dse error: dse has wrong rank'
  rank (0)
  self = ed_Dse
  rank (1)                      !all orbitals and sites
  if(.not.present(Nlat))stop "ed_get_dse error: Nlat require for rank-2 array"
  call assert_shape(self,[Nlat],'ed_get_dse','dse')
  self = dd_ineq(:,3)
  end select
end subroutine ed_get_dse


subroutine ed_get_dph(self,Nlat)
  real(8),dimension(..) :: self
  integer,optional      :: Nlat
  !
  if(present(Nlat))then
     if(.not.allocated(dd_ineq))stop "ed_get_dph error: dd_ineq not allocated"
     if(Nlat>size(dd_ineq,1))stop "ed_get_dph error: required N_sites > evaluated N_sites"
  endif
    select rank(self)
  rank default;stop 'ed_get_dph error: dph has wrong rank'
  rank (0)
  self = ed_Dph
  rank (1)                      !all orbitals and sites
  if(.not.present(Nlat))stop "ed_get_dph error: Nlat require for rank-2 array"
  call assert_shape(self,[Nlat],'ed_get_dph','dph')
  self = dd_ineq(:,4)
  end select
end subroutine ed_get_dph





subroutine ed_get_neigen_total(nlii,Nlat) 
  integer                      :: Nlat
  integer,dimension(Nlat) :: nlii
  nlii=0d0
  if(allocated(neigen_total_ineq))then
     if(Nlat>size(neigen_total_ineq)) stop "ed_get_docc error: required N_sites > evaluated N_sites"
     nlii=neigen_total_ineq
  endif
end subroutine ed_get_neigen_total


