subroutine ed_get_doubles_(docc)
  real(8),dimension(4) :: docc
  docc = [ed_Dust,ed_Dund,ed_Dse,ed_Dph]
end subroutine ed_get_doubles_

subroutine ed_get_dust_(docc)
  real(8) :: docc
  docc = ed_Dust
end subroutine ed_get_dust_

subroutine ed_get_dund_(docc)
  real(8) :: docc
  docc = ed_Dund
end subroutine ed_get_dund_

subroutine ed_get_dse_(docc)
  real(8) :: docc
  docc = ed_Dse
end subroutine ed_get_dse_

subroutine ed_get_dph_(docc)
  real(8) :: docc
  docc = ed_Dph
end subroutine ed_get_dph_

subroutine ed_get_doubles_lattice(yii,Nlat)
  integer                      :: Nlat
  real(8),dimension(Nlat,4)    :: yii
  yii=0d0
  if(allocated(dd_ineq))then
     if(Nlat>size(dd_ineq,1)) stop "ed_get_doubles error: required N_sites > evaluated N_sites"
     yii=dd_ineq(:,:)
  endif
end subroutine ed_get_doubles_lattice

subroutine ed_get_dust_lattice(yii,Nlat)
  integer                 :: Nlat
  real(8),dimension(Nlat) :: yii
  yii=0d0
  if(allocated(dd_ineq))then
     if(Nlat>size(dd_ineq,1)) stop "ed_get_dust error: required N_sites > evaluated N_sites"
     yii=dd_ineq(:,1)
  endif
end subroutine ed_get_dust_lattice

subroutine ed_get_dund_lattice(yii,Nlat)
  integer                 :: Nlat
  real(8),dimension(Nlat) :: yii
  yii=0d0
  if(allocated(dd_ineq))then
     if(Nlat>size(dd_ineq,1)) stop "ed_get_dund error: required N_sites > evaluated N_sites"
     yii=dd_ineq(:,2)
  endif
end subroutine ed_get_dund_lattice

subroutine ed_get_dse_lattice(yii,Nlat)
  integer                 :: Nlat
  real(8),dimension(Nlat) :: yii
  yii=0d0
  if(allocated(dd_ineq))then
     if(Nlat>size(dd_ineq,1)) stop "ed_get_dse error: required N_sites > evaluated N_sites"
     yii=dd_ineq(:,3)
  endif
end subroutine ed_get_dse_lattice

subroutine ed_get_dph_lattice(yii,Nlat)
  integer                 :: Nlat
  real(8),dimension(Nlat) :: yii
  yii=0d0
  if(allocated(dd_ineq))then
     if(Nlat>size(dd_ineq,1)) stop "ed_get_dph error: required N_sites > evaluated N_sites"
     yii=dd_ineq(:,4)
  endif
end subroutine ed_get_dph_lattice
