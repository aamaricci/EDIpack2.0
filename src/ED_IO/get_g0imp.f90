!NORMAL, MATSUBARA GREEN'S FUNCTIONS
subroutine ed_get_g0imp_matsubara_main(Gmats,bath)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Gmats
  real(8),dimension(:)                                            :: bath
  call allocate_grids
  Gmats(:,:,:,:,:) = ed_get_g0and_function(xi*wm,bath)
  call deallocate_grids
end subroutine ed_get_g0imp_matsubara_main


!NORMAL, REAL GREEN'S FUNCTION
subroutine ed_get_g0imp_realaxis_main(Greal,bath)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
  real(8),dimension(:)                                            :: bath
  call allocate_grids
  Greal(:,:,:,:,:) = ed_get_g0and_function(dcmplx(wr,eps),bath,axis='real')
  call deallocate_grids
end subroutine ed_get_g0imp_realaxis_main


!ANOMALous, MATSUBARA GREEN'S FUNCTION
subroutine ed_get_f0imp_matsubara_main(Fmats,bath)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Fmats
  real(8),dimension(:)                                            :: bath
  call allocate_grids
  Fmats(:,:,:,:,:) = ed_get_f0and_function(xi*wm,bath)
  call deallocate_grids
end subroutine ed_get_f0imp_matsubara_main


!ANOMALous, REAL GREEN'S FUNCTION
subroutine ed_get_f0imp_realaxis_main(Freal,bath)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Freal
  real(8),dimension(:)                                            :: bath
  call allocate_grids
  Freal(:,:,:,:,:) = ed_get_f0and_function(dcmplx(wr,eps),bath,axis='real')
  call deallocate_grids
end subroutine ed_get_f0imp_realaxis_main
