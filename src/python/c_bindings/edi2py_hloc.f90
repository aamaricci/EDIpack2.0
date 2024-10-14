!ED_SET_HLOC:

subroutine init_Hreplica_direct_nn_c(Hloc) bind(c, name='init_Hreplica_direct_nn')
  use, intrinsic :: iso_c_binding
  real(c_double),dimension(Nspin,Nspin,Norb,Norb) :: Hloc
  call ed_set_Hreplica(hloc)
end subroutine init_Hreplica_direct_nn_c

subroutine init_Hreplica_direct_so_c(Hloc) bind(c, name='init_Hreplica_direct_so')
  use, intrinsic :: iso_c_binding
  real(c_double),dimension(Nspin*Norb,Nspin*Norb) :: Hloc 
  call ed_set_Hreplica(hloc)
end subroutine init_Hreplica_direct_so_c

subroutine init_Hreplica_symmetries_site_c(Hvec,lambdavec,Nsym) bind(c, name='init_Hreplica_symmetries_site')
  use, intrinsic :: iso_c_binding
  integer(c_int),value                                     :: Nsym
  real(c_double),dimension(Nspin,Nspin,Norb,Norb,Nsym)     :: Hvec
  real(c_double),dimension(Nsym)                           :: lambdavec
  call ed_set_Hreplica(Hvec,lambdavec)
end subroutine init_Hreplica_symmetries_site_C

subroutine init_Hreplica_symmetries_ineq_c(Hvec,lambdavec,Nlat,Nsym) bind(c, name='init_Hreplica_symmetries_ineq')
  use, intrinsic :: iso_c_binding
  integer(c_int),value                                     :: Nlat, Nsym
  real(c_double),dimension(Nspin,Nspin,Norb,Norb,Nsym)     :: Hvec
  real(c_double),dimension(Nlat,Nsym)                      :: lambdavec
  call ed_set_Hreplica(Hvec,lambdavec)
end subroutine init_Hreplica_symmetries_ineq_c
