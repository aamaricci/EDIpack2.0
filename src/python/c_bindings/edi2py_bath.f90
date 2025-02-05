!ED_BATH:
integer(c_int) function get_bath_dimension_c() result(Nb) bind(c, name='get_bath_dimension')
  use, intrinsic :: iso_c_binding
  Nb=ed_get_bath_dimension()
end function get_bath_dimension_c

!H_REPLICA SETUP
subroutine init_Hreplica_symmetries_d5_c(Hvec,d_hvec,lambdavec,d_lambdavec) bind(c, name='init_Hreplica_symmetries_d5')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                                                     :: d_hvec(5), d_lambdavec(2)
  complex(c_double_complex),dimension(d_hvec(1),d_hvec(2),d_hvec(3),d_hvec(4),d_hvec(5)) :: Hvec
  real(c_double),dimension(d_lambdavec(1),d_lambdavec(2))                                :: lambdavec
  call ed_set_Hreplica(Hvec,lambdavec)
end subroutine init_Hreplica_symmetries_d5_c

subroutine init_Hreplica_symmetries_d3_c(Hvec,d_hvec,lambdavec,d_lambdavec) bind(c, name='init_Hreplica_symmetries_d3')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                                 :: d_hvec(3), d_lambdavec(2)
  complex(c_double_complex),dimension(d_hvec(1),d_hvec(2),d_hvec(3)) :: Hvec
  real(c_double),dimension(d_lambdavec(1),d_lambdavec(2))            :: lambdavec
  call ed_set_Hreplica(Hvec,lambdavec)
end subroutine init_Hreplica_symmetries_d3_c

subroutine init_Hreplica_symmetries_lattice_d5_c(Hvec,d_hvec,lambdavec,d_lambdavec) bind(c, name='init_Hreplica_symmetries_lattice_d5')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                                                     :: d_hvec(5), d_lambdavec(3)
  complex(c_double_complex),dimension(d_hvec(1),d_hvec(2),d_hvec(3),d_hvec(4),d_hvec(5)) :: Hvec
  real(c_double),dimension(d_lambdavec(1),d_lambdavec(2),d_lambdavec(3))                 :: lambdavec
  call ed_set_Hreplica(Hvec,lambdavec)
end subroutine init_Hreplica_symmetries_lattice_d5_c

subroutine init_Hreplica_symmetries_lattice_d3_c(Hvec,d_hvec,lambdavec,d_lambdavec) bind(c, name='init_Hreplica_symmetries_lattice_d3')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                                                     :: d_hvec(3), d_lambdavec(3)
  complex(c_double_complex),dimension(d_hvec(1),d_hvec(2),d_hvec(3))                     :: Hvec
  real(c_double),dimension(d_lambdavec(1),d_lambdavec(2),d_lambdavec(3))                 :: lambdavec
  call ed_set_Hreplica(Hvec,lambdavec)
end subroutine init_Hreplica_symmetries_lattice_d3_c


!H_general SETUP
subroutine init_Hgeneral_symmetries_d5_c(Hvec,d_hvec,lambdavec,d_lambdavec) bind(c, name='init_Hgeneral_symmetries_d5')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                                                     :: d_hvec(5), d_lambdavec(2)
  complex(c_double_complex),dimension(d_hvec(1),d_hvec(2),d_hvec(3),d_hvec(4),d_hvec(5)) :: Hvec
  real(c_double),dimension(d_lambdavec(1),d_lambdavec(2))                                :: lambdavec
  call ed_set_Hgeneral(Hvec,lambdavec)
end subroutine init_Hgeneral_symmetries_d5_c

subroutine init_Hgeneral_symmetries_d3_c(Hvec,d_hvec,lambdavec,d_lambdavec) bind(c, name='init_Hgeneral_symmetries_d3')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                                 :: d_hvec(3), d_lambdavec(2)
  complex(c_double_complex),dimension(d_hvec(1),d_hvec(2),d_hvec(3)) :: Hvec
  real(c_double),dimension(d_lambdavec(1),d_lambdavec(2))            :: lambdavec
  call ed_set_Hgeneral(Hvec,lambdavec)
end subroutine init_Hgeneral_symmetries_d3_c

subroutine init_Hgeneral_symmetries_lattice_d5_c(Hvec,d_hvec,lambdavec,d_lambdavec) bind(c, name='init_Hgeneral_symmetries_lattice_d5')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                                                     :: d_hvec(5), d_lambdavec(3)
  complex(c_double_complex),dimension(d_hvec(1),d_hvec(2),d_hvec(3),d_hvec(4),d_hvec(5)) :: Hvec
  real(c_double),dimension(d_lambdavec(1),d_lambdavec(2),d_lambdavec(3))                 :: lambdavec
  call ed_set_Hgeneral(Hvec,lambdavec)
end subroutine init_Hgeneral_symmetries_lattice_d5_c

subroutine init_Hgeneral_symmetries_lattice_d3_c(Hvec,d_hvec,lambdavec,d_lambdavec) bind(c, name='init_Hgeneral_symmetries_lattice_d3')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                                                     :: d_hvec(3), d_lambdavec(3)
  complex(c_double_complex),dimension(d_hvec(1),d_hvec(2),d_hvec(3))                     :: Hvec
  real(c_double),dimension(d_lambdavec(1),d_lambdavec(2),d_lambdavec(3))                 :: lambdavec
  call ed_set_Hgeneral(Hvec,lambdavec)
end subroutine init_Hgeneral_symmetries_lattice_d3_c



!BREAK_SYMMETRY_BATH
subroutine break_symmetry_bath_site_c(bath,dim_bath,field,sgn,sav) bind(c, name='break_symmetry_bath_site')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                    :: dim_bath(1)
  real(c_double),dimension(dim_bath(1)) :: bath
  real(c_double),value                  :: field
  real(c_double),value                  :: sgn
  integer(c_int),value                  :: sav
  call ed_break_symmetry_bath(bath,field,sgn,i2l(sav))
end subroutine break_symmetry_bath_site_c
!
subroutine break_symmetry_bath_ineq_c(bath,dim_bath,field,sgn,sav) bind(c, name='break_symmetry_bath_ineq')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                :: dim_bath(2)
  real(c_double),dimension(dim_bath(1),dim_bath(2)) :: bath
  real(c_double),value                              :: field
  real(c_double),value                              :: sgn
  integer(c_int),value                              :: sav
  call ed_break_symmetry_bath(bath,field,sgn,i2l(sav))
end subroutine break_symmetry_bath_ineq_c



!SPIN_SYMMETRIZE_BATH
subroutine spin_symmetrize_bath_site_c(bath,dim_bath,sav) bind(c, name='spin_symmetrize_bath_site')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                    :: dim_bath(1)
  integer(c_int),value                  :: sav
  real(c_double),dimension(dim_bath(1)) :: bath
  call ed_spin_symmetrize_bath(bath,i2l(sav))
end subroutine spin_symmetrize_bath_site_c
!
subroutine spin_symmetrize_bath_ineq_c(bath,dim_bath,sav) bind(c, name='spin_symmetrize_bath_ineq')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                :: dim_bath(2)
  integer(c_int),value                              :: sav
  real(c_double),dimension(dim_bath(1),dim_bath(2)) :: bath
  call ed_spin_symmetrize_bath(bath,i2l(sav))
end subroutine spin_symmetrize_bath_ineq_c



!ORB_SYMMETRIZE_BATH
subroutine orb_symmetrize_bath_site_c(bath,dim_bath,orb1,orb2,sav) bind(c, name='orb_symmetrize_bath_site')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                    :: dim_bath(1)
  integer(c_int),value                  :: sav,orb1,orb2
  real(c_double),dimension(dim_bath(1)) :: bath
  call ed_orb_symmetrize_bath(bath,orb1,orb2,i2l(sav))
end subroutine orb_symmetrize_bath_site_c
!
subroutine orb_symmetrize_bath_ineq_c(bath,dim_bath,orb1,orb2,sav) bind(c, name='orb_symmetrize_bath_ineq')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                :: dim_bath(2)
  integer(c_int),value                              :: sav,orb1,orb2
  real(c_double),dimension(dim_bath(1),dim_bath(2)) :: bath
  call ed_orb_symmetrize_bath(bath,orb1,orb2,i2l(sav))
end subroutine orb_symmetrize_bath_ineq_c


!ORB_EQUALITY_BATH
subroutine orb_equality_bath_site_c(bath,dim_bath,indx,sav) bind(c, name='orb_equality_bath_site')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                       :: dim_bath(1)
  real(c_double),dimension(dim_bath(1))    :: bath
  integer(c_int),value                     :: indx,sav
  call ed_orb_equality_bath(bath,indx,i2l(sav))
end subroutine orb_equality_bath_site_c
!
subroutine orb_equality_bath_ineq_c(bath,dim_bath,indx,sav) bind(c, name='orb_equality_bath_ineq')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                   :: dim_bath(2)
  real(c_double),dimension(dim_bath(1),dim_bath(2))    :: bath
  integer(c_int),value                                 :: indx,sav
  call ed_orb_equality_bath(bath,indx,i2l(sav))
end subroutine orb_equality_bath_ineq_c





!PH_SYMMETRIZE_BATH
subroutine ph_symmetrize_bath_site_c(bath,dim_bath,sav) bind(c, name='ph_symmetrize_bath_site')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                    :: dim_bath(1)
  real(c_double),dimension(dim_bath(1)) :: bath
  integer(c_int),value                  :: sav
  call ed_ph_symmetrize_bath(bath,i2l(sav))
end subroutine ph_symmetrize_bath_site_c
!
subroutine ph_symmetrize_bath_ineq_c(bath,dim_bath,sav) bind(c, name='ph_symmetrize_bath_ineq')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                :: dim_bath(2)
  real(c_double),dimension(dim_bath(1),dim_bath(2)) :: bath
  integer(c_int),value                              :: sav
  call ed_ph_symmetrize_bath(bath,i2l(sav))
end subroutine ph_symmetrize_bath_ineq_c


!save array as bath
subroutine save_array_as_bath_site_c(bath,dim_bath) bind(c, name='save_array_as_bath_site')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                :: dim_bath(1)
  real(c_double),dimension(dim_bath(1))             :: bath
  call ed_save_array_as_bath(bath)
end subroutine save_array_as_bath_site_c

subroutine save_array_as_bath_ineq_c(bath,dim_bath) bind(c, name='save_array_as_bath_ineq')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                    :: dim_bath(2)
  real(c_double),dimension(dim_bath(1),dim_bath(2))     :: bath
  call ed_save_array_as_bath(bath)
end subroutine save_array_as_bath_ineq_c


