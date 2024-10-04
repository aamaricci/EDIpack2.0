!ED_MAIN:
subroutine init_solver_site_c(bath,dim_bath) bind(c, name='init_solver_site')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t),dimension(1),intent(in)                 :: dim_bath
  real(c_double),dimension(dim_bath(1)),intent(inout)        :: bath
  call ed_init_solver(bath)
end subroutine init_solver_site_c
!
subroutine init_solver_ineq_c(bath,dim_bath) bind(c, name='init_solver_ineq')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t),dimension(2),intent(in)                           :: dim_bath
  real(c_double),dimension(dim_bath(1),dim_bath(2)),intent(inout)      :: bath
  call ed_init_solver(bath)
end subroutine init_solver_ineq_c


subroutine solve_site_c(bath,dim_bath,sflag,fmpi) bind(c, name='solve_site')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t),dimension(1),intent(in)                                             :: dim_bath
  integer(c_int),value                                                                   :: sflag,fmpi
  real(c_double),dimension(dim_bath(1)),intent(in)                                       :: bath
  call ed_solve(bath,sflag=i2l(sflag),fmpi=i2l(fmpi))
end subroutine solve_site_c
!
subroutine solve_ineq_c(bath,dim_bath,mpi_lanc,iflag) bind(c, name='solve_ineq')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t),dimension(2),intent(in)                                                           :: dim_bath
  integer(c_int),value                                                                                 :: mpi_lanc,iflag
  real(c_double),dimension(dim_bath(1),dim_bath(2)),intent(in)                                         :: bath
  integer                                                                                              :: Nineq
  Nineq = size(bath,1)
  call ed_solve(bath,mpi_lanc=i2l(mpi_lanc),iflag=i2l(iflag))
end subroutine solve_ineq_c

subroutine finalize_solver_c(Nineq) bind(c, name='finalize_solver')
  use, intrinsic :: iso_c_binding
  integer(c_int),value                            :: Nineq
  if (Nineq == 0) then
    call ed_finalize_solver()
  else
    call ed_finalize_solver(Nineq)
  endif
end subroutine finalize_solver_c
