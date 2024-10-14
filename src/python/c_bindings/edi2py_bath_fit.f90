!USE ED_BATH_FIT:

!single normal
subroutine chi2_fitgf_single_normal_n3_c(g,dim_g,bath,dim_bath,ispin,iorb,fmpi) bind(c, name='chi2_fitgf_single_normal_n3')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                                                  :: dim_g(3),dim_bath(1)
  complex(c_double_complex),dimension(dim_g(1),dim_g(2),dim_g(3))                     :: g
  real(c_double),dimension(dim_bath(1))                                               :: bath
  integer(c_int),value                                                                :: ispin,iorb,fmpi
  if(iorb>0)then
     call ed_chi2_fitgf(g,bath,ispin,iorb=iorb,fmpi=i2l(fmpi))
  else
     call ed_chi2_fitgf(g,bath,ispin,fmpi=i2l(fmpi))
  endif
end subroutine chi2_fitgf_single_normal_n3_c

subroutine chi2_fitgf_single_normal_n5_c(g,dim_g,bath,dim_bath,ispin,iorb,fmpi) bind(c, name='chi2_fitgf_single_normal_n5')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                                                           :: dim_g(5),dim_bath(1)
  complex(c_double_complex),dimension(dim_g(1),dim_g(2),dim_g(3),dim_g(4),dim_g(5))            :: g
  real(c_double),dimension(dim_bath(1))                                                        :: bath
  integer(c_int),value                                                                         :: ispin,iorb,fmpi
  if(iorb>0)then
     call ed_chi2_fitgf(g,bath,ispin,iorb=iorb,fmpi=i2l(fmpi))
  else
     call ed_chi2_fitgf(g,bath,ispin,fmpi=i2l(fmpi))
  endif
end subroutine chi2_fitgf_single_normal_n5_c

!single superc
subroutine chi2_fitgf_single_superc_n3_c(g,dim_g,f,dim_f,bath,dim_bath,ispin,iorb,fmpi) bind(c, name='chi2_fitgf_single_superc_n3')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                                                  :: dim_g(3),dim_f(3),dim_bath(1)
  complex(c_double_complex),dimension(dim_g(1),dim_g(2),dim_g(3))                     :: g
  complex(c_double_complex),dimension(dim_f(1),dim_f(2),dim_f(3))                     :: f
  real(c_double),dimension(dim_bath(1))                                               :: bath
  integer(c_int),value                                                                :: ispin,iorb,fmpi
  if(iorb>0)then
     call ed_chi2_fitgf(g,f,bath,ispin,iorb=iorb,fmpi=i2l(fmpi))
  else
     call ed_chi2_fitgf(g,f,bath,ispin,fmpi=i2l(fmpi))
  endif
end subroutine chi2_fitgf_single_superc_n3_c

subroutine chi2_fitgf_single_superc_n5_c(g,dim_g,f,dim_f,bath,dim_bath,ispin,iorb,fmpi) bind(c, name='chi2_fitgf_single_superc_n5')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                                                  :: dim_g(5),dim_f(5),dim_bath(1)
  complex(c_double_complex),dimension(dim_g(1),dim_g(2),dim_g(3),dim_g(4),dim_g(5))   :: g
  complex(c_double_complex),dimension(dim_f(1),dim_f(2),dim_f(3),dim_f(4),dim_f(5))   :: f
  real(c_double),dimension(dim_bath(1))                                               :: bath
  integer(c_int),value                                                                :: ispin,iorb,fmpi
  if(iorb>0)then
     call ed_chi2_fitgf(g,f,bath,ispin,iorb=iorb,fmpi=i2l(fmpi))
  else
     call ed_chi2_fitgf(g,f,bath,ispin,fmpi=i2l(fmpi))
  endif
end subroutine chi2_fitgf_single_superc_n5_c

!lattice normal
subroutine chi2_fitgf_lattice_normal_n3_c(g,dim_g,bath,dim_bath,ispin) bind(c, name='chi2_fitgf_lattice_normal_n3')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                                                  :: dim_g(3),dim_bath(2)
  complex(c_double_complex),dimension(dim_g(1),dim_g(2),dim_g(3))                     :: g
  real(c_double),dimension(dim_bath(1),dim_bath(2))                                   :: bath
  integer(c_int),value                                                                :: ispin
  call ed_chi2_fitgf(g,bath,ispin)
end subroutine chi2_fitgf_lattice_normal_n3_c

subroutine chi2_fitgf_lattice_normal_n4_c(g,dim_g,bath,dim_bath,ispin) bind(c, name='chi2_fitgf_lattice_normal_n4')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                                                  :: dim_g(4),dim_bath(2)
  complex(c_double_complex),dimension(dim_g(1),dim_g(2),dim_g(3),dim_g(4))            :: g
  real(c_double),dimension(dim_bath(1),dim_bath(2))                                   :: bath
  integer(c_int),value                                                                :: ispin
  call ed_chi2_fitgf(g,bath,ispin)
end subroutine chi2_fitgf_lattice_normal_n4_c

subroutine chi2_fitgf_lattice_normal_n6_c(g,dim_g,bath,dim_bath,ispin) bind(c, name='chi2_fitgf_lattice_normal_n6')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                                                                  :: dim_g(6),dim_bath(2)
  complex(c_double_complex),dimension(dim_g(1),dim_g(2),dim_g(3),dim_g(4),dim_g(5),dim_g(6))          :: g
  real(c_double),dimension(dim_bath(1),dim_bath(2))                                                   :: bath
  integer(c_int),value                                                                                :: ispin
  call ed_chi2_fitgf(g,bath,ispin)
end subroutine chi2_fitgf_lattice_normal_n6_c

!lattice superc
subroutine chi2_fitgf_lattice_superc_n3_c(g,dim_g,f,dim_f,bath,dim_bath,ispin) bind(c, name='chi2_fitgf_lattice_superc_n3')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                                                  :: dim_g(3),dim_f(3),dim_bath(2)
  complex(c_double_complex),dimension(dim_g(1),dim_g(2),dim_g(3))                     :: g
  complex(c_double_complex),dimension(dim_f(1),dim_f(2),dim_f(3))                     :: f
  real(c_double),dimension(dim_bath(1),dim_bath(2))                                   :: bath
  integer(c_int),value                                                                :: ispin
  call ed_chi2_fitgf(g,f,bath,ispin)
end subroutine chi2_fitgf_lattice_superc_n3_c

subroutine chi2_fitgf_lattice_superc_n4_c(g,dim_g,f,dim_f,bath,dim_bath,ispin) bind(c, name='chi2_fitgf_lattice_superc_n4')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                                                  :: dim_g(4),dim_f(4),dim_bath(2)
  complex(c_double_complex),dimension(dim_g(1),dim_g(2),dim_g(3),dim_g(4))            :: g
  complex(c_double_complex),dimension(dim_f(1),dim_f(2),dim_f(3),dim_f(4))            :: f
  real(c_double),dimension(dim_bath(1),dim_bath(2))                                   :: bath
  integer(c_int),value                                                                :: ispin
  call ed_chi2_fitgf(g,f,bath,ispin)
end subroutine chi2_fitgf_lattice_superc_n4_c

subroutine chi2_fitgf_lattice_superc_n6_c(g,dim_g,f,dim_f,bath,dim_bath,ispin) bind(c, name='chi2_fitgf_lattice_superc_n6')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                                                                  :: dim_g(6),dim_f(6),dim_bath(2)
  complex(c_double_complex),dimension(dim_g(1),dim_g(2),dim_g(3),dim_g(4),dim_g(5),dim_g(6))          :: g
  complex(c_double_complex),dimension(dim_f(1),dim_f(2),dim_f(3),dim_f(4),dim_f(5),dim_f(6))          :: f
  real(c_double),dimension(dim_bath(1),dim_bath(2))                                                   :: bath
  integer(c_int),value                                                                                :: ispin
  call ed_chi2_fitgf(g,f,bath,ispin)
end subroutine chi2_fitgf_lattice_superc_n6_c
