function ed_get_g0and_function(x,bath_,axis) result(G0and)
  complex(8),dimension(:),intent(in)                  :: x
  type(effective_bath)                                :: dmft_bath_
  character(len=*),optional                           :: axis
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
  real(8),dimension(:)                                :: bath_
  logical                                             :: check
  character(len=4)                                    :: axis_
  axis_='mats';if(present(axis))axis_=axis
  check= check_bath_dimension(bath_)
  if(.not.check)stop "g0and_bath_mats_main_ error: wrong bath dimensions"
  call allocate_dmft_bath(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  G0and = g0and_bath_function(x,dmft_bath_)
  call deallocate_dmft_bath(dmft_bath_)
end function ed_get_g0and_function


function ed_get_f0and_function(x,bath_,axis) result(f0and)
  complex(8),dimension(:),intent(in)                  :: x
  type(effective_bath)                                :: dmft_bath_
  character(len=*),optional                           :: axis
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: f0and
  real(8),dimension(:)                                :: bath_
  logical                                             :: check
  character(len=4)                                    :: axis_
  axis_='mats';if(present(axis))axis_=axis
  check= check_bath_dimension(bath_)
  if(.not.check)stop "g0and_bath_mats_main_ error: wrong bath dimensions"
  call allocate_dmft_bath(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  F0and = f0and_bath_function(x,dmft_bath_)
  call deallocate_dmft_bath(dmft_bath_)
end function ed_get_f0and_function



function ed_get_delta_function(x,bath_,axis) result(Delta)
  complex(8),dimension(:),intent(in)                  :: x
  real(8),dimension(:)                                :: bath_
  character(len=*),optional                           :: axis
  type(effective_bath)                                :: dmft_bath_
  character(len=4)                                    :: axis_
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Delta
  logical                                             :: check
  axis_='mats';if(present(axis))axis_=axis
  check= check_bath_dimension(bath_)
  if(.not.check)stop "delta_bath_mats_main_ error: wrong bath dimensions"
  call allocate_dmft_bath(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  Delta = delta_bath_function(x,dmft_bath_,axis_)
  call deallocate_dmft_bath(dmft_bath_)
end function ed_get_delta_function

function ed_get_fdelta_function(x,bath_,axis) result(Delta)
  complex(8),dimension(:),intent(in)                  :: x
  real(8),dimension(:)                                :: bath_
  character(len=*),optional                           :: axis
  type(effective_bath)                                :: dmft_bath_
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Delta
  logical                                             :: check
  character(len=4)                                    :: axis_
  axis_='mats';if(present(axis))axis_=axis
  check= check_bath_dimension(bath_)
  if(.not.check)stop "delta_bath_mats_main_ error: wrong bath dimensions"
  call allocate_dmft_bath(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  Delta = fdelta_bath_function(x,dmft_bath_,axis_)
  call deallocate_dmft_bath(dmft_bath_)
end function ed_get_fdelta_function





