function invf0_bath_array_hyrege(x,dmft_bath_,axis) result(F0and)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
  complex(8),dimension(:),intent(in)                  :: x          !complex  array for the frequency
  type(effective_bath)                                :: dmft_bath_ !the current :f:var:`effective_bath` instance
  character(len=*),optional                           :: axis       !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis    
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: F0and
  character(len=4)                                    :: axis_
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Fdelta
  integer                                             :: iorb,jorb,ispin,L
  !
  axis_="mats";if(present(axis))axis_=str(axis)
  !
  F0and=zero
  !
  L = size(x)
  !
  Fdelta= fdelta_bath_array(x,dmft_bath_,axis_)
  do ispin=1,Nspin
     do iorb=1,Norb
        do jorb=1,Norb
           F0and(ispin,ispin,iorb,jorb,:) = -Fdelta(ispin,ispin,iorb,jorb,:)
        enddo
     enddo
  enddo
  !
end function invf0_bath_array_hyrege