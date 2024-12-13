function fdelta_bath_array_normal(x,dmft_bath_,axis) result(Fdelta)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
  complex(8),dimension(:),intent(in)                                :: x !complex  array for the frequency
  type(effective_bath)                                              :: dmft_bath_ !the current :f:var:`effective_bath` instance
  character(len=*),optional                                         :: axis    !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x))               :: Fdelta
  character(len=4)                                                  :: axis_
  integer                                                           :: iorb,ispin,jorb
  real(8),dimension(Nbath)                                          :: eps,dps,vps
  integer                                                           :: i,L
  !
  axis_="mats";if(present(axis))axis_=str(axis)
  !
  Fdelta=zero
  !
  L = size(x)
  !
  do ispin=1,Nspin
     do iorb=1,Norb
        eps = dmft_bath_%e(ispin,iorb,1:Nbath)
        dps = dmft_bath_%d(ispin,iorb,1:Nbath)
        vps = dmft_bath_%v(ispin,iorb,1:Nbath)
        select case(axis_)
        case default         !mats
           do i=1,L
              Fdelta(ispin,ispin,iorb,iorb,i) = sum( dps(:)*vps(:)*vps(:)/(dimag(x(i))**2 + eps(:)**2 + dps(:)**2) )
           enddo
        case ("real")
           do i=1,L
              Fdelta(ispin,ispin,iorb,iorb,i) = sum( dps(:)*vps(:)*vps(:)/( x(i)*(-x(i)) + eps(:)**2 + dps(:)**2) )
           enddo
        end select
     enddo
  enddo
  !
end function fdelta_bath_array_normal