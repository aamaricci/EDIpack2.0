function f0and_bath_array_hyrege(x,dmft_bath_,axis) result(F0and)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
  complex(8),dimension(:),intent(in)                  :: x !complex  array for the frequency
  type(effective_bath)                                :: dmft_bath_ !the current :f:var:`effective_bath` instance
  character(len=*),optional                           :: axis!string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis    
  character(len=4)                                    :: axis_
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: F0and
  !
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Delta,Fdelta
  integer                                             :: iorb,jorb,ispin,i,L
  complex(8),dimension(Nnambu*Nspin*Norb,size(x))     :: z
  complex(8),dimension(:,:),allocatable               :: fgorb,zeta
  !
  axis_="mats";if(present(axis))axis_=str(axis)
  !
  F0and=zero
  !
  L = size(x)
  !
  allocate(fgorb(2*Norb,2*Norb),zeta(2*Norb,2*Norb))
  Delta =  delta_bath_array(x,dmft_bath_,axis_)
  Fdelta= fdelta_bath_array(x,dmft_bath_,axis_)
  z     = zeta_superc(x,xmu,axis_)
  !
  select case(axis_)
  case default;stop "f0and_bath_array error: axis_ not support"
  case ("mats")
     do ispin=1,Nspin   !==1
        do i=1,L
           zeta = diag(z(:,i))
           fgorb= zero
           do iorb=1,Norb
              do jorb=1,Norb
                 fgorb(iorb,jorb)           = zeta(iorb,jorb)           - impHloc(ispin,ispin,iorb,jorb)  - Delta(ispin,ispin,iorb,jorb,i)
                 fgorb(iorb,jorb+Norb)      = zeta(iorb,jorb+Norb)                                        - Fdelta(ispin,ispin,iorb,jorb,i)
                 fgorb(iorb+Norb,jorb)      = zeta(iorb+Norb,jorb)                                        - Fdelta(ispin,ispin,iorb,jorb,i)
                 fgorb(iorb+Norb,jorb+Norb) = zeta(iorb+Norb,jorb+Norb) + conjg(impHloc(ispin,ispin,iorb,jorb)) + conjg( Delta(ispin,ispin,iorb,jorb,i) )
              enddo
           enddo
           call inv(fgorb)
           F0and(ispin,ispin,:,:,i) = fgorb(1:Norb,1+Norb:Norb+Norb)
        enddo
     enddo
  case("real")
     do ispin=1,Nspin   !==1
        do i=1,L
           zeta = diag(z(:,i))
           fgorb= zero
           do iorb=1,Norb
              do jorb=1,Norb
                 fgorb(iorb,jorb)           = zeta(iorb,jorb)           - impHloc(ispin,ispin,iorb,jorb)  - Delta(ispin,ispin,iorb,jorb,i)
                 fgorb(iorb,jorb+Norb)      = zeta(iorb,jorb+Norb)                                        - Fdelta(ispin,ispin,iorb,jorb,i)
                 fgorb(iorb+Norb,jorb)      = zeta(iorb+Norb,jorb)                                        - Fdelta(ispin,ispin,iorb,jorb,i)
                 fgorb(iorb+Norb,jorb+Norb) = zeta(iorb+Norb,jorb+Norb) + conjg(impHloc(ispin,ispin,iorb,jorb))  + conjg( Delta(ispin,ispin,iorb,jorb,L-i+1) )
              enddo
           enddo
           call inv(fgorb)
           F0and(ispin,ispin,:,:,i) = fgorb(1:Norb,1+Norb:Norb+Norb)
        enddo
     enddo
  end select
  deallocate(fgorb,zeta)
  !
end function f0and_bath_array_hyrege
