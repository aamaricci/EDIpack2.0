function delta_bath_array_hybrid(x,dmft_bath_,axis) result(Delta)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
  complex(8),dimension(:),intent(in)                                :: x          !complex  array for the frequency
  type(effective_bath)                                              :: dmft_bath_ !the current :f:var:`effective_bath` instance
  character(len=*),optional                                         :: axis       !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x))               :: Delta
  !
  integer                                                           :: i,ih,L
  integer                                                           :: iorb,jorb,ispin,jspin,ibath
  integer                                                           :: io,jo
  real(8),dimension(Nbath)                                          :: eps,dps,vps
  real(8),dimension(Norb,Nbath)                                     :: vops
  complex(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb)         :: Vk
  !
  real(8),dimension(Nspin,Nbath)                                    :: ehel
  real(8),dimension(Nspin,Nspin,Nbath)                              :: whel
  real(8),dimension(Nspin,Nspin,Norb,Nbath)                         :: wohel
  !
  complex(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb,size(x)) :: zeta
  complex(8),dimension(Nnambu*Nspin*Norb,size(x))                   :: z
  real(8),dimension(Nspin*Norb)                                     :: v
  complex(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb)         :: Hk
  complex(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb)         :: invH_k
  complex(8),dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb)         :: invH_knn
  complex(8),dimension(Nnambu*Norb,Nnambu*Norb)                     :: JJ
  character(len=1)                                                  :: axis_
  !
  axis_="m";if(present(axis))axis_=str(to_lower(axis))
  !
  Delta=zero
  !
  L = size(x)
  !
  select case(ed_mode)
  case default;stop "delta_bath_array_hybrid error: ed_mode not supported"
  case ("normal")
     do ispin=1,Nspin
        eps  = dmft_bath_%e(ispin,1     ,1:Nbath)
        vops = dmft_bath_%v(ispin,1:Norb,1:Nbath)
        do iorb=1,Norb
           do jorb=1,Norb
              do i=1,L
                 Delta(ispin,ispin,iorb,jorb,i) = sum( vops(iorb,:)*vops(jorb,:)/(x(i) - eps(:)) )
              enddo
           enddo
        enddo
     enddo
     !
  case ("superc")
     do ispin=1,Nspin
        eps  = dmft_bath_%e(ispin,1     ,1:Nbath)
        dps  = dmft_bath_%d(ispin,1     ,1:Nbath)
        vops = dmft_bath_%v(ispin,1:Norb,1:Nbath)
        do iorb=1,Norb
           do jorb=1,Norb
              select case(axis_)
              case default ;stop "delta_bath_array_hybrid error: axis not supported"         !mats
              case ("m")
                 do i=1,L
                    Delta(ispin,ispin,iorb,jorb,i) = &
                         -sum( vops(iorb,:)*vops(jorb,:)*(x(i) + eps(:))/(dimag(x(i))**2 + eps(:)**2 + dps(:)**2) )
                 enddo
              case ("r")
                 do i=1,L
                    Delta(ispin,ispin,iorb,jorb,i) = &
                         -sum( vops(iorb,:)*vops(jorb,:)*(x(i) + eps(:))/((x(i)*(-x(i)) + eps(:)**2 + dps(:)**2)) )
                 enddo
              end select
           enddo
        enddo
     enddo
     !
  case ("nonsu2")
     ehel  = dmft_bath_%e(1:Nspin,1,1:Nbath)
     wohel = get_Whyb_matrix(dmft_bath_%v(1:Nspin,1:Norb,1:Nbath),dmft_bath_%u(1:Nspin,1:Norb,1:Nbath))
     do iorb=1,Norb
        do jorb=1,Norb
           do ispin=1,Nspin
              do jspin=1,Nspin
                 do i=1,L
                    do ih=1,Nspin
                       Delta(ispin,jspin,iorb,jorb,i) = Delta(ispin,jspin,iorb,jorb,i) + &
                            sum( wohel(ispin,ih,iorb,:)*wohel(jspin,ih,jorb,:)/(x(i) - ehel(ih,:)) )
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
     !
  end select
end function delta_bath_array_hybrid
