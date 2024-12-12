function delta_bath_array_general(x,dmft_bath_,axis) result(Delta)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
  complex(8),dimension(:),intent(in)                                :: x          !complex  array for the frequency
  type(effective_bath)                                              :: dmft_bath_ !the current :f:var:`effective_bath` instance
  character(len=*),optional                                         :: axis       !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x))               :: Delta
  integer                                                           :: i,ih,L
  integer                                                           :: iorb,jorb,ispin,jspin,ibath
  integer                                                           :: io,jo
  complex(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb)         :: Vk
  complex(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb,size(x)) :: zeta
  complex(8),dimension(Nnambu*Nspin*Norb,size(x))                   :: z
  real(8),dimension(Nspin*Norb)                                     :: v
  complex(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb)         :: Hk
  complex(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb)         :: invH_k
  complex(8),dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb)         :: invH_knn
  character(len=4)                                                  :: axis_
  !
  axis_="mats";if(present(axis))axis_=str(axis)
  !
  Delta=zero
  !
  L = size(x)
  !
  select case(ed_mode)
  case default;stop "delta_bath_array_general error: ed_mode not supported"
  case ("normal","nonsu2")             !normal OR nonsu2
     invH_k=zero
     do ibath=1,Nbath
        Vk = one*diag(dmft_bath_%item(ibath)%vg(:))
        Hk = nn2so_reshape(Hgeneral_build(dmft_bath_%item(ibath)%lambda),Nspin,Norb)
        do i=1,L
           invH_k   = zeye(Nspin*Norb)*x(i) - Hk
           call inv(invH_k)
           invH_k   = matmul(matmul(Vk,invH_k),Vk)
           invH_knn = so2nn_reshape(invH_k,Nspin,Norb)
           Delta(:,:,:,:,i)=Delta(:,:,:,:,i) + invH_knn
        enddo
     enddo
  case ("superc")
     !
     Z = zeta_superc(x,0d0,axis_)
     !
     do ibath=1,Nbath
        v  = dmft_bath_%item(ibath)%vg
        Hk = nn2so_reshape(Hreplica_build(dmft_bath_%item(ibath)%lambda),Nnambu*Nspin,Norb)
        Vk = kron( pauli_sigma_z, one*diag(v) )
        do i=1,L
           invH_k   = diag(Z(:,i)) - Hk
           call inv(invH_k)
           invH_k   = matmul(matmul(Vk,invH_k),Vk)
           invH_knn = so2nn_reshape(invH_k,Nnambu*Nspin,Norb)
           Delta(1,1,:,:,i)=Delta(1,1,:,:,i) + invH_knn(1,1,:,:)
        enddo
     enddo
  end select
  !
end function delta_bath_array_general
