MODULE ED_CHI_DENS
  !Evaluates the impurity density susceptibility.
  !
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: str,reg,txtfy
  USE SF_LINALG,  only: inv,eigh,eye
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE
  USE ED_BATH
  USE ED_SETUP
  USE ED_SECTOR
  USE ED_HAMILTONIAN_NORMAL
  USE ED_AUX_FUNX

  implicit none
  private


  public :: build_chi_dens_normal

  integer                          :: istate,iorb,jorb,ispin,jspin
  integer                          :: isector
  real(8),allocatable              :: vvinit(:),vI(:),vJ(:)
  real(8),allocatable              :: alfa_(:),beta_(:)
  integer                          :: ialfa
  integer                          :: jalfa
  integer                          :: ipos,jpos
  integer                          :: i,j
  integer                          :: iph,i_el
  real(8)                          :: sgn,norm2
  real(8),dimension(:),allocatable :: v_state
  real(8)                          :: e_state


contains


  !+------------------------------------------------------------------+
  !                            DENS
  !PURPOSE  : Evaluate the Dens susceptibility \Chi_dens for a 
  ! \chi_ab = <n_a(\tau)n_b(0)>
  !+------------------------------------------------------------------+
  subroutine build_chi_dens_normal()
    !
    ! Evaluates the impurity Spin susceptibility :math:`\chi^{n}=\langle T_\tau n_a(\tau) n_b\rangle` in the Matsubara :math:`i\omega_n` and Real :math:`\omega` frequency axis as well as imaginary time :math:`\tau`.
    !
    ! As for the Green's function, the off-diagonal component of the the susceptibility is determined using an algebraic manipulation to ensure use of Hermitian operator in the dynamical Lanczos.
    !
    write(LOGfile,"(A)")"Get impurity dens Chi:"
    if(MPIMASTER)call start_timer(unit=LOGfile)
    do iorb=1,Norb
       call lanc_ed_build_densChi_diag(iorb)
    enddo
    !
    if(Norb>1)then
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             call lanc_ed_build_densChi_mix(iorb,jorb)
          end do
       end do
       !
       !
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             densChi_w(iorb,jorb,:)   = 0.5d0*(densChi_w(iorb,jorb,:) - densChi_w(iorb,iorb,:) - densChi_w(jorb,jorb,:))
             densChi_tau(iorb,jorb,:) = 0.5d0*(densChi_tau(iorb,jorb,:) - densChi_tau(iorb,iorb,:) - densChi_tau(jorb,jorb,:))
             densChi_iv(iorb,jorb,:)  = 0.5d0*(densChi_iv(iorb,jorb,:) - densChi_iv(iorb,iorb,:) - densChi_iv(jorb,jorb,:))
             !
             densChi_w(jorb,iorb,:)   = densChi_w(iorb,jorb,:)
             densChi_tau(jorb,iorb,:) = densChi_tau(iorb,jorb,:)
             densChi_iv(jorb,iorb,:)  = densChi_iv(iorb,jorb,:)
          enddo
       enddo
    endif
    !
    if(MPIMASTER)call stop_timer
    !
  end subroutine build_chi_dens_normal




  





  subroutine lanc_ed_build_densChi_diag(iorb)
    integer                     :: iorb
    type(sector)                :: sectorI,sectorJ
    !
    write(LOGfile,"(A)")"Get Chi_dens_l"//reg(txtfy(iorb))
    !
    do istate=1,state_list%size
       isector    =  es_return_sector(state_list,istate)
       e_state    =  es_return_energy(state_list,istate)
       v_state    =  es_return_dvec(state_list,istate)
       !
       if(MpiMaster)call build_sector(isector,sectorI)
       vvinit = apply_op_N(v_state,iorb,sectorI)
       call tridiag_Hv_sector_normal(isector,vvinit,alfa_,beta_,norm2)
       call add_to_lanczos_densChi(norm2,e_state,alfa_,beta_,iorb,iorb)
       deallocate(alfa_,beta_,vvinit)
       if(MpiMaster)call delete_sector(sectorI)
    enddo
    !
    if(allocated(v_state))deallocate(v_state)
    return
  end subroutine lanc_ed_build_densChi_diag




  subroutine lanc_ed_build_densChi_mix(iorb,jorb)
    integer                     :: iorb,jorb
    type(sector)                :: sectorI,sectorJ
    real(8)                     :: Niorb,Njorb
    !
    write(LOGfile,"(A)")"Get Chi_dens_mix_l"//reg(txtfy(iorb))//reg(txtfy(jorb))
    !
    do istate=1,state_list%size
       isector    =  es_return_sector(state_list,istate)
       e_state    =  es_return_energy(state_list,istate)
       v_state    =  es_return_dvec(state_list,istate)
       !
       !EVALUATE (N_jorb + N_iorb)|gs> = N_jorb|gs> + N_iorb|gs>
       if(MpiMaster)call build_sector(isector,sectorI)
       vI = apply_op_N(v_state,iorb,sectorI)
       vJ = apply_op_N(v_state,jorb,sectorI)
       call tridiag_Hv_sector_normal(isector,vI+vJ,alfa_,beta_,norm2)
       call add_to_lanczos_densChi(norm2,e_state,alfa_,beta_,iorb,jorb)
       deallocate(alfa_,beta_,vI,vJ)
       if(MpiMaster)call delete_sector(sectorI)
    enddo
    if(allocated(v_state))deallocate(v_state)
    return
  end subroutine lanc_ed_build_densChi_mix




  subroutine add_to_lanczos_densChi(vnorm2,Ei,alanc,blanc,iorb,jorb)
    integer                                    :: iorb,jorb
    real(8)                                    :: pesoF,pesoAB,pesoBZ,peso,vnorm2  
    real(8)                                    :: Ei,Ej,Egs,de
    integer                                    :: nlanc
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
    complex(8)                                 :: iw,chisp
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")&
         "DEBUG add_to_lanczos_densChi: add-up to GF"
#endif
    !
    Egs = state_list%emin       !get the gs energy
    !
    Nlanc = size(alanc)
    !
    pesoF  = vnorm2/zeta_function 
    pesoBZ = 1d0
    if(finiteT)pesoBZ = exp(-beta*(Ei-Egs))
    !
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,alanc)
       call Bcast_MPI(MpiComm,blanc)
    endif
#endif
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A)")&
         "DEBUG add_to_lanczos_densChi: LApack tridiagonalization"
#endif
    call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
    !
    do j=1,nlanc
       Ej     = diag(j)
       dE     = Ej-Ei
       pesoAB = Z(1,j)*Z(1,j)
       peso   = pesoF*pesoAB*pesoBZ
       ! the correct behavior for beta*dE << 1 is recovered only by assuming that v_n is still finite
       ! beta*dE << v_n for v_n--> 0 slower. First limit beta*dE--> 0 and only then v_n -->0.
       ! This ensures that the correct null contribution is obtained.
       ! So we impose that: if (beta*dE is larger than a small qty) we sum up the contribution, else
       ! we do not include the contribution (because we are in the situation described above).
       ! For the real-axis case this problem is circumvented by the usual i*0+ = xi*eps
       if(beta*dE > 1d-3)densChi_iv(iorb,jorb,0)=densChi_iv(iorb,jorb,0) + peso*2*(1d0-exp(-beta*dE))/dE 
       do i=1,Lmats
          densChi_iv(iorb,jorb,i)=densChi_iv(iorb,jorb,i) + peso*(1d0-exp(-beta*dE))*2d0*dE/(vm(i)**2+dE**2)
       enddo
       do i=0,Ltau
          densChi_tau(iorb,jorb,i)=densChi_tau(iorb,jorb,i) + exp(-tau(i)*dE)*peso
       enddo
       do i=1,Lreal
          densChi_w(iorb,jorb,i)=densChi_w(iorb,jorb,i) - &
               peso*(1d0-exp(-beta*dE))*(1d0/(dcmplx(vr(i),eps) - dE) - 1d0/(dcmplx(vr(i),eps) + dE))
       enddo
    enddo
    !
  end subroutine add_to_lanczos_densChi






END MODULE ED_CHI_DENS
























