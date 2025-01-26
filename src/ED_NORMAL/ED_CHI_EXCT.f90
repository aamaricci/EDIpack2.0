MODULE ED_CHI_EXCT
  !Evaluates the impurity excitonc susceptibility.
  !
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: str,free_unit,reg,free_units,txtfy
  USE SF_LINALG,  only: inv,eigh,eye
  USE SF_SP_LINALG, only: sp_lanc_tridiag
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_IO                     !< this contains the routine to print GF,Sigma and G0
  USE ED_EIGENSPACE
  USE ED_BATH
  USE ED_SETUP
  USE ED_SECTOR
  USE ED_HAMILTONIAN_NORMAL
  USE ED_AUX_FUNX

  implicit none
  private


  public :: build_chi_exct_normal

  integer                          :: istate,iorb,jorb,ispin,jspin
  integer                          :: isector,jsector,ksector
  real(8),allocatable              :: vvinit(:)
  real(8),allocatable              :: alfa_(:),beta_(:)
  integer                          :: ipos,jpos
  integer                          :: i,j,k
  real(8)                          :: sgn,norm2
  real(8),dimension(:),allocatable :: v_state
  real(8)                          :: e_state

contains


  !+------------------------------------------------------------------+
  !                            EXCITON
  !PURPOSE  : Evaluate the Exciton susceptibility \Chi_exct for a 
  ! \chi_ab = <O*_a(\tau)O_b(0)>
  ! a/=b
  ! Singlet: \sum_\sigma <C^+_{a\sigma}C_{b\sigma} 
  ! Triplet: \sum_{\sigma\rho} C^+_{a\sigma} \tau_{\sigma\rho} C_{b\rho}
  !+------------------------------------------------------------------+
  subroutine build_chi_exct_normal()
    !
    !
    ! Evaluates the impurity exciton-exciton susceptibility :math:`\chi^{X}_{ab}=\langle T_\tau X^\dagger_{ab}(\tau) X_{ab}\rangle` in the Matsubara :math:`i\omega_n` and Real :math:`\omega` frequency axis, the imaginary time :math:`\tau` as well as the singlet and triplet components of the operator. 
    !
    ! As for the Green's function, the off-diagonal component of the the susceptibility is determined using an algebraic manipulation to ensure use of Hermitian operator in the dynamical Lanczos. 
    !
    if(Norb>1)then
       write(LOGfile,"(A)")"Get impurity exciton Chi:"
       if(MPIMASTER)call start_timer(unit=LOGfile)
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             call lanc_ed_build_exctChi_singlet(iorb,jorb)
             call lanc_ed_build_exctChi_tripletXY(iorb,jorb)
             call lanc_ed_build_exctChi_tripletZ(iorb,jorb)
             !
             exctChi_w(0:,jorb,iorb,:)   = exctChi_w(0:,iorb,jorb,:)
             exctChi_tau(0:,jorb,iorb,:) = exctChi_tau(0:,iorb,jorb,:)
             exctChi_iv(0:,jorb,iorb,:)  = exctChi_iv(0:,iorb,jorb,:)
          end do
       end do
       if(MPIMASTER)call stop_timer
    endif
  end subroutine build_chi_exct_normal






  
  ! \chi_ab  = <Delta*_ab(\tau)Delta_ab(0)>
  !\Delta_ab = \sum_\sigma C^+_{a\sigma}C_{b\sigma}
  subroutine lanc_ed_build_exctChi_singlet(iorb,jorb)
    integer      :: iorb,jorb
    type(sector) :: sectorI
    real(8),dimension(:),allocatable :: vup,vdw,vtmp
    !
    write(LOGfile,"(A)")"Get singlet Chi_exct_l"//reg(txtfy(iorb))//reg(txtfy(jorb))
    !
    do istate=1,state_list%size
       isector    =  es_return_sector(state_list,istate)
       e_state    =  es_return_energy(state_list,istate)
       v_state    =  es_return_dvec(state_list,istate)
       !
       ksector = getCsector(1,2,isector)       
       if(ksector/=0)then
          !C_b,dw|gs>=|tmp>
          vtmp = apply_op_C(v_state,jorb,2,isector,ksector)
          !C^+_a,dw|tmp>=|vvinit>
          vdw  = apply_op_CDG(vtmp,iorb,2,ksector,isector)
       endif
       ksector = getCsector(1,1,isector)
       if(ksector/=0)then
          !C_b,up|gs>=|tmp>
          vtmp = apply_op_C(v_state,jorb,1,isector,ksector)
          !C^+_a,up|tmp>=|vvinit>
          vup  = apply_op_CDG(vtmp,iorb,1,ksector,isector)
       endif
       call tridiag_Hv_sector_normal(isector,vup+vdw,alfa_,beta_,norm2)
       call add_to_lanczos_exctChi(norm2,e_state,alfa_,beta_,iorb,jorb,0)
       deallocate(alfa_,beta_,vup,vdw)
       if(allocated(v_state))deallocate(v_state)
    enddo
    return
  end subroutine lanc_ed_build_exctChi_singlet







  ! \chi_ab  = <Z_ab(\tau)Z_ab(0)>
  !Z_ab = \sum_sp C^+_{as}.tau^z_{sp}.C_{bp}
  subroutine lanc_ed_build_exctChi_tripletZ(iorb,jorb)
    integer                          :: iorb,jorb
    real(8),dimension(:),allocatable :: vup,vdw,vtmp
    !
    write(LOGfile,"(A)")"Get triplet Z Chi_exct_l"//reg(txtfy(iorb))//reg(txtfy(jorb))
    !
    !
    do istate=1,state_list%size
       isector    =  es_return_sector(state_list,istate)
       e_state    =  es_return_energy(state_list,istate)
       v_state    =  es_return_dvec(state_list,istate)
       !
       !Z - Component:
       !Z_{ab}= C^+_{a,up}C_{b,up} - C^+_{a,dw}C_{b,dw}
       ksector = getCsector(1,2,isector)
       if(ksector/=0)then
          !C_b,dw  |gs> =|tmp>
          vtmp = apply_op_C(v_state,jorb,2,isector,ksector)
          !C^+_a,dw|tmp>=|vvinit>
          vdw  = apply_op_CDG(vtmp,iorb,2,ksector,isector)
       endif
       ksector = getCsector(1,1,isector)
       if(ksector/=0)then
          !C_b,up  |gs> =|tmp>
          vtmp = apply_op_C(v_state,jorb,1,isector,ksector)
          !C^+_a,up|tmp>=|vvinit>
          vup  = apply_op_CDG(vtmp,iorb,1,ksector,isector)
       endif
       call tridiag_Hv_sector_normal(isector,vup-vdw,alfa_,beta_,norm2)
       call add_to_lanczos_exctChi(norm2,e_state,alfa_,beta_,iorb,jorb,2)
       deallocate(alfa_,beta_,vup,vdw)
       if(allocated(v_state))deallocate(v_state)
    enddo
    return
  end subroutine lanc_ed_build_exctChi_tripletZ





  ! \chi_ab  = <O_ab(\tau)O_ab(0)>
  ! O_ab = \sum_sp C^+_{as}.tau^o_{sp}.C_{bp} with o=X,Y
  ! O_ab|0> X:=   [C^+_{a,up}C_{b,dw} + C^+_{a,dw}C_{b,up}]|0>
  !         Y:= -i[C^+_{a,up}C_{b,dw} - C^+_{a,dw}C_{b,up}]|0>
  !         X:=   [P_{up,dw} +  P_{dw,up}]|0> = |v> + |w> 
  !         X:= -i[P_{up,dw} -  P_{dw,up}]|0> = |v> + |w>
  ! If |0>\in\SS(N_up,N_dw) => |v>\in\SS(N_up+1,N_dw-1), |w>\in\SS(N_up-1,N_dw+1)
  ! so that the sum |v> + |w> can not be accumulated onto a single vector |vvinit>
  ! Yet, we can recast the \Chi_ab expression in:
  ! \chi_ab = (<w| + <v|)[z-H]^{-1}(|v> + |w>)
  ! the direct terms: <v|[z-H]^{-1}|v> and <w|[z-H]^{-1}|w> are evaluated as usual.
  ! the mixed terms: <v|[z-H]^{-1}|w> and <w|[z-H]^{-1}|v> are indeed null.
  ! Proof:
  ! |v> and |w> belong to different sectors. H has a sector-block structure and so
  ! does its inverse H^{-1} made of the inverse of each block. 
  ! The expected values <v|H^{-1}|w> are taken across different sectors, but unless
  ! spin-conservation is broken these sectors are not connected and, as such, these
  ! numbers have to be zero.
  ! Note that while H can have a sparse nature, its inverse in general does not.
  ! For this reason the same argument as above DOES NOT apply to the Z case as
  ! in that case |v> and |w> belong to the same sector (the same as |0>) and the
  ! mixed term is in general non null.
  subroutine lanc_ed_build_exctChi_tripletXY(iorb,jorb)
    integer                          :: iorb,jorb
    real(8),dimension(:),allocatable :: vtmp
    !
    write(LOGfile,"(A)")"Get triplet XY Chi_exct_l"//reg(txtfy(iorb))//reg(txtfy(jorb))
    !
    do istate=1,state_list%size
       isector    =  es_return_sector(state_list,istate)
       e_state    =  es_return_energy(state_list,istate)
       v_state    =  es_return_dvec(state_list,istate)
       !
       !X - Component == Y -Component 
       !X_{ab}= C^+_{a,up}C_{b,dw} + C^+_{a,dw}C_{b,up}
       !
       !C^+_{a,dw}C_{b,up}:
       ksector = getCsector(1,1,isector)
       if(ksector/=0)then
          jsector = getCDGsector(1,2,ksector)
          if(jsector/=0)then
             !C_{b,up}|gs>   =|tmp>
             vtmp   = apply_op_C(v_state,jorb,1,isector,ksector)
             !C^+_{a,dw}|tmp>=|vvinit>
             vvinit = apply_op_CDG(vtmp,iorb,2,ksector,jsector)
             call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
             call add_to_lanczos_exctChi(norm2,e_state,alfa_,beta_,iorb,jorb,1)
             deallocate(alfa_,beta_,vtmp,vvinit)
          endif
       endif
       !
       !C^+_{a,up}C_{b,dw}:
       ksector = getCsector(1,2,isector)
       if(ksector/=0)then
          jsector = getCDGsector(1,1,ksector)
          if(jsector/=0)then
             !C_{b,dw}|gs>   =|tmp>
             vtmp   = apply_op_C(v_state,jorb,2,isector,ksector)
             !C^+_{a,up}|tmp>=|vvinit>
             vvinit = apply_op_CDG(vtmp,iorb,1,ksector,jsector)
             call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
             call add_to_lanczos_exctChi(norm2,e_state,alfa_,beta_,iorb,jorb,1)
             deallocate(alfa_,beta_,vtmp,vvinit)
          endif
       endif
       if(allocated(v_state))deallocate(v_state)
       !
    enddo
    return
  end subroutine lanc_ed_build_exctChi_tripletXY






  subroutine add_to_lanczos_exctChi(vnorm2,Ei,alanc,blanc,iorb,jorb,indx)
    integer                                    :: iorb,jorb,indx
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
         "DEBUG add_to_lanczos_exctChi: add-up to GF"
#endif
    !
    if(vnorm2==0)return
    !
    Egs = state_list%emin       !get the gs energy
    !
    Nlanc = size(alanc)
    !
    if(.not.any([0,1,2]==indx))stop "add_to_lanczos_exctChi ERROR: indx/=any[0,1,2]"
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
         "DEBUG add_to_lanczos_exctChi: LApack tridiagonalization"
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
       if(beta*dE > 1d-3)exctChi_iv(indx,iorb,jorb,0)=exctChi_iv(indx,iorb,jorb,0) + peso*2*(1d0-exp(-beta*dE))/dE 
       do i=1,Lmats
          exctChi_iv(indx,iorb,jorb,i)=exctChi_iv(indx,iorb,jorb,i) + &
               peso*(1d0-exp(-beta*dE))*2d0*dE/(vm(i)**2+dE**2)
       enddo
       do i=0,Ltau
          exctChi_tau(indx,iorb,jorb,i)=exctChi_tau(indx,iorb,jorb,i) + exp(-tau(i)*dE)*peso
       enddo
       do i=1,Lreal
          exctChi_w(indx,iorb,jorb,i)=exctChi_w(indx,iorb,jorb,i) - &
               peso*(1d0-exp(-beta*dE))*(1d0/(dcmplx(vr(i),eps) - dE) - 1d0/(dcmplx(vr(i),eps) + dE))
       enddo
    enddo
    !
  end subroutine add_to_lanczos_exctChi



END MODULE ED_CHI_EXCT




























