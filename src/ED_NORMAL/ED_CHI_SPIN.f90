MODULE ED_CHI_SPIN
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: str,reg,txtfy
  USE SF_LINALG,  only: inv,eigh,eye
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_EIGENSPACE
  USE ED_BATH
  USE ED_BATH_FUNCTIONS
  USE ED_SETUP
  USE ED_SECTOR
  USE ED_HAMILTONIAN_NORMAL
  implicit none
  private


  public :: build_chi_spin_normal

  integer                      :: istate,iorb,jorb,ispin
  integer                      :: isector
  real(8),allocatable          :: vvinit(:)
  real(8),allocatable          :: alfa_(:),beta_(:)
  integer                      :: ialfa
  integer                      :: jalfa
  integer                      :: ipos,jpos
  integer                      :: i,j
  real(8)                      :: sgn,norm2
  real(8),dimension(:),pointer :: state_dvec
  real(8)                      :: state_e




contains


  !+------------------------------------------------------------------+
  !                            SPIN
  !PURPOSE  : Evaluate the Spin susceptibility \Chi_spin for a 
  ! single orbital: \chi = <S_a(\tau)S_a(0)>
  ! note: as S_a is hermitian particle and holes contributions
  ! are identical so work out only one lanczos tridiag. work out the 
  ! reduction for both values of isign in the same call.
  !+------------------------------------------------------------------+
  subroutine build_chi_spin_normal()
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")&
         "DEBUG build_Chi_spin_normal: build spin-Chi"
#endif
    write(LOGfile,"(A)")"Get impurity spin Chi:"
    do iorb=1,Norb
       write(LOGfile,"(A)")"Get Chi_spin_l"//reg(txtfy(iorb))
       if(MPIMASTER)call start_timer()
       call lanc_ed_build_spinChi_diag(iorb)
       if(MPIMASTER)call stop_timer(unit=LOGfile)
#ifdef _DEBUG
       if(ed_verbose>1)write(Logfile,"(A)")""
#endif
    enddo
    !
    if(Norb>1)then
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             write(LOGfile,"(A)")"Get Chi_spin_mix_l"//reg(txtfy(iorb))//reg(txtfy(jorb))
             if(MPIMASTER)call start_timer()
             call lanc_ed_build_spinChi_mix(iorb,jorb)
             if(MPIMASTER)call stop_timer(unit=LOGfile)
#ifdef _DEBUG
             if(ed_verbose>1)write(Logfile,"(A)")""
#endif
          end do
       end do
       !
       !
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             select case(ed_diag_type)
             case default
                spinChi_w(iorb,jorb,:)   = 0.5d0*(spinChi_w(iorb,jorb,:) - spinChi_w(iorb,iorb,:) - spinChi_w(jorb,jorb,:))
                spinChi_tau(iorb,jorb,:) = 0.5d0*(spinChi_tau(iorb,jorb,:) - spinChi_tau(iorb,iorb,:) - spinChi_tau(jorb,jorb,:))
                spinChi_iv(iorb,jorb,:)  = 0.5d0*(spinChi_iv(iorb,jorb,:) - spinChi_iv(iorb,iorb,:) - spinChi_iv(jorb,jorb,:))
                !
             case ("full")
                ! The previous calculation is not needed in the FULL ED case
             end select
             !
             spinChi_w(jorb,iorb,:)   = spinChi_w(iorb,jorb,:)
             spinChi_tau(jorb,iorb,:) = spinChi_tau(iorb,jorb,:)
             spinChi_iv(jorb,iorb,:)  = spinChi_iv(iorb,jorb,:)
          enddo
       enddo
    endif
    !
  end subroutine build_chi_spin_normal






  !################################################################
  !################################################################
  !################################################################
  !################################################################






  subroutine lanc_ed_build_spinChi_diag(iorb)
    integer                     :: iorb
    type(sector)                :: sectorI,sectorJ
    !
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")&
         "DEBUG lanc_ed_build_spinChi diag: Lanczos build spin Chi l"//str(iorb)
#endif
    !
    if(ed_total_ud)then
       ialfa = 1
       ipos  = iorb
    else
       ialfa = iorb
       ipos  = 1
    endif
    !
    do istate=1,state_list%size
       isector    =  es_return_sector(state_list,istate)
       state_e    =  es_return_energy(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          state_dvec => es_return_dvector(MpiComm,state_list,istate)
       else
          state_dvec => es_return_dvector(state_list,istate)
       endif
#else
       state_dvec => es_return_dvector(state_list,istate)
#endif
       !
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          if(ed_verbose>=3)write(LOGfile,"(A20,I6,20I4)")&
               'Apply Sz',isector,sectorI%Nups,sectorI%Ndws
          allocate(vvinit(sectorI%Dim)) ; vvinit=zero
          do i=1,sectorI%Dim
             call apply_op_Sz(i,sgn,ipos,ialfa,sectorI)            
             vvinit(i) = sgn*state_dvec(i)
          enddo
          call delete_sector(sectorI)
       else
          allocate(vvinit(1));vvinit=0.d0
       endif
       !
       call tridiag_Hv_sector_normal(isector,vvinit,alfa_,beta_,norm2)
       call add_to_lanczos_spinChi(norm2,state_e,alfa_,beta_,iorb,iorb)
       deallocate(alfa_,beta_)
       if(allocated(vvinit))deallocate(vvinit)
       !
#ifdef _MPI
       if(MpiStatus)then
          if(associated(state_dvec))deallocate(state_dvec)
       else
          if(associated(state_dvec))nullify(state_dvec)
       endif
#else
       if(associated(state_dvec))nullify(state_dvec)
#endif
    enddo
    return
  end subroutine lanc_ed_build_spinChi_diag




  !################################################################




  subroutine lanc_ed_build_spinChi_mix(iorb,jorb)
    integer                     :: iorb,jorb
    type(sector)                :: sectorI,sectorJ
    real(8)                     :: Siorb,Sjorb
    !
    !
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")&
         "DEBUG lanc_ed_build_spinChi mix: Lanczos build spin Chi l"//str(iorb)//",m"//str(jorb)
#endif
    !    
    if(ed_total_ud)then
       ialfa = 1
       jalfa = 1
       ipos  = iorb
       jpos  = jorb
    else
       ialfa = iorb
       jalfa = jorb
       ipos  = 1
       jpos  = 1
    endif
    !
    do istate=1,state_list%size
       isector    =  es_return_sector(state_list,istate)
       state_e    =  es_return_energy(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          state_dvec => es_return_dvector(MpiComm,state_list,istate)
       else
          state_dvec => es_return_dvector(state_list,istate)
       endif
#else
       state_dvec => es_return_dvector(state_list,istate)
#endif
       !
       !EVALUATE (Sz_jorb + Sz_iorb)|gs> = Sz_jorb|gs> + Sz_iorb|gs>
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          if(ed_verbose>=3)write(LOGfile,"(A20,I6,20I4)")&
               'From sector',isector,sectorI%Nups,sectorI%Ndws
          if(ed_verbose>=3)write(LOGfile,"(A20,I15)")'Apply (Sz_a + Sz_b)',isector
          allocate(vvinit(sectorI%Dim)) ; vvinit=zero
          do i=1,sectorI%Dim
             call apply_op_Sz(i,Siorb,ipos,ialfa,sectorI)
             call apply_op_Sz(i,Sjorb,jpos,jalfa,sectorI)
             sgn       = Siorb + Sjorb
             vvinit(i) = sgn*state_dvec(i)
          enddo
          call delete_sector(sectorI)
       else
          allocate(vvinit(1));vvinit=0.d0
       endif
       !
       call tridiag_Hv_sector_normal(isector,vvinit,alfa_,beta_,norm2)
       call add_to_lanczos_spinChi(norm2,state_e,alfa_,beta_,iorb,iorb)
       deallocate(alfa_,beta_)
       if(allocated(vvinit))deallocate(vvinit)
       !
#ifdef _MPI
       if(MpiStatus)then
          if(associated(state_dvec))deallocate(state_dvec)
       else
          if(associated(state_dvec))nullify(state_dvec)
       endif
#else
       if(associated(state_dvec))nullify(state_dvec)
#endif
    enddo
    return
  end subroutine lanc_ed_build_spinChi_mix





  !################################################################


  subroutine add_to_lanczos_spinChi(vnorm2,Ei,alanc,blanc,iorb,jorb)
    real(8)                                    :: vnorm2,Ei,Ej,Egs,pesoF,pesoAB,pesoBZ,de,peso
    integer                                    :: nlanc
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    integer                                    :: iorb,jorb
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
    complex(8)                                 :: iw,chisp
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")&
         "DEBUG add_to_lanczos_spinChi: add-up to GF"
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
         "DEBUG add_to_lanczos_spinChi: LApack tridiagonalization"
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
       if(beta*dE > 1d-3)spinChi_iv(iorb,jorb,0)=spinChi_iv(iorb,jorb,0) + peso*2*(1d0-exp(-beta*dE))/dE 
       do i=1,Lmats
          spinChi_iv(iorb,jorb,i)=spinChi_iv(iorb,jorb,i) + peso*(1d0-exp(-beta*dE))*2d0*dE/(vm(i)**2+dE**2)
       enddo
       do i=0,Ltau
          spinChi_tau(iorb,jorb,i)=spinChi_tau(iorb,jorb,i) + exp(-tau(i)*dE)*peso
       enddo
       do i=1,Lreal
          spinChi_w(iorb,jorb,i)=spinChi_w(iorb,jorb,i) - peso*(1d0-exp(-beta*dE))*(1d0/(dcmplx(vr(i),eps) - dE) - 1d0/(dcmplx(vr(i),eps) + dE))
       enddo
    enddo
  end subroutine add_to_lanczos_spinChi



END MODULE ED_CHI_SPIN
























