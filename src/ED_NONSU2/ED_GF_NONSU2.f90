MODULE ED_GF_NONSU2
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: str,reg,txtfy
  USE SF_LINALG,  only: inv,eigh,eye
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_EIGENSPACE
  USE ED_BATH
  USE ED_SETUP
  USE ED_SECTOR
  USE ED_HAMILTONIAN_NONSU2
  implicit none
  private


  public :: build_gf_nonsu2
  public :: build_sigma_nonsu2
  public :: rebuild_gf_nonsu2



  integer                               :: istate
  integer                               :: isector,jsector,ksector
  integer                               :: idim,jdim
  integer                               :: isz,jsz
  integer                               :: ialfa,jalfa
  integer                               :: i,j,m,r
  integer                               :: iph,i_el
  real(8)                               :: sgn,norm2
  complex(8)                            :: cnorm2

  complex(8),allocatable                :: vvinit(:)
  real(8),allocatable                   :: alfa_(:),beta_(:)



  !Lanczos shared variables
  !=========================================================
  complex(8),dimension(:),allocatable   :: v_state
  real(8)                               :: e_state



contains



  !+------------------------------------------------------------------+
  !                            NONSU2
  !+------------------------------------------------------------------+
  subroutine build_gf_nonsu2()
    !
    !
    !Evaluates the impurity electrons Green's functions :math:`G(z)` using dynamical Lanczos method. The result is stored in rank-5 arrays :f:var:`impgmats`, :f:var:`impgreal` , :f:var:`impfmats` , :f:var:`impfreal` of dimensions [ |Nspin| , |Nspin| , |Norb| , |Norb| , :f:var:`Lmats` / :f:var:`Lreal` ]
    !
    !The off-diagonal components of :math:`G_{ab\sigma\sigma'}` with :math:`a \neq b` and :math:`\sigma,\sigma'=\uparrow, \downarrow` are obtained using algebraic manipulation to ensure working with hermitian conjugate operators in the dynamical Lanczos procedure.  
    !
    !The weights and the poles obtained in this procedure are saved in a hierarchical data structure (for every state, every channel (creation or annihilation of excitations, normal or anomalous) and every degree of freedom) :f:var:`impgmatrix` of type :f:var:`gfmatrix`. 
    !
    ! .. _j.cpc.2021.108261: https://doi.org/10.1016/j.cpc.2021.108261
    !
    integer                                     :: iorb,jorb,ispin,jspin,i,io,jo
    logical                                     :: MaskBool
    logical(8),dimension(Nspin,Nspin,Norb,Norb) :: Hmask
    !
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG build_gf NONSU2: build GFs"
#endif
    !    

    Hmask=.true. !Aren't we sure about hermiticity? -> Hmask=Hreplica_mask(wdiag=.false.,uplo=.true.)
    if(.not.ed_all_g)then
       select case(bath_type)
       case("replica");Hmask=Hreplica_mask(wdiag=.true.,uplo=.false.)
       case("general");Hmask=Hgeneral_mask(wdiag=.true.,uplo=.false.)
       case default;stop "ERROR: ED_ALL_G=FALSE AND BATH_TYPE!=REPLICA/GENERAL"
       end select
    end if
    write(LOGfile,"(A)")"Get mask(G):"
    do ispin=1,Nspin
       do iorb=1,Norb
          write(LOGfile,*)((Hmask(ispin,jspin,iorb,jorb),jorb=1,Norb),jspin=1,Nspin)
       enddo
    enddo
    !
    if(MPIMASTER)call start_timer(unit=LOGfile)
    !same orbital, same spin GF: G_{aa}^{ss}(z)    
    do ispin=1,Nspin
       do iorb=1,Norb
          call allocate_GFmatrix(impGmatrix(ispin,ispin,iorb,iorb),Nstate=state_list%size)
          call lanc_build_gf_nonsu2_diagOrb_diagSpin(iorb,ispin)
       enddo
    enddo
    !
    !
    !same orbital, different spin GF: G_{aa}^{ss'}(z)
    do ispin=1,Nspin
       do jspin=1,Nspin
          if(ispin==jspin)cycle
          do iorb=1,Norb
             MaskBool=.true.   
             if(bath_type=="replica".or.bath_type=="general")MaskBool=Hmask(ispin,jspin,iorb,iorb)
             if(.not.MaskBool)cycle
             call allocate_GFmatrix(impGmatrix(ispin,jspin,iorb,iorb),Nstate=state_list%size)
             call lanc_build_gf_nonsu2_mixOrb_mixSpin(iorb,iorb,ispin,jspin)
          enddo
       enddo
    enddo
    do ispin=1,Nspin
       do jspin=1,Nspin
          if(ispin==jspin)cycle
          do iorb=1,Norb
             !
             impGmats(ispin,jspin,iorb,iorb,:) = 0.5d0*(impGmats(ispin,jspin,iorb,iorb,:) &
                  - (one-xi)*impGmats(ispin,ispin,iorb,iorb,:) &
                  - (one-xi)*impGmats(jspin,jspin,iorb,iorb,:))
             !
             impGreal(ispin,jspin,iorb,iorb,:) = 0.5d0*(impGreal(ispin,jspin,iorb,iorb,:) &
                  - (one-xi)*impGreal(ispin,ispin,iorb,iorb,:) &
                  - (one-xi)*impGreal(jspin,jspin,iorb,iorb,:))
             !
          enddo
       enddo
    enddo
    !
    !
    !
    !different orbital, same spin GF: G_{ab}^{ss}(z)
    select case(bath_type)
    case default;
    case("hybrid","replica","general")
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                if(iorb==jorb)cycle
                MaskBool=.true.   
                if(bath_type=="replica".or.bath_type=="general")MaskBool=Hmask(ispin,ispin,iorb,jorb)
                if(.not.MaskBool)cycle
                !
                call allocate_GFmatrix(impGmatrix(ispin,ispin,iorb,jorb),Nstate=state_list%size)
                call lanc_build_gf_nonsu2_mixOrb_mixSpin(iorb,jorb,ispin,ispin)
             enddo
          enddo
       enddo
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                if(iorb==jorb)cycle
                impGmats(ispin,ispin,iorb,jorb,:) = 0.5d0*(impGmats(ispin,ispin,iorb,jorb,:) &
                     - (one-xi)*impGmats(ispin,ispin,iorb,iorb,:) &
                     - (one-xi)*impGmats(ispin,ispin,jorb,jorb,:))
                !
                impGreal(ispin,ispin,iorb,jorb,:) = 0.5d0*(impGreal(ispin,ispin,iorb,jorb,:) &
                     - (one-xi)*impGreal(ispin,ispin,iorb,iorb,:) &
                     - (one-xi)*impGreal(ispin,ispin,jorb,jorb,:))
                !
             enddo
          enddo
       enddo
       !
       !
       !
       !different orbital, different spin GF: G_{ab}^{ss'}(z)
       do ispin=1,Nspin
          do jspin=1,Nspin
             if(ispin==jspin)cycle
             do iorb=1,Norb
                do jorb=1,Norb
                   if(iorb==jorb)cycle
                   MaskBool=.true.   
                   if(bath_type=="replica".or.bath_type=="general")MaskBool=Hmask(ispin,jspin,iorb,jorb)
                   if(.not.MaskBool)cycle
                   !
                   call allocate_GFmatrix(impGmatrix(ispin,jspin,iorb,jorb),Nstate=state_list%size)
                   call lanc_build_gf_nonsu2_mixOrb_mixSpin(iorb,jorb,ispin,jspin)
                enddo
             enddo
          enddo
       enddo
       do ispin=1,Nspin
          do jspin=1,Nspin
             if(ispin==jspin)cycle
             do iorb=1,Norb
                do jorb=1,Norb
                   if(iorb==jorb)cycle
                   impGmats(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGmats(ispin,jspin,iorb,jorb,:) &
                        - (one-xi)*impGmats(ispin,ispin,iorb,iorb,:) &
                        - (one-xi)*impGmats(jspin,jspin,jorb,jorb,:))
                   !
                   impGreal(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGreal(ispin,jspin,iorb,jorb,:) &
                        - (one-xi)*impGreal(ispin,ispin,iorb,iorb,:) &
                        - (one-xi)*impGreal(jspin,jspin,jorb,jorb,:))
                enddo
             enddo
          enddo
       enddo
    end select
    !
  end subroutine build_gf_nonsu2













  !PURPOSE: Evaluate the same orbital IORB, same spin ISPIN impurity GF.
  subroutine lanc_build_gf_nonsu2_diagOrb_diagSpin(iorb,ispin)
    integer      :: iorb,ispin,isite
    type(sector) :: sectorI,sectorJ
    integer :: ib(2*Ns)
    !
    write(LOGfile,"(A)")"Get G_l"//str(iorb)//str(iorb)//"_s"//str(ispin)//str(ispin)
    !
    do istate=1,state_list%size
       !
       call allocate_GFmatrix(impGmatrix(ispin,ispin,iorb,iorb),istate,Nchan=2) !2*[c,cdg]
       !
       isector    =  es_return_sector(state_list,istate)
       e_state    =  es_return_energy(state_list,istate)
       v_state    =  es_return_cvec(state_list,istate)
       !
       !EVALUATE c^+_{a,s}|v> --> G^>_{s,s;a,a}
       jsector = getCDGsector(1,ispin,isector)
       if(Jz_basis)jsector = getCDGsector_Jz(iorb,ispin,isector)
       if(getN(isector)/=Nlevels.and.jsector>=0)then
          vvinit = apply_op_CDG(v_state,iorb,ispin,isector,jsector)
          call tridiag_Hv_sector_nonsu2(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_nonsu2(one*norm2,e_state,alfa_,beta_,1,iorb,iorb,ispin,ispin,1,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(ispin,ispin,iorb,iorb),istate,1,Nexc=0)
       endif
       !
       !EVALUATE c_{a,s}|v> --> G^<_{s,s;a,a}
       jsector = getCsector(1,ispin,isector)
       if(Jz_basis)jsector = getCsector_Jz(iorb,ispin,isector)
       if(getN(isector)/=0.and.jsector>=0)then
          vvinit =  apply_op_C(v_state,iorb,ispin,isector,jsector)
          call tridiag_Hv_sector_nonsu2(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_nonsu2(one*norm2,e_state,alfa_,beta_,-1,iorb,iorb,ispin,ispin,2,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(ispin,ispin,iorb,iorb),istate,2,Nexc=0)
       endif
       !
       if(allocated(v_state))deallocate(v_state)
       !
    enddo
    return
  end subroutine lanc_build_gf_nonsu2_diagOrb_diagSpin




  !PURPOSE: Evaluate the  different orbital IORB,JORB, different spin ISPIN,JSPIN impurity GF.
  subroutine lanc_build_gf_nonsu2_mixOrb_mixSpin(iorb,jorb,ispin,jspin)
    integer      :: iorb,jorb,ispin,jspin
    type(sector) :: sectorI,sectorJ
    !
    write(LOGfile,"(A)")"Get G_l"//str(iorb)//str(iorb)//"_s"//str(ispin)//str(jspin)
    !
    do istate=1,state_list%size
       !
       call allocate_GFmatrix(impGmatrix(ispin,jspin,iorb,jorb),istate,Nchan=4) !2*[aux1,aux2]
       !
       isector    =  es_return_sector(state_list,istate)
       e_state    =  es_return_energy(state_list,istate)
       v_state    =  es_return_cvec(state_list,istate)
       !
       !
       !APPLY (c^+_{jorb,jspin} + c^+_{iorb,ispin})|gs> = [1,1].[C_{+1},C_{+1}].[iorb,jorb].[ispin,jspin]
       jsector = getCDGsector(1,ispin,isector)
       if(Jz_basis)then
          jsector = getCDGsector_Jz(iorb,ispin,isector)
          ksector = getCDGsector_Jz(jorb,jspin,isector)
          if(getDim(jsector)/=getDim(ksector))stop "lanczos builgf dimensional error"
       endif
       if(getN(isector)/=Nlevels.and.jsector>=0)then
          vvinit = apply_Cops(v_state,[one,one],[1,1],[iorb,jorb],[ispin,jspin],isector,jsector)
          call tridiag_Hv_sector_nonsu2(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_nonsu2(one*norm2,e_state,alfa_,beta_,1,iorb,jorb,ispin,jspin,1,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(ispin,jspin,iorb,jorb),istate,1,Nexc=0)
       endif
       !
       !
       !APPLY (c_{jorb,jspin} + c_{iorb,ispin})|gs> =  [1,1].[C_{-1},C_{-1}].[iorb,jorb].[ispin,jspin]
       jsector = getCsector(1,ispin,isector)
       if(Jz_basis)then
          jsector = getCsector_Jz(iorb,ispin,isector)
          ksector = getCsector_Jz(jorb,jspin,isector)
          if(getDim(jsector)/=getDim(ksector))stop "lanczos builgf dimensional error"
       endif
       if(getN(isector)/=0.and.jsector>=0)then
          vvinit = apply_Cops(v_state,[one,one],[-1,-1],[iorb,jorb],[ispin,jspin],isector,jsector)
          call tridiag_Hv_sector_nonsu2(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_nonsu2(one*norm2,e_state,alfa_,beta_,-1,iorb,jorb,ispin,jspin,2,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(ispin,jspin,iorb,jorb),istate,2,Nexc=0)
       endif
       !
       !
       !APPLY (+i*c^+_{jorb,jspin} + c^+_{iorb,ispin})|gs> =  [1,1j].[C_{+1},C_{+1}].[iorb,jorb].[ispin,jspin]
       jsector = getCDGsector(1,ispin,isector)
       if(Jz_basis)then
          jsector = getCDGsector_Jz(iorb,ispin,isector)
          ksector = getCDGsector_Jz(jorb,jspin,isector)
          if(getDim(jsector)/=getDim(ksector))stop "lanczos builgf dimensional error"
       endif
       !
       if(getN(isector)/=Nlevels.and.jsector>=0)then
          vvinit = apply_Cops(v_state,[one,xi],[1,1],[iorb,jorb],[ispin,jspin],isector,jsector)          
          call tridiag_Hv_sector_nonsu2(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_nonsu2(-xi*norm2,e_state,alfa_,beta_,1,iorb,jorb,ispin,jspin,3,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(ispin,jspin,iorb,jorb),istate,3,Nexc=0)
       endif
       !
       !
       !APPLY (-xi*c_{jorb,jspin} + c_{iorb,ispin})|gs> =  [1,-1j].[C_{-1},C_{-1}].[iorb,jorb].[ispin,jspin]
       jsector = getCsector(1,ispin,isector)
       if(Jz_basis)then
          jsector = getCsector_Jz(iorb,ispin,isector)
          ksector = getCsector_Jz(jorb,jspin,isector)
          if(getDim(jsector)/=getDim(ksector))stop "lanczos builgf dimensional error"
       endif
       if(getN(isector)/=0.and.jsector>=0)then
          vvinit = apply_Cops(v_state,[one,-xi],[-1,-1],[iorb,jorb],[ispin,jspin],isector,jsector)
          call tridiag_Hv_sector_nonsu2(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_nonsu2(-xi*norm2,e_state,alfa_,beta_,-1,iorb,jorb,ispin,jspin,4,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(ispin,jspin,iorb,jorb),istate,4,Nexc=0)
       endif
       !
       if(allocated(v_state))deallocate(v_state)
       !
    enddo
    return
  end subroutine lanc_build_gf_nonsu2_mixOrb_mixSpin






  subroutine add_to_lanczos_gf_nonsu2(vnorm2,Ei,alanc,blanc,isign,iorb,jorb,ispin,jspin,ichan,istate)
    complex(8)                                 :: vnorm2,pesoBZ,peso
    real(8)                                    :: Ei,Egs,de
    integer                                    :: nlanc,itype
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    integer                                    :: isign,iorb,jorb,ispin,jspin,ichan,istate
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
    complex(8)                                 :: iw
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")&
         "DEBUG add_to_lanczos_GF NONSU2: add-up to GF. istate:"//str(istate)
#endif
    !
    Egs = state_list%emin       !get the gs energy
    !
    Nlanc = size(alanc)
    !
    if((finiteT).and.(beta*(Ei-Egs).lt.200))then
       pesoBZ = vnorm2*exp(-beta*(Ei-Egs))/zeta_function
    elseif(.not.finiteT)then
       pesoBZ = vnorm2/zeta_function
    else
       pesoBZ=0.d0
    endif
    !pesoBZ = vnorm2/zeta_function
    !if(finiteT)pesoBZ = vnorm2*exp(-beta*(Ei-Egs))/zeta_function
    !
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,alanc)
       call Bcast_MPI(MpiComm,blanc)
    endif
#endif
    !
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A)")&
         "DEBUG add_to_lanczos_GF NONSU2: LApack tridiagonalization"
#endif
    call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
    !
    call allocate_GFmatrix(impGmatrix(ispin,jspin,iorb,jorb),istate,ichan,Nlanc)
    !    
    do j=1,nlanc
       de = diag(j)-Ei
       peso = pesoBZ*Z(1,j)*Z(1,j)
       !
       impGmatrix(ispin,jspin,iorb,jorb)%state(istate)%channel(ichan)%weight(j) = peso
       impGmatrix(ispin,jspin,iorb,jorb)%state(istate)%channel(ichan)%poles(j)  = isign*de
       !       
       do i=1,Lmats
          iw=xi*wm(i)
          impGmats(ispin,jspin,iorb,jorb,i)=impGmats(ispin,jspin,iorb,jorb,i) + peso/(iw-isign*de)
       enddo
       do i=1,Lreal
          iw=dcmplx(wr(i),eps)
          impGreal(ispin,jspin,iorb,jorb,i)=impGreal(ispin,jspin,iorb,jorb,i) + peso/(iw-isign*de)
       enddo
    enddo
  end subroutine add_to_lanczos_gf_nonsu2




  !################################################################
  !################################################################
  !################################################################
  !################################################################




  !+------------------------------------------------------------------+
  !PURPOSE  : Build the Self-energy functions, NONSU2 case
  !+------------------------------------------------------------------+
  subroutine build_sigma_nonsu2
    !
    ! Obtains the self-energy function :math:`\Sigma` on the current Matsubara and Real-axis intervals using impurity Dyson equation  :math:`\hat{\Sigma}(z) = \hat{G}^{-1}_0(z) - \hat{G}^{-1}(z)`
    !
    !
    integer                                           :: i,j,isign,unit(7),iorb,jorb,ispin,jspin,io,jo
    complex(8)                                        :: fg0
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: invG0mats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal) :: invG0real
    complex(8),dimension(Nspin*Norb,Nspin*Norb)       :: invGimp,invG0imp
    character(len=20)                                 :: suffix
    !
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")&
         "DEBUG build_sigma NONSU2: get Self-energy"
#endif
    !
    impG0mats = zero
    impG0real = zero
    invG0mats = zero
    invG0real = zero
    impSmats  = zero
    impSreal  = zero
    !
    !
    !Get G0^-1
    invG0mats = invg0_bath_function(dcmplx(0d0,wm(:)),dmft_bath)
    invG0real = invg0_bath_function(dcmplx(wr(:),eps),dmft_bath)
    !
    !Get Gimp^-1 - Matsubara freq.
    !Get Gimp^-1 - Real freq.
    do i=1,Lmats
       invGimp  = nn2so_reshape(impGmats(:,:,:,:,i),Nspin,Norb)
       call inv(invGimp)
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   impSmats(ispin,jspin,iorb,jorb,i) = invG0mats(ispin,jspin,iorb,jorb,i) - invGimp(io,jo)
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    do i=1,Lreal
       invGimp  = nn2so_reshape(impGreal(:,:,:,:,i),Nspin,Norb)
       call inv(invGimp)
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   impSreal(ispin,jspin,iorb,jorb,i) = invG0real(ispin,jspin,iorb,jorb,i) - invGimp(io,jo)
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    !Get G0and:
    impG0mats(:,:,:,:,:) = g0and_bath_function(dcmplx(0d0,wm(:)),dmft_bath)
    impG0real(:,:,:,:,:) = g0and_bath_function(dcmplx(wr(:),eps),dmft_bath)
    !
  end subroutine build_sigma_nonsu2





  !################################################################
  !################################################################
  !################################################################
  !################################################################


  subroutine rebuild_gf_nonsu2()
    !
    ! Reconstructs the system impurity electrons Green's functions using :f:var:`impgmatrix` to retrieve weights and poles.
    !
    integer                                     :: iorb,jorb,ispin,jspin,i,io,jo
    logical                                     :: MaskBool
    logical(8),dimension(Nspin,Nspin,Norb,Norb) :: Hmask
    !
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG rebuild_gf NONSU2: rebuild GFs"
#endif
    !
    Hmask= .true.
    if(.not.ed_all_g)then
       select case(bath_type)
       case("replica");Hmask=Hreplica_mask(wdiag=.true.,uplo=.false.)
       case("general");Hmask=Hgeneral_mask(wdiag=.true.,uplo=.false.)
       case default;stop "ERROR: ED_ALL_G=FALSE AND BATH_TYPE!=REPLICA/GENERAL"
       end select
    endif
    do ispin=1,Nspin
       do iorb=1,Norb
          write(LOGfile,*)((Hmask(ispin,jspin,iorb,jorb),jorb=1,Norb),jspin=1,Nspin)
       enddo
    enddo
    !
    !same orbital, same spin GF: G_{aa}^{ss}(z)
    do ispin=1,Nspin
       do iorb=1,Norb
          call rebuild_gf_nonsu2_main(iorb,iorb,ispin,ispin)
       enddo
    enddo
    !
    !
    !same orbital, different spin GF: G_{aa}^{ss'}(z)
    do ispin=1,Nspin
       do jspin=1,Nspin
          if(ispin==jspin)cycle
          do iorb=1,Norb
             MaskBool=.true.   
             if(bath_type=="replica".or.bath_type=="general")MaskBool=Hmask(ispin,jspin,iorb,iorb)
             if(.not.MaskBool)cycle
             !
             call rebuild_gf_nonsu2_main(iorb,iorb,ispin,jspin)
          enddo
       enddo
    enddo
    do ispin=1,Nspin
       do jspin=1,Nspin
          if(ispin==jspin)cycle
          do iorb=1,Norb
             impGmats(ispin,jspin,iorb,iorb,:) = 0.5d0*(impGmats(ispin,jspin,iorb,iorb,:) &
                  - (one-xi)*impGmats(ispin,ispin,iorb,iorb,:) &
                  - (one-xi)*impGmats(jspin,jspin,iorb,iorb,:))
             !
             impGreal(ispin,jspin,iorb,iorb,:) = 0.5d0*(impGreal(ispin,jspin,iorb,iorb,:) &
                  - (one-xi)*impGreal(ispin,ispin,iorb,iorb,:) &
                  - (one-xi)*impGreal(jspin,jspin,iorb,iorb,:))
             !
          enddo
       enddo
    enddo
    !
    !
    !
    select case(bath_type)
    case default;
    case("hybrid","replica","general")
       !different orbital, same spin GF: G_{ab}^{ss}(z)
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                if(iorb==jorb)cycle
                MaskBool=.true.   
                if(bath_type=="replica".or.bath_type=="general")MaskBool=Hmask(ispin,ispin,iorb,jorb)
                if(.not.MaskBool)cycle
                !
                call rebuild_gf_nonsu2_main(iorb,jorb,ispin,ispin)
             enddo
          enddo
       enddo
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                if(iorb==jorb)cycle
                impGmats(ispin,ispin,iorb,jorb,:) = 0.5d0*(impGmats(ispin,ispin,iorb,jorb,:) &
                     - (one-xi)*impGmats(ispin,ispin,iorb,iorb,:) &
                     - (one-xi)*impGmats(ispin,ispin,jorb,jorb,:))
                !
                impGreal(ispin,ispin,iorb,jorb,:) = 0.5d0*(impGreal(ispin,ispin,iorb,jorb,:) &
                     - (one-xi)*impGreal(ispin,ispin,iorb,iorb,:) &
                     - (one-xi)*impGreal(ispin,ispin,jorb,jorb,:))
                !
             enddo
          enddo
       enddo
       !
       !
       !
       !different orbital, different spin GF: G_{ab}^{ss'}(z)
       do ispin=1,Nspin
          do jspin=1,Nspin
             if(ispin==jspin)cycle
             do iorb=1,Norb
                do jorb=1,Norb
                   if(iorb==jorb)cycle
                   MaskBool=.true.   
                   if(bath_type=="replica".or.bath_type=="general")MaskBool=Hmask(ispin,jspin,iorb,jorb)
                   if(.not.MaskBool)cycle
                   !
                   call rebuild_gf_nonsu2_main(iorb,jorb,ispin,jspin)
                enddo
             enddo
          enddo
       enddo
       do ispin=1,Nspin
          do jspin=1,Nspin
             if(ispin==jspin)cycle
             do iorb=1,Norb
                do jorb=1,Norb
                   if(iorb==jorb)cycle
                   impGmats(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGmats(ispin,jspin,iorb,jorb,:) &
                        - (one-xi)*impGmats(ispin,ispin,iorb,iorb,:) &
                        - (one-xi)*impGmats(jspin,jspin,jorb,jorb,:))
                   !
                   impGreal(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGreal(ispin,jspin,iorb,jorb,:) &
                        - (one-xi)*impGreal(ispin,ispin,iorb,iorb,:) &
                        - (one-xi)*impGreal(jspin,jspin,jorb,jorb,:))
                enddo
             enddo
          enddo
       enddo
    end select
    !
  end subroutine rebuild_gf_nonsu2



  subroutine rebuild_gf_nonsu2_main(iorb,jorb,ispin,jspin)
    integer,intent(in) :: iorb,jorb,ispin,jspin
    integer            :: Nstates,istate
    integer            :: Nchannels,ichan
    integer            :: Nexcs,iexc
    complex(8)         :: peso
    real(8)            :: de
    !
    write(LOGfile,"(A)")"Get G_l"//str(iorb)//str(jorb)//"_s"//str(ispin)//str(jspin)
    !
    if(.not.allocated(impGmatrix(ispin,jspin,iorb,jorb)%state)) then
       print*, "ED_GF_NONSU2 WARNING: impGmatrix%state not allocated. Nothing to do"
       return
    endif
    !
    Nstates = size(impGmatrix(ispin,jspin,iorb,jorb)%state)
    do istate=1,Nstates
       Nchannels = size(impGmatrix(ispin,jspin,iorb,jorb)%state(istate)%channel)
       do ichan=1,Nchannels
          Nexcs  = size(impGmatrix(ispin,jspin,iorb,jorb)%state(istate)%channel(ichan)%poles)
          if(Nexcs==0)cycle
          do iexc=1,Nexcs
             peso  = impGmatrix(ispin,jspin,iorb,jorb)%state(istate)%channel(ichan)%weight(iexc)
             de    = impGmatrix(ispin,jspin,iorb,jorb)%state(istate)%channel(ichan)%poles(iexc)
             impGmats(ispin,jspin,iorb,jorb,:)=impGmats(ispin,jspin,iorb,jorb,:) + &
                  peso/(dcmplx(0d0,wm(i))-de)
             impGreal(ispin,jspin,iorb,jorb,:)=impGreal(ispin,jspin,iorb,jorb,:) + &
                  peso/(dcmplx(wr(i),eps)-de)
          enddo
       enddo
    enddo
    return
  end subroutine rebuild_gf_nonsu2_main











END MODULE ED_GF_NONSU2
