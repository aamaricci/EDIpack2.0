MODULE ED_GF_SUPERC
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER
  USE SF_IOTOOLS, only: str,reg,txtfy      , save_array
  USE SF_LINALG,  only: inv,eigh,eye
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_EIGENSPACE
  USE ED_BATH
  USE ED_SETUP
  USE ED_SECTOR
  USE ED_HAMILTONIAN_SUPERC
  implicit none
  private

  public :: build_gf_superc
  public :: rebuild_gf_superc
  public :: build_sigma_superc

  integer                               :: istate
  integer                               :: isector,jsector,ksector
  integer                               :: idim,jdim
  integer                               :: isz,jsz
  integer                               :: ialfa,jalfa
  integer                               :: i,j,m,r
  integer                               :: iph,i_el
  real(8)                               :: sgn,norm2
  !
  complex(8),allocatable                :: vvinit(:)
  real(8),allocatable                   :: alfa_(:),beta_(:)

  !Lanczos shared variables
  !=========================================================
  complex(8),dimension(:),allocatable   :: v_state
  real(8)                               :: e_state

  !AUX GF
  !=========================================================
  complex(8),allocatable,dimension(:,:) :: auxGmats,auxGreal

contains



  !+------------------------------------------------------------------+
  !                        SUPERC
  !+------------------------------------------------------------------+
  subroutine build_gf_superc()
    !
    !
    !Evaluates the impurity electrons Green's functions :math:`G(z)` and :math:`F(z)` and the phonons one :math:`D(z)` using dynamical Lanczos method. The result is stored in rank-5 arrays :f:var:`impgmats`, :f:var:`impgreal` , :f:var:`impfmats` , :f:var:`impfreal` of dimensions [ |Nspin| , |Nspin| , |Norb| , |Norb| , :f:var:`Lmats` / :f:var:`Lreal` ] and rank-1 array :f:var:`impdmats`, :f:var:`impdreal`.    
    !
    !The off-diagonal components of :math:`G_{ab}` with :math:`a \neq b` as well as the anomalous Green's functions :math:`F_{ab}(z)\, \forall a,b` are obtained using algebraic manipulation to ensure working with hermitian conjugate operators in the dynamical Lanczos procedure.  
    !
    !The weights and the poles obtained in this procedure are saved in a hierarchical data structure (for every state, every channel (creation or annihilation of excitations, normal or anomalous) and every degree of freedom) :f:var:`impgmatrix` of type :f:var:`gfmatrix`. 
    !
    ! .. _j.cpc.2021.108261: https://doi.org/10.1016/j.cpc.2021.108261
    !
    !
    integer    :: iorb,jorb,ispin,i,isign
    complex(8) :: barGmats(Norb,Lmats),barGreal(Norb,Lreal)
    !
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG build_gf SUPERC: build GFs"
#endif
    !
    if(.not.allocated(auxGmats))allocate(auxGmats(4,Lmats))
    if(.not.allocated(auxGreal))allocate(auxGreal(4,Lreal))
    auxgmats=zero
    auxGreal=zero
    barGmats=zero
    barGreal=zero
    !
    ispin=1                       !in this channel Nspin=2 is forbidden. check in ED_AUX_FUNX.
    !
    !1: Get G^{aa}_{upup}-> saved in impG(1,1,a,a) and \bar{G}^{aa}_{dwdw}-> saved in barG
    ! (ispin,ispin) to components:
    ! (1,1) -> diag normal (diag e offdiag)
    ! (1,2) -> anomalous (diag e offdiag)
    ! (2,2) -> bar       
    if(MPIMASTER)call start_timer(unit=LOGfile)    
    do iorb=1,Norb
       auxGmats=zero
       auxGreal=zero
       call allocate_GFmatrix(impGmatrix(1,1,iorb,iorb),Nstate=state_list%size)
       call allocate_GFmatrix(impGmatrix(2,2,iorb,iorb),Nstate=state_list%size)
       call lanc_build_gf_superc_Gdiag(iorb)
       impGmats(ispin,ispin,iorb,iorb,:) = auxGmats(1,:) !this is G_{up,up;iorb,iorb} <CCdg>
       impGreal(ispin,ispin,iorb,iorb,:) = auxGreal(1,:)  
       barGmats(                 iorb,:) = auxGmats(2,:) !this is \bar{G}_{dw,dw;iorb,iorb} <CdgC>
       barGreal(                 iorb,:) = auxGreal(2,:)
    enddo
    if(MPIMASTER)call stop_timer
    !
    !
    if(MPIMASTER)call start_timer(unit=LOGfile)
    select case(bath_type)
    case default
       !get F^{aa}_{updw} = <adg_up . adg_dw>
       do iorb=1,Norb
          auxGmats=zero
          auxGreal=zero
          call allocate_GFmatrix(impGmatrix(1,2,iorb,iorb),Nstate=state_list%size)
          call lanc_build_gf_superc_Fmix(iorb,iorb)!<=auxG(4,:)= <O.Odg> -xi*<P.Pdg>, w O=(a_up+adg_dw), P=(a_up+xi*adg_dw)
          impFmats(ispin,ispin,iorb,iorb,:) = 0.5d0*(auxGmats(4,:)-(one-xi)*(impGmats(ispin,ispin,iorb,iorb,:)+barGmats(iorb,:)) )
          impFreal(ispin,ispin,iorb,iorb,:) = 0.5d0*(auxGreal(4,:)-(one-xi)*(impGreal(ispin,ispin,iorb,iorb,:)+barGreal(iorb,:)) )
       enddo
       !
    case ('hybrid','replica','general')
       !Get G^{ab}_{upup} --> saved in impG(1,1,a,b)
       do iorb=1,Norb
          do jorb=1,Norb
             if(iorb==jorb)cycle
             auxGmats=zero
             auxGreal=zero
             call allocate_GFmatrix(impGmatrix(1,1,iorb,jorb),Nstate=state_list%size)
             call lanc_build_gf_superc_Gmix(iorb,jorb)!<= auxG(3,:)= <Q.Qdg> -xi*<R.Rdg>, w Q=(a_up+b_up), R=(a_up+xi.b_up)
             impGmats(ispin,ispin,iorb,jorb,:) = 0.5d0*(auxGmats(3,:)-(one-xi)*(impGmats(ispin,ispin,iorb,iorb,:)+impGmats(ispin,ispin,jorb,jorb,:)))
             impGreal(ispin,ispin,iorb,jorb,:) = 0.5d0*(auxGreal(3,:)-(one-xi)*(impGreal(ispin,ispin,iorb,iorb,:)+impGreal(ispin,ispin,jorb,jorb,:)))
          enddo
       enddo
       !
       !get F^{ab}_{updw} = <adg_up . bdg_dw>
       do iorb=1,Norb
          do jorb=1,Norb
             auxGmats=zero
             auxGreal=zero
             call allocate_GFmatrix(impGmatrix(1,2,iorb,jorb),Nstate=state_list%size)
             call lanc_build_gf_superc_Fmix(iorb,jorb)!<=auxG(4,:)= <O.Odg> -xi*<P.Pdg>, w O=(a_up+bdg_dw), P=(a_up+xi*bdg_dw)
             impFmats(ispin,ispin,iorb,jorb,:) = 0.5d0*(auxGmats(4,:)-(one-xi)*(impGmats(ispin,ispin,iorb,iorb,:)+barGmats(jorb,:)) )
             impFreal(ispin,ispin,iorb,jorb,:) = 0.5d0*(auxGreal(4,:)-(one-xi)*(impGreal(ispin,ispin,iorb,iorb,:)+barGreal(jorb,:)) )
          enddo
       enddo
    end select
    !
    !PHONONS
    if(DimPh>1)then
       call lanc_build_gf_phonon_main()
    endif
    !
    if(MPIMASTER)call stop_timer
    !
    !
    deallocate(auxGmats,auxGreal)
    !
  end subroutine build_gf_superc















  subroutine lanc_build_gf_superc_Gdiag(iorb)
    integer      :: iorb
    !
    write(LOGfile,"(A)")"Get G & barG_l"//str(iorb)
    !
    do istate=1,state_list%size
       !
       call allocate_GFmatrix(impGmatrix(1,1,iorb,iorb),istate,Nchan=2) ![c+_up,c_up]
       call allocate_GFmatrix(impGmatrix(2,2,iorb,iorb),istate,Nchan=2) ![c_dw,c+_dw]
       !
       isector    =  es_return_sector(state_list,istate)
       e_state    =  es_return_energy(state_list,istate)
       v_state    =  es_return_cvec(state_list,istate)
       !
       !EVALUATE c^+_{up,iorb}|v> --> Gaux(1) = G^>_{up,up;iorb,iorb}
       jsector = getCDGsector(1,1,isector)
       if(jsector/=0)then 
          vvinit = apply_op_CDG(v_state,iorb,1,isector,jsector)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,e_state,alfa_,beta_,1,1,iorb,iorb,1,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,iorb,iorb),istate,1,Nexc=0)
       endif
       !
       !EVALUATE c_{up,iorb}|v> --> Gaux(1) = G^<_{up,up;iorb,iorb}
       jsector = getCsector(1,1,isector)
       if(jsector/=0)then
          vvinit = apply_op_C(v_state,iorb,1,isector,jsector)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,e_state,alfa_,beta_,-1,1,iorb,iorb,2,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,iorb,iorb),istate,2,Nexc=0)
       endif
       !
       !
       !
       !EVALUATE c_{dw,iorb}|v> --> Gaux(2) = G^>_{dw,dw;iorb,iorb}
       jsector = getCsector(1,2,isector)
       if(jsector/=0)then
          vvinit = apply_op_C(v_state,iorb,2,isector,jsector)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,e_state,alfa_,beta_,1,2,iorb,iorb,1,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(2,2,iorb,iorb),istate,1,Nexc=0)
       endif
       !
       !EVALUATE c^+_{dw,iorb}|v> --> Gaux(2) = G^<_{dw,dw;iorb,iorb}
       jsector = getCDGsector(1,2,isector)
       if(jsector/=0)then
          vvinit = apply_op_CDG(v_state,iorb,2,isector,jsector)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,e_state,alfa_,beta_,-1,2,iorb,iorb,2,istate)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       else
          call allocate_GFmatrix(impGmatrix(2,2,iorb,iorb),istate,2,Nexc=0)
       endif
       !
       if(allocated(v_state))deallocate(v_state)
       !
    enddo
    return
  end subroutine lanc_build_gf_superc_Gdiag










  subroutine lanc_build_gf_superc_Gmix(iorb,jorb)
    integer      :: iorb,jorb
    !
    write(LOGfile,"(A)")"Get G_l"//str(iorb)//"_m"//str(jorb)
    !
    do istate=1,state_list%size
       !
       call allocate_GFmatrix(impGmatrix(1,1,iorb,jorb),istate,Nchan=4) !2*[QdgQ,RdgR]
       !
       isector    =  es_return_sector(state_list,istate)
       e_state    =  es_return_energy(state_list,istate)
       v_state    =  es_return_cvec(state_list,istate)
       !
       !EVALUATE (cdg_iorb + cdg_jorb)|gs> = Qdg|gs> += <Q.Qdg>
       ! =[1,1].[C_{+1},C_{+1}].[iorb,jorb].[up,up]
       jsector = getCDGsector(1,1,isector)
       if(jsector/=0)then
          vvinit = apply_Cops(v_state,[one,one],[1,1],[iorb,jorb],[1,1],isector,jsector)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,e_state,alfa_,beta_,1,3,iorb,jorb,1,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,iorb,jorb),istate,1,Nexc=0)
       endif
       !
       !EVALUATE (c_iorb + c_jorb)|gs> = Q|gs> += <Qdg.Q>
       !=[1,1].[C_{-1},C_{-1}].[iorb,jorb].[up,up]
       jsector = getCsector(1,1,isector)
       if(jsector/=0)then
          vvinit = apply_Cops(v_state,[one,one],[-1,-1],[iorb,jorb],[1,1],isector,jsector)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,e_state,alfa_,beta_,-1,3,iorb,jorb,2,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,iorb,jorb),istate,2,Nexc=0)
       endif
       !
       !
       !
       !EVALUATE (cdg_iorb + xi cdg_jorb)|gs> = Rdg|gs> += -xi*<Rdg.R>
       !=[1,1j].[C_{+1},C_{+1}].[iorb,jorb].[up,up]
       jsector = getCDGsector(1,1,isector)
       if(jsector/=0)then
          vvinit =  apply_Cops(v_state,[one,xi],[1,1],[iorb,jorb],[1,1],isector,jsector)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(-xi*norm2,e_state,alfa_,beta_,1,3,iorb,jorb,3,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,iorb,jorb),istate,3,Nexc=0)
       endif
       !
       !EVALUATE (c_iorb - xi c_jorb)|gs>  = R|gs> += -xi*<R.Rdg>
       !=[1,-1j].[C_{-1},C_{-1}].[iorb,jorb].[up,up]
       jsector = getCsector(1,1,isector)
       if(jsector/=0)then
          vvinit =  apply_Cops(v_state,[one,-xi],[-1,-1],[iorb,jorb],[1,1],isector,jsector)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(-xi*norm2,e_state,alfa_,beta_,-1,3,iorb,jorb,4,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,iorb,jorb),istate,4,Nexc=0)
       endif
       !
       if(allocated(v_state))deallocate(v_state)
       !
    enddo
    return
  end subroutine lanc_build_gf_superc_Gmix







  subroutine lanc_build_gf_superc_Fmix(iorb,jorb)
    integer                :: iorb,jorb
    type(sector)           :: sectorI,sectorJ
    !
    write(LOGfile,"(A)")"Get F_l"//str(iorb)//"_m"//str(jorb)
    !
    do istate=1,state_list%size
       call allocate_GFmatrix(impGmatrix(1,2,iorb,jorb),istate,Nchan=4) !2*[Odg.O,Pdg.P]
       !
       isector    =  es_return_sector(state_list,istate)
       e_state    =  es_return_energy(state_list,istate)
       v_state    =  es_return_cvec(state_list,istate)
       !
       !EVALUATE [cdg_{up,iorb} + c_{dw,jorb}]|gs> = Odg|gs> += <O.Odg>
       ! =[1,1].[C_{+1},C_{-1}].[iorb,jorb].[up,dw]
       isz = getsz(isector)
       if(isz<Ns)then
          jsector = getsector(isz+1,1)
          vvinit = apply_Cops(v_state,[one,one],[1,-1],[iorb,jorb],[1,2],isector,jsector)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,e_state,alfa_,beta_,1,4,iorb,jorb,1,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,2,iorb,jorb),istate,1,Nexc=0)
       endif
       !
       !EVALUATE [c_{up,iorb} + cdg_{dw,jorb}]|gs> = O|gs> += <Odg.O>
       ! =[1,1].[C_{-1},C_{+1}].[iorb,jorb].[up,dw]
       isz = getsz(isector)
       if(isz>-Ns)then
          jsector = getsector(isz-1,1)
          vvinit = apply_Cops(v_state,[one,one],[-1,1],[iorb,jorb],[1,2],isector,jsector)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,e_state,alfa_,beta_,-1,4,iorb,jorb,2,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,2,iorb,jorb),istate,2,Nexc=0)
       endif
       !
       !
       !EVALUATE [c^+_{up,iorb} + xi*c_{dw,jorb}]|gs> = Pdg|gs> += -xi*<P.Pdg>
       ! =[1,1j].[C_{+1},C_{-1}].[iorb,jorb].[up,dw]
       isz = getsz(isector)
       if(isz<Ns)then
          jsector = getsector(isz+1,1)
          vvinit = apply_Cops(v_state,[one,xi],[1,-1],[iorb,jorb],[1,2],isector,jsector)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(-xi*norm2,e_state,alfa_,beta_,1,4,iorb,jorb,3,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,2,iorb,jorb),istate,3,Nexc=0)
       endif
       !
       !EVALUATE [c_{up,iorb} - xi*c^+_{dw,jorb}]|gs> = P|gs> += -xi*<Pdg.P>
       ! =[1,-1j].[C_{-1},C_{+1}].[iorb,jorb].[up,dw]
       isz = getsz(isector)
       if(isz>-Ns)then
          jsector = getsector(isz-1,1)
          vvinit = apply_Cops(v_state,[one,-xi],[-1,+1],[iorb,jorb],[1,2],isector,jsector)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(-xi*norm2,e_state,alfa_,beta_,-1,4,iorb,jorb,4,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,2,iorb,jorb),istate,4,Nexc=0)
       endif
       !
       if(allocated(v_state))deallocate(v_state)
       !
    enddo
    !
    return
  end subroutine lanc_build_gf_superc_Fmix





  subroutine add_to_lanczos_gf_superc(vnorm2,Ei,alanc,blanc,isign,ichan,iorb,jorb,ic,istate)
    complex(8)                                 :: vnorm2,pesoBZ,peso
    real(8)                                    :: Ei,Egs,de
    integer                                    :: nlanc,itype
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    integer                                    :: isign,ichan,iorb,jorb,istate,ic
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr,chanI,chanJ
    complex(8)                                 :: iw
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")&
         "DEBUG add_to_lanczos_GF SUPERC: add-up to GF. istate:"//str(istate)
#endif
    !
    Egs = state_list%emin       !get the gs energy
    !
    Nlanc = size(alanc)
    !
    ! if((finiteT).and.(beta*(Ei-Egs).lt.200))then
    !    pesoBZ = vnorm2*exp(-beta*(Ei-Egs))/zeta_function
    ! elseif(.not.finiteT)then
    !    pesoBZ = vnorm2/zeta_function
    ! else
    !    pesoBZ=0.d0
    ! endif
    !
    pesoBZ = vnorm2/zeta_function
    if(finiteT)pesoBZ = vnorm2*exp(-beta*(Ei-Egs))/zeta_function
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
         "DEBUG add_to_lanczos_GF SUPERC: LApack tridiagonalization"
#endif
    call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
    !
    select case(ichan)
    case(1,3)
       chanI=1
       chanJ=1
    case(2)
       chanI=2
       chanJ=2
    case(4)
       chanI=1
       chanJ=2
    case default
       STOP "add_to_lanczos_gf_superc: wrong ichan"
    end select
    !
    call allocate_GFmatrix(impGmatrix(chanI,chanJ,iorb,jorb),istate,ic,Nlanc)
    !
    do j=1,nlanc
       de = diag(j)-Ei
       peso = pesoBZ*Z(1,j)*Z(1,j)
       !
       impGmatrix(chanI,chanJ,iorb,jorb)%state(istate)%channel(ic)%weight(j) = peso
       impGmatrix(chanI,chanJ,iorb,jorb)%state(istate)%channel(ic)%poles(j)  = isign*de
       !
       do i=1,Lmats
          iw=xi*wm(i)
          auxGmats(ichan,i)=auxGmats(ichan,i) + peso/(iw-isign*de)
       enddo
       do i=1,Lreal
          iw=dcmplx(wr(i),eps)
          auxGreal(ichan,i)=auxGreal(ichan,i) + peso/(iw-isign*de)
       enddo
    enddo
  end subroutine add_to_lanczos_gf_superc







  !################################################################
  !################################################################
  !################################################################
  !################################################################






  !+------------------------------------------------------------------+
  !PURPOSE  : Build the Self-energy functions, SUPERC case
  !+------------------------------------------------------------------+
  subroutine build_sigma_superc
    !
    ! Obtains the self-energy function :math:`\Sigma` on the current Matsubara and Real-axis intervals using impurity Dyson equation  :math:`\hat{\hat{\Sigma}}(z) = \hat{\hat{G}}^{-1}_0(z) - \hat{\hat{G}}^{-1}(z)`. 
    !
    !
    integer                                               :: i,ispin,iorb
    real(8)                                               :: det_mats(Lmats)
    complex(8)                                            :: det_real(Lreal)
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats)     :: invG0mats,invF0mats,invGmats,invFmats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal)     :: invG0real,invF0real,invGreal,invFreal
    complex(8),dimension(2*Nspin*Norb,2*Nspin*Norb)       :: invGimp
    !
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")&
         "DEBUG build_sigma SUPERC: get Self-energy"
#endif
    !
    invG0mats = zero
    invF0mats = zero
    invGmats  = zero
    invFmats  = zero
    invG0real = zero
    invF0real = zero
    invGreal  = zero
    invFreal  = zero
    impSmats  = zero
    impSAmats = zero
    impSreal  = zero
    impSAreal = zero
    !
    !Get G0^-1,F0^-1
    ispin=1
    invG0mats(:,:,:,:,:) = invg0_bath_function(dcmplx(0d0,wm(:)),dmft_bath)
    invF0mats(:,:,:,:,:) = invf0_bath_function(dcmplx(0d0,wm(:)),dmft_bath)
    invG0real(:,:,:,:,:) = invg0_bath_function(dcmplx(wr(:),eps),dmft_bath,axis='real')
    invF0real(:,:,:,:,:) = invf0_bath_function(dcmplx(wr(:),eps),dmft_bath,axis='real')
    !
    !Get Gimp^-1
    select case(bath_type)
       !
    case default
       do iorb=1,Norb
          det_mats  =  abs(impGmats(ispin,ispin,iorb,iorb,:))**2 + (impFmats(ispin,ispin,iorb,iorb,:))**2
          invGmats(ispin,ispin,iorb,iorb,:) = conjg(impGmats(ispin,ispin,iorb,iorb,:))/det_mats
          invFmats(ispin,ispin,iorb,iorb,:) = impFmats(ispin,ispin,iorb,iorb,:)/det_mats
          !
          det_real  = -impGreal(ispin,ispin,iorb,iorb,:)*conjg(impGreal(ispin,ispin,iorb,iorb,Lreal:1:-1)) - impFreal(ispin,ispin,iorb,iorb,:)**2
          invGreal(ispin,ispin,iorb,iorb,:) =  -conjg(impGreal(ispin,ispin,iorb,iorb,Lreal:1:-1))/det_real(:)
          invFreal(ispin,ispin,iorb,iorb,:) =  -impFreal(ispin,ispin,iorb,iorb,:)/det_real(:)
       enddo
       do iorb=1,Norb
          impSmats(ispin,ispin,iorb,iorb,:)  = invG0mats(ispin,ispin,iorb,iorb,:) - invGmats(ispin,ispin,iorb,iorb,:)
          impSAmats(ispin,ispin,iorb,iorb,:) = invF0mats(ispin,ispin,iorb,iorb,:) - invFmats(ispin,ispin,iorb,iorb,:)
          impSreal(ispin,ispin,iorb,iorb,:)  = invG0real(ispin,ispin,iorb,iorb,:) - invGreal(ispin,ispin,iorb,iorb,:)
          impSAreal(ispin,ispin,iorb,iorb,:) = invF0real(ispin,ispin,iorb,iorb,:) - invFreal(ispin,ispin,iorb,iorb,:)
       enddo
       !
       !
    case ("hybrid","replica","general")
       do i=1,Lmats
          invGimp=zero
          invGimp(1     :Norb  ,     1:Norb  ) = impGmats(ispin,ispin,:,:,i)
          invGimp(1     :Norb  ,Norb+1:2*Norb) = impFmats(ispin,ispin,:,:,i)
          invGimp(Norb+1:2*Norb,     1:Norb  ) = conjg(impFmats(ispin,ispin,:,:,i)) !this is real so conjg does none, but it shouldn't be there
          invGimp(Norb+1:2*Norb,Norb+1:2*Norb) =-conjg(impGmats(ispin,ispin,:,:,i))
          call inv(invGimp)
          invGmats(ispin,ispin,:,:,i) = invGimp(1:Norb,     1:Norb  )
          invFmats(ispin,ispin,:,:,i) = invGimp(1:Norb,Norb+1:2*Norb)
       enddo
       do i=1,Lreal
          invGimp=zero
          invGimp(1:Norb       ,     1:Norb)   = impGreal(ispin,ispin,:,:,i)
          invGimp(1:Norb       ,Norb+1:2*Norb) = impFreal(ispin,ispin,:,:,i)
          invGimp(Norb+1:2*Norb,     1:Norb)   = impFreal(ispin,ispin,:,:,i)
          invGimp(Norb+1:2*Norb,Norb+1:2*Norb) =-conjg(impGreal(ispin,ispin,:,:,Lreal-i+1))
          call inv(invGimp)
          invGreal(ispin,ispin,:,:,i) =  invGimp(1:Norb,     1:Norb  )
          invFreal(ispin,ispin,:,:,i) =  invGimp(1:Norb,Norb+1:2*Norb)
       enddo
       impSmats(ispin,ispin,:,:,:)  = invG0mats(ispin,ispin,:,:,:) - invGmats(ispin,ispin,:,:,:)
       impSAmats(ispin,ispin,:,:,:) = invF0mats(ispin,ispin,:,:,:) - invFmats(ispin,ispin,:,:,:)
       impSreal(ispin,ispin,:,:,:)  = invG0real(ispin,ispin,:,:,:) - invGreal(ispin,ispin,:,:,:)
       impSAreal(ispin,ispin,:,:,:) = invF0real(ispin,ispin,:,:,:) - invFreal(ispin,ispin,:,:,:)
       !
       !
    end select
    !
    !Get G0and:
    impG0mats(:,:,:,:,:) = g0and_bath_function(dcmplx(0d0,wm(:)),dmft_bath)
    impF0mats(:,:,:,:,:) = f0and_bath_function(dcmplx(0d0,wm(:)),dmft_bath)
    impG0real(:,:,:,:,:) = g0and_bath_function(dcmplx(wr(:),eps),dmft_bath,axis='real')
    impF0real(:,:,:,:,:) = f0and_bath_function(dcmplx(wr(:),eps),dmft_bath,axis='real')
    !
    !
  end subroutine build_sigma_superc




  !################################################################
  !################################################################
  !################################################################
  !################################################################




  subroutine rebuild_gf_superc()
    !
    ! Reconstructs the system impurity electrons Green's functions using :f:var:`impgmatrix` to retrieve weights and poles.
    !
    integer    :: iorb,jorb,ispin,i,isign
    complex(8) :: barGmats(Norb,Lmats),barGreal(Norb,Lreal)
    !
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG rebuild_gf SUPERC: rebuild GFs"
#endif
    !
    if(.not.allocated(auxGmats))allocate(auxGmats(4,Lmats))
    if(.not.allocated(auxGreal))allocate(auxGreal(4,Lreal))
    auxgmats=zero
    auxGreal=zero
    barGmats=zero
    barGreal=zero
    !
    ispin=1                       !in this channel Nspin=2 is forbidden. check in ED_AUX_FUNX.
    !
    do iorb=1,Norb
       auxGmats=zero
       auxGreal=zero
       write(LOGfile,"(A)")"Get G_l"//str(iorb) !//"_s"//str(ispin)
       call rebuild_gf_superc_Gdiag(iorb)
       !
       impGmats(ispin,ispin,iorb,iorb,:) = auxGmats(1,:) !this is G_{up,up;iorb,iorb}
       impGreal(ispin,ispin,iorb,iorb,:) = auxGreal(1,:)  
       barGmats(                 iorb,:) = auxGmats(2,:) !this is G_{dw,dw;iorb,iorb}
       barGreal(                 iorb,:) = auxGreal(2,:)
       !                                   auxGmats(3,:) !this is G_oo -xi*G_pp = A_aa -xi*B_aa
       !                                   auxGreal(3,:)
       !impFmats(ispin,ispin,iorb,iorb,:) = 0.5d0*(auxGmats(3,:)-(one-xi)*auxGmats(1,:)-(one-xi)*auxGmats(2,:))
       !impFreal(ispin,ispin,iorb,iorb,:) = 0.5d0*(auxGreal(3,:)-(one-xi)*auxGreal(1,:)-(one-xi)*auxGreal(2,:))
    enddo
    !
    !now we add the other mixed/anomalous GF in for the bath_type="hybrid" case
    select case(bath_type)
    case default
       ! Get F^{aa}_{updw} = <adg_up . adg_dw>
       do iorb=1,Norb
          write(LOGfile,"(A)")"Get F_l"//str(iorb)
          call rebuild_gf_superc_Fmix(iorb,iorb)
#ifdef _DEBUG
          write(Logfile,"(A)")""
#endif
          impFmats(ispin,ispin,iorb,iorb,:) = 0.5d0*(auxGmats(4,:)-(one-xi)*(impGmats(ispin,ispin,iorb,iorb,:)+barGmats(iorb,:)) )
          impFreal(ispin,ispin,iorb,iorb,:) = 0.5d0*(auxGreal(4,:)-(one-xi)*(impGreal(ispin,ispin,iorb,iorb,:)+barGreal(iorb,:)) )
       enddo
       !
    case ('hybrid','replica','general')
       !
       ! Get G^{ab}_{upup} = <adg_up . b_up>
       do iorb=1,Norb
          do jorb=1,Norb
             if(iorb==jorb)cycle
             call rebuild_gf_superc_Gmix(iorb,jorb)
             ! auxG{mats,real}(3,:) !this is <Qdg.Q> -xi*<Rdg.R>, w Q=(iorb+jorb), R=(iorb+xi.jorb)
             impGmats(ispin,ispin,iorb,jorb,:) = 0.5d0*(auxGmats(3,:)-(one-xi)*(impGmats(ispin,ispin,iorb,iorb,:)+impGmats(ispin,ispin,jorb,jorb,:)))
             impGreal(ispin,ispin,iorb,jorb,:) = 0.5d0*(auxGreal(3,:)-(one-xi)*(impGreal(ispin,ispin,iorb,iorb,:)+impGreal(ispin,ispin,jorb,jorb,:)))
          enddo
       enddo
       !
       ! Get F^{ab}_{updw} = <adg_up . bdg_dw>
       do iorb=1,Norb
          do jorb=1,Norb
             write(LOGfile,"(A)")"Get F_l"//str(iorb)//"_m"//str(jorb)//"_s"//str(ispin)
             call rebuild_gf_superc_Fmix(iorb,jorb)
             impFmats(ispin,ispin,iorb,jorb,:) = 0.5d0*( impGmats(ispin,ispin,iorb,jorb,:) - &
                  (one-xi)*impGmats(ispin,ispin,iorb,iorb,:) - (one-xi)*barGmats(jorb,:) )
             impFreal(ispin,ispin,iorb,jorb,:) = 0.5d0*( impGreal(ispin,ispin,iorb,jorb,:) - &
                  (one-xi)*impGreal(ispin,ispin,iorb,iorb,:) - (one-xi)*barGreal(jorb,:) )
          enddo
       enddo
    end select
    deallocate(auxGmats,auxGreal)
  end subroutine rebuild_gf_superc



  subroutine rebuild_gf_superc_Gdiag(iorb)
    integer,intent(in) :: iorb
    integer            :: Nstates,istate
    integer            :: Nchannels,ic,ichan
    integer            :: Nexcs,iexc
    complex(8)         :: peso
    real(8)            :: de
    !
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")&
         "DEBUG rebuild_gf SUPERC: reconstruct G_diagonal impurity GFs"
#endif
    !
    if(.not.allocated(impGmatrix(1,1,iorb,iorb)%state)) then
       print*, "ED_GF_SUPERC WARNING: impGmatrix%state not allocated. Nothing to do"
       return
    endif
    !
    Nstates = size(impGmatrix(1,1,iorb,iorb)%state)
    do istate=1,Nstates
       Nchannels = size(impGmatrix(1,1,iorb,iorb)%state(istate)%channel)     
       do ic=1,Nchannels        
          Nexcs  = size(impGmatrix(1,1,iorb,iorb)%state(istate)%channel(ic)%poles)
          if(Nexcs==0)cycle
          select case(ic)
          case(1,2);ichan=1
          case(3,4);ichan=2
          end select
          do iexc=1,Nexcs
             peso  = impGmatrix(1,1,iorb,iorb)%state(istate)%channel(ic)%weight(iexc)
             de    = impGmatrix(1,1,iorb,iorb)%state(istate)%channel(ic)%poles(iexc)
             auxGmats(ichan,:)=auxGmats(ichan,:) + peso/(dcmplx(0d0,wm(i))-de)
             auxGreal(ichan,:)=auxGreal(ichan,:) + peso/(dcmplx(wr(i),eps)-de)
          enddo
       enddo
    enddo
    return
  end subroutine rebuild_gf_superc_Gdiag

  subroutine rebuild_gf_superc_Gmix(iorb,jorb)
    integer,intent(in) :: iorb,jorb
    integer            :: Nstates,istate
    integer            :: Nchannels,ic,ichan
    integer            :: Nexcs,iexc
    complex(8)         :: peso
    real(8)            :: de
    !
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")&
         "DEBUG rebuild_gf SUPERC: reconstruct G_mixed impurity GFs"
#endif
    !
    !
    if(.not.allocated(impGmatrix(1,1,iorb,jorb)%state)) then
       print*, "ED_GF_SUPERC WARNING: impGmatrix%state not allocated. Nothing to do"
       return
    endif
    !
    ichan = 3
    !
    Nstates = size(impGmatrix(1,1,iorb,jorb)%state)
    do istate=1,Nstates
       Nchannels = size(impGmatrix(1,1,iorb,jorb)%state(istate)%channel)     
       do ic=1,Nchannels        
          Nexcs  = size(impGmatrix(1,1,iorb,jorb)%state(istate)%channel(ic)%poles)
          if(Nexcs==0)cycle
          do iexc=1,Nexcs
             peso  = impGmatrix(1,1,iorb,jorb)%state(istate)%channel(ic)%weight(iexc)
             de    = impGmatrix(1,1,iorb,jorb)%state(istate)%channel(ic)%poles(iexc)
             auxGmats(ichan,:)=auxGmats(ichan,:) + peso/(dcmplx(0d0,wm(i))-de)
             auxGreal(ichan,:)=auxGreal(ichan,:) + peso/(dcmplx(wr(i),eps)-de)
          enddo
       enddo
    enddo
    return
  end subroutine rebuild_gf_superc_Gmix


  subroutine rebuild_gf_superc_Fmix(iorb,jorb)
    integer,intent(in) :: iorb,jorb
    integer            :: Nstates,istate
    integer            :: Nchannels,ic,ichan
    integer            :: Nexcs,iexc
    complex(8)         :: peso
    real(8)            :: de
    !
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")&
         "DEBUG rebuild_gf SUPERC: reconstruct F_mixed impurity GFs"
#endif
    !
    !
    if(.not.allocated(impGmatrix(Nnambu,Nnambu,iorb,jorb)%state)) then
       print*, "ED_GF_SUPERC WARNING: impGmatrix%state not allocated. Nothing to do"
       return
    endif
    !
    ichan = 4
    !
    Nstates = size(impGmatrix(Nnambu,Nnambu,iorb,jorb)%state)
    do istate=1,Nstates
       Nchannels = size(impGmatrix(Nnambu,Nnambu,iorb,jorb)%state(istate)%channel)     
       do ic=1,Nchannels        
          Nexcs  = size(impGmatrix(Nnambu,Nnambu,iorb,jorb)%state(istate)%channel(ic)%poles)
          if(Nexcs==0)cycle
          do iexc=1,Nexcs
             peso  = impGmatrix(Nnambu,Nnambu,iorb,jorb)%state(istate)%channel(ic)%weight(iexc)
             de    = impGmatrix(Nnambu,Nnambu,iorb,jorb)%state(istate)%channel(ic)%poles(iexc)
             auxGmats(ichan,:)=auxGmats(ichan,:) + peso/(dcmplx(0d0,wm(i))-de)
             auxGreal(ichan,:)=auxGreal(ichan,:) + peso/(dcmplx(wr(i),eps)-de)
          enddo
       enddo
    enddo
    return
  end subroutine rebuild_gf_superc_Fmix





  !################################################################
  !################################################################
  !################################################################
  !################################################################








  subroutine lanc_build_gf_phonon_main()
    type(sector)                :: sectorI
    integer,dimension(Ns_Ud)    :: iDimUps,iDimDws
    integer                     :: Nups(Ns_Ud)
    integer                     :: Ndws(Ns_Ud)
    !
    write(LOGfile,"(A)")"Get phonon Green function:"
    do istate=1,state_list%size
       !
       ! call allocate_GFmatrix(impDmatrix,istate=istate,Nchan=1)
       !
       isector    =  es_return_sector(state_list,istate)
       e_state    =  es_return_energy(state_list,istate)
       v_state    =  es_return_cvec(state_list,istate)
       !
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          if(ed_verbose>=3)write(LOGfile,"(A,I6,I6)")&
               'From sector  :',isector,sectorI%Sz
       endif
       !
       if(MpiMaster)then
          if(ed_verbose>=3)write(LOGfile,"(A20,I12)")'Apply x',isector
          !
          allocate(vvinit(sectorI%Dim));vvinit=0d0
          !
          do i=1,sectorI%Dim
             iph = (i-1)/(sectorI%DimEl) + 1
             i_el = mod(i-1,sectorI%DimEl) + 1
             !
             !apply destruction operator
             if(iph>1) then
                j = i_el + ((iph-1)-1)*sectorI%DimEl
                vvinit(j) = vvinit(j) + sqrt(dble(iph-1))*v_state(i)
             endif
             !
             !apply creation operator
             if(iph<DimPh) then
                j = i_el + ((iph+1)-1)*sectorI%DimEl
                vvinit(j) = vvinit(j) + sqrt(dble(iph))*v_state(i)
             endif
          enddo
       else
          allocate(vvinit(1));vvinit=0.d0
       endif
       !
       call tridiag_Hv_sector_superc(isector,vvinit,alfa_,beta_,norm2)
       call add_to_lanczos_phonon(norm2,e_state,alfa_,beta_,istate)
       deallocate(alfa_,beta_)
       if(allocated(vvinit))deallocate(vvinit)
       if(allocated(v_state))deallocate(v_state)
    end do
    return
  end subroutine lanc_build_gf_phonon_main


  subroutine add_to_lanczos_phonon(vnorm2,Ei,alanc,blanc,istate)
    real(8)                                    :: vnorm2,Ei,Ej,Egs,pesoF,pesoAB,pesoBZ,de,peso
    integer                                    :: nlanc
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,istate
    complex(8)                                 :: iw
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
    call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
    !
    ! call allocate_GFmatrix(impDmatrx,istate,1,Nexc=Nlanc)
    !
    do j=1,nlanc
       Ej     = diag(j)
       dE     = Ej-Ei
       pesoAB = Z(1,j)*Z(1,j)
       peso   = pesoF*pesoAB*pesoBZ
       !COPYPASTE FROM NORMAL, TO CHECK
       ! impDmatrix%state(istate)%channel(ichan)%weight(j) = peso
       ! impDmatrix%state(istate)%channel(ichan)%poles(j)  = de
       !
       ! the correct behavior for beta*dE << 1 is recovered only by assuming that v_n is still finite
       ! beta*dE << v_n for v_n--> 0 slower. First limit beta*dE--> 0 and only then v_n -->0.
       ! This ensures that the correct null contribution is obtained.
       ! So we impose that: if (beta*dE is larger than a small qty) we sum up the contribution, else
       ! we do not include the contribution (because we are in the situation described above).
       ! For the real-axis case this problem is circumvented by the usual i*0+ = xi*eps
       if(beta*dE > 1d-3)impDmats_ph(0)=impDmats_ph(0) - peso*2*(1d0-exp(-beta*dE))/dE
       do i=1,Lmats
          impDmats_ph(i)=impDmats_ph(i) - peso*(1d0-exp(-beta*dE))*2d0*dE/(vm(i)**2+dE**2)
       enddo
       do i=1,Lreal
          impDreal_ph(i)=impDreal_ph(i) + peso*(1d0-exp(-beta*dE))*(1d0/(dcmplx(vr(i),eps) - dE) - 1d0/(dcmplx(vr(i),eps) + dE))
       enddo
    enddo
  end subroutine add_to_lanczos_phonon






  !################################################################
  !################################################################
  !################################################################
  !################################################################









END MODULE ED_GF_SUPERC
