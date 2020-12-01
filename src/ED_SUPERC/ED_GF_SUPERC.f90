MODULE ED_GF_SUPERC
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
  USE ED_HAMILTONIAN_SUPERC
  implicit none
  private

  public :: build_gf_superc
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
  complex(8),dimension(:),pointer       :: state_cvec
  real(8)                               :: state_e

  !AUX GF
  !=========================================================
  complex(8),allocatable,dimension(:,:) :: auxGmats,auxGreal

contains



  !+------------------------------------------------------------------+
  !                        SUPERC
  !+------------------------------------------------------------------+
  subroutine build_gf_superc()
    integer    :: iorb,jorb,ispin,i,isign
    complex(8) :: barGmats(Norb,Lmats),barGreal(Norb,Lreal)
    !
    !
    if(.not.allocated(auxGmats))allocate(auxGmats(3,Lmats))
    if(.not.allocated(auxGreal))allocate(auxGreal(3,Lreal))
    auxgmats=zero
    auxGreal=zero
    barGmats=zero
    barGreal=zero
    !
    ispin=1                       !in this channel Nspin=2 is forbidden. check in ED_AUX_FUNX.
    do iorb=1,Norb
       auxGmats=zero
       auxGreal=zero
       write(LOGfile,"(A)")"Get G&F_l"//str(iorb)//"_s"//str(ispin)
       if(MPIMASTER)call start_timer()
       call lanc_build_gf_superc_c(iorb)
       if(MPIMASTER)call stop_timer(unit=logfile)
       !
       impGmats(ispin,ispin,iorb,iorb,:) = auxGmats(1,:) !this is G_{iorb,iorb} = G_{up,up;iorb,iorb}
       impGreal(ispin,ispin,iorb,iorb,:) = auxGreal(1,:)
       barGmats(                 iorb,:) = auxGmats(2,:) !this is \bar{G}_{iorb,iorb} = \bar{G}_{dw,dw;iorb,iorb}
       barGreal(                 iorb,:) = auxGreal(2,:)
       impFmats(ispin,ispin,iorb,iorb,:) = 0.5d0*(auxGmats(3,:)-auxGmats(1,:)-auxGmats(2,:))
       impFreal(ispin,ispin,iorb,iorb,:) = 0.5d0*(auxGreal(3,:)-auxGreal(1,:)-auxGreal(2,:))
       !
       ! Comment out this and following lines marked with >anomal to use the more general algorithm
       ! for the evaluation of the anomalous gf
       ! >ANOMAL
       ! impFmats(ispin,ispin,iorb,iorb,:) = 0.5d0*(auxGmats(3,:)-(one-xi)*auxGmats(1,:)-(one-xi)*auxGmats(2,:))
       ! impFreal(ispin,ispin,iorb,iorb,:) = 0.5d0*(auxGreal(3,:)-(one-xi)*auxGreal(1,:)-(one-xi)*auxGreal(2,:))
       ! <ANOMAL
    enddo
    !
    !now we add the other mixed/anomalous GF in for the bath_type="hybrid" case
    if(bath_type=='hybrid')then
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             write(LOGfile,"(A)")"Get G_l"//str(iorb)//"_m"//str(jorb)//"_s"//str(ispin)
             if(MPIMASTER)call start_timer()
             call lanc_build_gf_superc_mix_c(iorb,jorb)
             if(MPIMASTER)call stop_timer()
             impGmats(ispin,ispin,iorb,jorb,:) = auxGmats(3,:)
             impGreal(ispin,ispin,iorb,jorb,:) = auxGreal(3,:)
          enddo
       enddo
       !
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             impFmats(ispin,ispin,iorb,jorb,:) = 0.5d0*( impGmats(ispin,ispin,iorb,jorb,:) - &
                  (one-xi)*impGmats(ispin,ispin,iorb,iorb,:) - (one-xi)*barGmats(jorb,:) )
             impFreal(ispin,ispin,iorb,jorb,:) = 0.5d0*( impGreal(ispin,ispin,iorb,jorb,:) - &
                  (one-xi)*impGreal(ispin,ispin,iorb,iorb,:) - (one-xi)*barGreal(jorb,:) )
          enddo
       enddo
    endif
    deallocate(auxGmats,auxGreal)
  end subroutine build_gf_superc










  subroutine lanc_build_gf_superc_c(iorb)
    integer      :: iorb
    type(sector) :: sectorI,sectorJ
    !
    ialfa = 1
    !
    do istate=1,state_list%size
       isector    =  es_return_sector(state_list,istate)
       state_e    =  es_return_energy(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          state_cvec => es_return_cvector(MpiComm,state_list,istate) 
       else
          state_cvec => es_return_cvector(state_list,istate)
       endif
#else
       state_cvec => es_return_cvector(state_list,istate)
#endif
       !
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          if(ed_verbose>=3)write(LOGfile,"(A,I6,I6)")&
               'From sector  :',isector,sectorI%Sz
       endif
       !
       !EVALUATE c^+_{up,iorb}|v> --> Gaux(1) = G_{iorb,iorb}
       jsector = getCDGsector(ialfa,1,isector)
       if(jsector/=0)then 
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(ed_verbose>=3)write(LOGfile,"(A23,I3)")'apply c^+_a,up:',sectorJ%Sz
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             do i=1,sectorI%Dim
                call apply_op_CDG(i,j,sgn,iorb,ialfa,1,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*state_cvec(i)
             enddo
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=zero
          endif
          !
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,state_e,alfa_,beta_,1,ichan=1)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
       endif
       !
       !EVALUATE c_{up,iorb}|v> --> Gaux(1) = G_{iorb,iorb}
       jsector = getCsector(ialfa,1,isector)
       if(jsector/=0)then 
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(ed_verbose>=3)write(LOGfile,"(A23,I3)")'apply c_a,up:',sectorJ%Sz
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             do i=1,sectorI%Dim
                call apply_op_C(i,j,sgn,iorb,ialfa,1,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*state_cvec(i)
             enddo
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=zero
          endif
          !
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,state_e,alfa_,beta_,-1,ichan=1)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
       endif
       !
       !EVALUATE c_{dw,iorb}|v> --> Gaux(2) = barG_{iorb,iorb}
       jsector = getCsector(ialfa,2,isector)
       if(jsector/=0)then 
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(ed_verbose>=3)write(LOGfile,"(A23,I3)")'apply c_a,dw:',sectorJ%Sz
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             do i=1,sectorI%Dim
                call apply_op_C(i,j,sgn,iorb,ialfa,1,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*state_cvec(i)
             enddo
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=zero
          endif
          !
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,state_e,alfa_,beta_,1,ichan=2)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
       endif
       !
       !EVALUATE c^+_{dw,iorb}|v> --> Gaux(2) = barG_{iorb,iorb}
       jsector = getCDGsector(ialfa,2,isector)
       if(jsector/=0)then 
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(ed_verbose>=3)write(LOGfile,"(A23,I3)")'apply c^+_a,dw:',sectorJ%Sz
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             do i=1,sectorI%Dim
                call apply_op_CDG(i,j,sgn,iorb,ialfa,2,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*state_cvec(i)  !build the cdg_up|gs> state
             enddo
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=zero
          endif
          !
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,state_e,alfa_,beta_,-1,ichan=2)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
       endif
       !
       !
       !
       !EVALUATE [c^+_{up,iorb} + c_{dw,iorb}]|gs> --> A_{iorb,iorb}
       isz = getsz(isector)
       if(isz<Ns)then
          jsector = getsector(isz+1,1)
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(ed_verbose>=3)write(LOGfile,"(A23,I3)")'apply c^+_a,up + c_a,dw:',sectorJ%Sz
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             do i=1,sectorI%Dim
                call apply_op_CDG(i,j,sgn,iorb,ialfa,1,sectorI,sectorJ)!c^+_a,up
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*state_cvec(i)
             enddo
             do i=1,sectorI%Dim
                call apply_op_C(i,j,sgn,iorb,ialfa,2,sectorI,sectorJ) !c_a,dw
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = vvinit(j) + sgn*state_cvec(i)
             enddo
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=zero
          endif
          !
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,state_e,alfa_,beta_,1,ichan=3)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
       endif
       !
       !EVALUATE [c_{up,iorb} + c^+_{dw,iorb}]|gs>  --> A_{iorb,iorb}
       isz = getsz(isector)
       if(isz>-Ns)then
          jsector = getsector(isz-1,1)
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(ed_verbose>=3)write(LOGfile,"(A23,I3)")'apply c_a,up + c^+_a,dw:',sectorJ%Sz
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             do i=1,sectorI%Dim
                call apply_op_C(i,j,sgn,iorb,ialfa,1,sectorI,sectorJ) !c_a,up
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*state_cvec(i)
             enddo
             do i=1,sectorI%Dim
                call apply_op_CDG(i,j,sgn,iorb,ialfa,2,sectorI,sectorJ) !c^+_a,dw
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = vvinit(j) + sgn*state_cvec(i)
             enddo
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=zero
          endif
          !
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,state_e,alfa_,beta_,-1,ichan=3)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
       endif
       !
       if(MpiMaster)call delete_sector(sectorI)
#ifdef _MPI
       if(MpiStatus)then
          if(associated(state_cvec))deallocate(state_cvec)
       else
          if(associated(state_cvec))nullify(state_cvec)
       endif
#else
       if(associated(state_cvec))nullify(state_cvec)
#endif
       !
       !
    enddo
    return
  end subroutine lanc_build_gf_superc_c





  subroutine lanc_build_gf_superc_mix_c(iorb,jorb)
    integer                :: iorb,jorb
    type(sector)           :: sectorI,sectorJ
    !
    !
    ialfa = 1
    !
    !
    do istate=1,state_list%size
       isector    =  es_return_sector(state_list,istate)
       state_e    =  es_return_energy(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          state_cvec => es_return_cvector(MpiComm,state_list,istate) 
       else
          state_cvec => es_return_cvector(state_list,istate)
       endif
#else
       state_cvec => es_return_cvector(state_list,istate)
#endif
       !
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          if(ed_verbose>=3)write(LOGfile,"(A,I6,I6)")&
               'From sector  :',isector,sectorI%Sz
       endif
       !
       !EVALUATE [c^+_{up,iorb} + c_{dw,jorb}]|gs> --> A_{iorb,jorb}
       isz = getsz(isector)
       if(isz<Ns)then
          jsector = getsector(isz+1,1)
          !
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(ed_verbose>=3)write(LOGfile,"(A23,I3)")'apply c^+_a,up + c_b,dw:',sectorJ%Sz
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             do i=1,sectorI%Dim
                call apply_op_CDG(i,j,sgn,iorb,ialfa,1,sectorI,sectorJ)!c^+_a,up
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*state_cvec(i)
             enddo
             do i=1,sectorI%Dim
                call apply_op_C(i,j,sgn,jorb,ialfa,2,sectorI,sectorJ) !c_b,dw
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = vvinit(j) + sgn*state_cvec(i)
             enddo
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=zero
          endif
          !
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,state_e,alfa_,beta_,1,ichan=3)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
       endif
       !
       !EVALUATE [c_{up,iorb} + c^+_{dw,jorb}]|gs>  --> A_{iorb,jorb}
       isz = getsz(isector)
       if(isz>-Ns)then
          jsector = getsector(isz-1,1)
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(ed_verbose>=3)write(LOGfile,"(A23,I3)")'apply c_a,up + c^+_b,dw:',sectorJ%Sz
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             do i=1,sectorI%Dim
                call apply_op_C(i,j,sgn,iorb,ialfa,1,sectorI,sectorJ) !c_a,up
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*state_cvec(i)
             enddo
             do i=1,sectorI%Dim
                call apply_op_CDG(i,j,sgn,jorb,ialfa,2,sectorI,sectorJ) !c^+_b,dw
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = vvinit(j) + sgn*state_cvec(i)
             enddo
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=zero
          endif
          !
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,state_e,alfa_,beta_,-1,ichan=3)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
       endif
       !
       !EVALUATE [c^+_{up,iorb} + xi*c_{dw,jorb}]|gs> --> -xi*B_{iorb,jorb}
       isz = getsz(isector)
       if(isz<Ns)then
          jsector = getsector(isz+1,1)
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(ed_verbose>=3)write(LOGfile,"(A23,I3)")'apply c^+_a,up + xi*c_b,dw:',sectorJ%Sz
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             do i=1,sectorI%Dim
                call apply_op_CDG(i,j,sgn,iorb,ialfa,1,sectorI,sectorJ)!c^+_a,up
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*state_cvec(i)
             enddo
             do i=1,sectorI%Dim
                call apply_op_C(i,j,sgn,jorb,ialfa,2,sectorI,sectorJ) !xi*c_b,dw
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = vvinit(j) + xi*sgn*state_cvec(i)
             enddo
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=zero
          endif
          !
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(-xi*norm2,state_e,alfa_,beta_,1,ichan=3)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
       endif
       !
       !EVALUATE [c_{up,iorb} - xi*c^+_{dw,jorb}]|gs> --> -xi*B_{iorb,jorb}
       isz = getsz(isector)
       if(isz>-Ns)then
          jsector = getsector(isz-1,1)
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(ed_verbose>=3)write(LOGfile,"(A23,I3)")'apply c_a,up - xi*c^+_b,dw:',sectorJ%Sz
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             do i=1,sectorI%Dim
                call apply_op_C(i,j,sgn,iorb,ialfa,1,sectorI,sectorJ)!c_a,up
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*state_cvec(i)
             enddo
             do i=1,sectorI%Dim
                call apply_op_CDG(i,j,sgn,jorb,ialfa,2,sectorI,sectorJ) !-xi*c^+_b,dw
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = vvinit(j) - xi*sgn*state_cvec(i)
             enddo
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=zero
          endif
          !
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(-xi*norm2,state_e,alfa_,beta_,-1,ichan=3)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
       endif
       !
       if(MpiMaster)call delete_sector(sectorI)
#ifdef _MPI
       if(MpiStatus)then
          if(associated(state_cvec))deallocate(state_cvec)
       else
          if(associated(state_cvec))nullify(state_cvec)
       endif
#else
       if(associated(state_cvec))nullify(state_cvec)
#endif

       !
    enddo
    !
    return
  end subroutine lanc_build_gf_superc_mix_c






  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################








  subroutine add_to_lanczos_gf_superc(vnorm2,Ei,alanc,blanc,isign,ichan)
    complex(8)                                 :: vnorm2,pesoBZ,peso
    real(8)                                    :: Ei,Egs,de
    integer                                    :: nlanc,itype
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    integer                                    :: isign,ichan
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
    complex(8)                                 :: iw
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
    !
    ! pesoBZ = vnorm2/zeta_function
    ! if(finiteT)pesoBZ = vnorm2*exp(-beta*(Ei-Egs))/zeta_function
    !
    !Only the nodes in Mpi_Comm_Group did get the alanc,blanc.
    !However after delete_sectorHv MpiComm returns to be the global one
    !so we can safely Bcast the alanc,blanc (known only to the operative group)
    !to every nodes. The master is in charge of this (as a
    !participant of the operative group)
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,alanc)
       call Bcast_MPI(MpiComm,blanc)
    endif
#endif
    !
    ! itype=(3+isign)/2
    ! diag             = 0.d0
    ! subdiag          = 0.d0
    ! Z                = eye(Nlanc)
    ! diag(1:Nlanc)    = alanc(1:Nlanc)
    ! subdiag(2:Nlanc) = blanc(2:Nlanc)
    ! call tql2(Nlanc,diag,subdiag,Z,ierr)
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
    !
    do j=1,nlanc
       de = diag(j)-Ei
       peso = pesoBZ*Z(1,j)*Z(1,j)
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





  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################







  !+------------------------------------------------------------------+
  !PURPOSE  : Build the Self-energy functions, SUPERC case
  !+------------------------------------------------------------------+
  subroutine build_sigma_superc
    integer                                               :: i,ispin,iorb
    real(8)                                               :: det_mats(Lmats)
    complex(8)                                            :: det_real(Lreal)
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats)     :: invG0mats,invF0mats,invGmats,invFmats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal)     :: invG0real,invF0real,invGreal,invFreal
    complex(8),dimension(2*Nspin*Norb,2*Nspin*Norb)       :: invGimp
    !
    !
    invG0mats = zero
    invF0mats = zero
    invGmats  = zero
    invFmats  = zero
    invG0real = zero
    invF0real = zero
    invGreal  = zero
    invFreal  = zero
    !
    !Get G0^-1,F0^-1
    ispin=1
    invG0mats(:,:,:,:,:) = invg0_bath_function(dcmplx(0d0,wm(:)),dmft_bath)
    invF0mats(:,:,:,:,:) = invf0_bath_function(dcmplx(0d0,wm(:)),dmft_bath)
    !
    invG0real(:,:,:,:,:) = invg0_bath_function(dcmplx(wr(:),eps),dmft_bath,axis='real')
    invF0real(:,:,:,:,:) = invf0_bath_function(dcmplx(wr(:),eps),dmft_bath,axis='real')
    !
    select case(bath_type)
    case default
       !      
       !Get Gimp^-1
       do iorb=1,Norb
          det_mats  =  abs(impGmats(ispin,ispin,iorb,iorb,:))**2 + (impFmats(ispin,ispin,iorb,iorb,:))**2
          invGmats(ispin,ispin,iorb,iorb,:) = conjg(impGmats(ispin,ispin,iorb,iorb,:))/det_mats
          invFmats(ispin,ispin,iorb,iorb,:) = impFmats(ispin,ispin,iorb,iorb,:)/det_mats
          !
          det_real  = -impGreal(ispin,ispin,iorb,iorb,:)*conjg(impGreal(ispin,ispin,iorb,iorb,Lreal:1:-1)) - impFreal(ispin,ispin,iorb,iorb,:)**2
          invGreal(ispin,ispin,iorb,iorb,:) =  -conjg(impGreal(ispin,ispin,iorb,iorb,Lreal:1:-1))/det_real(:)
          invFreal(ispin,ispin,iorb,iorb,:) =  -impFreal(ispin,ispin,iorb,iorb,:)/det_real(:)
       enddo
       !Get Sigma functions: Sigma= G0^-1 - G^-1
       impSmats=zero
       impSAmats=zero
       impSreal=zero
       impSAreal=zero
       do iorb=1,Norb
          impSmats(ispin,ispin,iorb,iorb,:)  = invG0mats(ispin,ispin,iorb,iorb,:) - invGmats(ispin,ispin,iorb,iorb,:)
          impSAmats(ispin,ispin,iorb,iorb,:) = invF0mats(ispin,ispin,iorb,iorb,:) - invFmats(ispin,ispin,iorb,iorb,:)
          !
          impSreal(ispin,ispin,iorb,iorb,:)  = invG0real(ispin,ispin,iorb,iorb,:) - invGreal(ispin,ispin,iorb,iorb,:)
          impSAreal(ispin,ispin,iorb,iorb,:) = invF0real(ispin,ispin,iorb,iorb,:) - invFreal(ispin,ispin,iorb,iorb,:)
       enddo
       !
    case ("hybrid")
       !
       !Get Gimp^-1
       do i=1,Lmats
          invGimp=zero
          invGimp(1:Norb,1:Norb)               = impGmats(ispin,ispin,:,:,i)
          invGimp(1:Norb,Norb+1:2*Norb)        = impFmats(ispin,ispin,:,:,i)
          invGimp(Norb+1:2*Norb,1:Norb)        = impFmats(ispin,ispin,:,:,i)
          invGimp(Norb+1:2*Norb,Norb+1:2*Norb) =-conjg(impGmats(ispin,ispin,:,:,i))
          call inv(invGimp)
          invGmats(ispin,ispin,:,:,i) = invGimp(1:Norb,1:Norb)
          invFmats(ispin,ispin,:,:,i) = invGimp(1:Norb,Norb+1:2*Norb)
       enddo
       do i=1,Lreal
          invGimp=zero
          invGimp(1:Norb,1:Norb)               = impGreal(ispin,ispin,:,:,i)
          invGimp(1:Norb,Norb+1:2*Norb)        = impFreal(ispin,ispin,:,:,i)
          invGimp(Norb+1:2*Norb,1:Norb)        = impFreal(ispin,ispin,:,:,i)
          invGimp(Norb+1:2*Norb,Norb+1:2*Norb) =-conjg(impGreal(ispin,ispin,:,:,Lreal-i+1))
          call inv(invGimp)
          invGreal(ispin,ispin,:,:,i) =  invGimp(1:Norb,1:Norb)
          invFreal(ispin,ispin,:,:,i) =  invGimp(1:Norb,Norb+1:2*Norb)
       enddo
       !Get Sigma functions: Sigma= G0^-1 - G^-1
       impSmats=zero
       impSAmats=zero
       impSreal=zero
       impSAreal=zero
       !
       impSmats(ispin,ispin,:,:,:)  = invG0mats(ispin,ispin,:,:,:) - invGmats(ispin,ispin,:,:,:)
       impSAmats(ispin,ispin,:,:,:) = invF0mats(ispin,ispin,:,:,:) - invFmats(ispin,ispin,:,:,:)
       !
       impSreal(ispin,ispin,:,:,:)  = invG0real(ispin,ispin,:,:,:) - invGreal(ispin,ispin,:,:,:)
       impSAreal(ispin,ispin,:,:,:) = invF0real(ispin,ispin,:,:,:) - invFreal(ispin,ispin,:,:,:)
       !
    end select
    !
    !Get G0and:
    impG0mats(:,:,:,:,:) = g0and_bath_function(dcmplx(0d0,wm(:)),dmft_bath)
    impF0mats(:,:,:,:,:) = f0and_bath_function(dcmplx(0d0,wm(:)),dmft_bath)
    !
    impG0real(:,:,:,:,:) = g0and_bath_function(dcmplx(wr(:),eps),dmft_bath,axis='real')
    impF0real(:,:,:,:,:) = f0and_bath_function(dcmplx(wr(:),eps),dmft_bath,axis='real')
    !!
    !
    !
  end subroutine build_sigma_superc





END MODULE ED_GF_SUPERC
