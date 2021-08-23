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
    do iorb=1,Norb
       auxGmats=zero
       auxGreal=zero
       write(LOGfile,"(A)")"Get G&F_l"//str(iorb)//"_s"//str(ispin)
       if(MPIMASTER)call start_timer()
       call allocate_GFmatrix(impGmatrix(ispin,ispin,iorb,iorb),Nstate=state_list%size)
       call lanc_build_gf_superc_c(iorb)
       if(MPIMASTER)call stop_timer(unit=logfile)
#ifdef _DEBUG
       write(Logfile,"(A)")""
#endif
       !
       impGmats(ispin,ispin,iorb,iorb,:) = auxGmats(1,:) !this is G_{up,up;iorb,iorb}
       impGreal(ispin,ispin,iorb,iorb,:) = auxGreal(1,:)  
       barGmats(                 iorb,:) = auxGmats(2,:) !this is G_{dw,dw;iorb,iorb}
       barGreal(                 iorb,:) = auxGreal(2,:)
       !                                   auxGmats(3,:) !this is G_oo -xi*G_pp = A_aa -xi*B_aa
       !                                   auxGreal(3,:)
       impFmats(ispin,ispin,iorb,iorb,:) = 0.5d0*(auxGmats(3,:)-auxGmats(1,:)-auxGmats(2,:))
       impFreal(ispin,ispin,iorb,iorb,:) = 0.5d0*(auxGreal(3,:)-auxGreal(1,:)-auxGreal(2,:))
       !
       ! Comment out lines marked with >/<ANOMAL to use the more general algorithm to evaluate anomalous gf
       ! >ANOMAL
       ! impFmats(ispin,ispin,iorb,iorb,:) = 0.5d0*(auxGmats(3,:)-(one-xi)*auxGmats(1,:)-(one-xi)*auxGmats(2,:))
       ! impFreal(ispin,ispin,iorb,iorb,:) = 0.5d0*(auxGreal(3,:)-(one-xi)*auxGreal(1,:)-(one-xi)*auxGreal(2,:))
       ! <ANOMAL
    enddo
    !
    !now we add the other mixed/anomalous GF in for the bath_type="hybrid" case
    if(bath_type=='hybrid')then
       do iorb=1,Norb
          do jorb=1,Norb
             if(iorb==jorb)cycle
             write(LOGfile,"(A)")"Get G_l"//str(iorb)//"_m"//str(jorb)//"_s"//str(ispin)
             if(MPIMASTER)call start_timer()
             call allocate_GFmatrix(impGmatrix(ispin,ispin,iorb,jorb),Nstate=state_list%size)
             call lanc_build_gf_superc_mix_c(iorb,jorb)
             if(MPIMASTER)call stop_timer()
#ifdef _DEBUG
             write(Logfile,"(A)")""
#endif
             !
             impGmats(ispin,ispin,iorb,jorb,:) = auxGmats(4,:) !G_oo -xi*G_pp=A_ab-xi*B_ab
             impGreal(ispin,ispin,iorb,jorb,:) = auxGreal(4,:)
          enddo
       enddo
       !
       do iorb=1,Norb
          do jorb=1,Norb
             if(iorb==jorb)cycle
             impFmats(ispin,ispin,iorb,jorb,:) = 0.5d0*( impGmats(ispin,ispin,iorb,jorb,:) - &
                  (one-xi)*impGmats(ispin,ispin,iorb,iorb,:) - (one-xi)*barGmats(jorb,:) )
             impFreal(ispin,ispin,iorb,jorb,:) = 0.5d0*( impGreal(ispin,ispin,iorb,jorb,:) - &
                  (one-xi)*impGreal(ispin,ispin,iorb,iorb,:) - (one-xi)*barGreal(jorb,:) )
          enddo
       enddo
    endif
    deallocate(auxGmats,auxGreal)
  end subroutine build_gf_superc


  subroutine rebuild_gf_superc()
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
       write(LOGfile,"(A)")"Get G&F_l"//str(iorb)//"_s"//str(ispin)
       call rebuild_gf_superc_diag(iorb)
       !
       impGmats(ispin,ispin,iorb,iorb,:) = auxGmats(1,:) !this is G_{up,up;iorb,iorb}
       impGreal(ispin,ispin,iorb,iorb,:) = auxGreal(1,:)  
       barGmats(                 iorb,:) = auxGmats(2,:) !this is G_{dw,dw;iorb,iorb}
       barGreal(                 iorb,:) = auxGreal(2,:)
       !                                   auxGmats(3,:) !this is G_oo -xi*G_pp = A_aa -xi*B_aa
       !                                   auxGreal(3,:)
       impFmats(ispin,ispin,iorb,iorb,:) = 0.5d0*(auxGmats(3,:)-auxGmats(1,:)-auxGmats(2,:))
       impFreal(ispin,ispin,iorb,iorb,:) = 0.5d0*(auxGreal(3,:)-auxGreal(1,:)-auxGreal(2,:))
       !
       ! Comment out lines marked with >/<ANOMAL to use the more general algorithm to evaluate anomalous gf
       ! >ANOMAL
       ! impFmats(ispin,ispin,iorb,iorb,:) = 0.5d0*(auxGmats(3,:)-(one-xi)*auxGmats(1,:)-(one-xi)*auxGmats(2,:))
       ! impFreal(ispin,ispin,iorb,iorb,:) = 0.5d0*(auxGreal(3,:)-(one-xi)*auxGreal(1,:)-(one-xi)*auxGreal(2,:))
       ! <ANOMAL
    enddo
    !
    !now we add the other mixed/anomalous GF in for the bath_type="hybrid" case
    if(bath_type=='hybrid')then
       do iorb=1,Norb
          do jorb=1,Norb
             if(iorb==jorb)cycle
             write(LOGfile,"(A)")"Get G_l"//str(iorb)//"_m"//str(jorb)//"_s"//str(ispin)
             call rebuild_gf_superc_mix(iorb,jorb)
             impGmats(ispin,ispin,iorb,jorb,:) = auxGmats(4,:) !This is G_oo -xi*G_pp = A_ab -xi*B_ab
             impGreal(ispin,ispin,iorb,jorb,:) = auxGreal(4,:)
          enddo
       enddo
       !
       do iorb=1,Norb
          do jorb=1,Norb
             if(iorb==jorb)cycle
             impFmats(ispin,ispin,iorb,jorb,:) = 0.5d0*( impGmats(ispin,ispin,iorb,jorb,:) - &
                  (one-xi)*impGmats(ispin,ispin,iorb,iorb,:) - (one-xi)*barGmats(jorb,:) )
             impFreal(ispin,ispin,iorb,jorb,:) = 0.5d0*( impGreal(ispin,ispin,iorb,jorb,:) - &
                  (one-xi)*impGreal(ispin,ispin,iorb,iorb,:) - (one-xi)*barGreal(jorb,:) )
          enddo
       enddo
    endif
    deallocate(auxGmats,auxGreal)
  end subroutine rebuild_gf_superc






  !################################################################
  !################################################################
  !################################################################
  !################################################################




  subroutine lanc_build_gf_superc_c(iorb)
    integer      :: iorb
    type(sector) :: sectorI,sectorJ
    !
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")&
         "DEBUG lanc_build_gf SUPERC DIAG: build diagonal GF l"//str(iorb)
#endif
    !
    ialfa = 1
    !
    do istate=1,state_list%size
       !
       call allocate_GFmatrix(impGmatrix(1,1,iorb,iorb),istate,Nchan=6) !2*[uu,dd,aux]
       !<ANOMAL
       ! call allocate_GFmatrix(impGmatrix(1,1,iorb,iorb),istate,Nchan=8) !2*[uu,dd,aux1,aux2]
       !>ANOMAL
       !
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
       !EVALUATE c^+_{up,iorb}|v> --> Gaux(1) = G^>_{up,up;iorb,iorb}
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
          call add_to_lanczos_gf_superc(one*norm2,state_e,alfa_,beta_,1,1,iorb,iorb,1,istate)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,iorb,iorb),istate,1,Nexc=0)
       endif
       !
       !EVALUATE c_{up,iorb}|v> --> Gaux(1) = G^<_{up,up;iorb,iorb}
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
          call add_to_lanczos_gf_superc(one*norm2,state_e,alfa_,beta_,-1,1,iorb,iorb,2,istate)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,iorb,iorb),istate,2,Nexc=0)
       endif
       !
       !
       !
       !EVALUATE c_{dw,iorb}|v> --> Gaux(2) = G^>_{dw,dw;iorb,iorb}
       jsector = getCsector(ialfa,2,isector)
       if(jsector/=0)then 
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(ed_verbose>=3)write(LOGfile,"(A23,I3)")'apply c_a,dw:',sectorJ%Sz
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             do i=1,sectorI%Dim
                call apply_op_C(i,j,sgn,iorb,ialfa,2,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*state_cvec(i)
             enddo
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=zero
          endif
          !
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,state_e,alfa_,beta_,1,2,iorb,iorb,3,istate)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,iorb,iorb),istate,3,Nexc=0)
       endif
       !
       !EVALUATE c^+_{dw,iorb}|v> --> Gaux(2) = G^<_{dw,dw;iorb,iorb}
       jsector = getCDGsector(ialfa,2,isector)
       if(jsector/=0)then 
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(ed_verbose>=3)write(LOGfile,"(A23,I3)")'apply c^+_a,dw:',sectorJ%Sz
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             do i=1,sectorI%Dim
                call apply_op_CDG(i,j,sgn,iorb,ialfa,2,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*state_cvec(i)
             enddo
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=zero
          endif
          !
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,state_e,alfa_,beta_,-1,2,iorb,iorb,4,istate)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,iorb,iorb),istate,4,Nexc=0)
       endif
       !
       !
       !
       !EVALUATE [c^+_{up,iorb} + c_{dw,iorb}]|gs> --> A^>_{iorb,iorb}
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
          call add_to_lanczos_gf_superc(one*norm2,state_e,alfa_,beta_,1,3,iorb,iorb,5,istate)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,iorb,iorb),istate,5,Nexc=0)
       endif
       !
       !EVALUATE [c_{up,iorb} + c^+_{dw,iorb}]|gs>  --> A^<_{iorb,iorb}
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
          call add_to_lanczos_gf_superc(one*norm2,state_e,alfa_,beta_,-1,3,iorb,iorb,6,istate)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,iorb,iorb),istate,6,Nexc=0)
       endif
       !
       !
       !
       ! !<ANOMAL
       ! !EVALUATE [c^+_{up,iorb} + xi*c_{dw,iorb}]|gs> --> -xi*B^>_{iorb,iorb}
       ! isz = getsz(isector)
       ! if(isz<Ns)then
       !    jsector = getsector(isz+1,1)
       !    if(MpiMaster)then
       !       call build_sector(jsector,sectorJ)
       !       if(ed_verbose>=3)write(LOGfile,"(A23,I3)")'apply c^+_a,up + xi*c_a,dw:',sectorJ%Sz
       !       allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
       !       do i=1,sectorI%Dim
       !          call apply_op_CDG(i,j,sgn,iorb,ialfa,1,sectorI,sectorJ)!c^+_a,up
       !          if(sgn==0d0.OR.j==0)cycle
       !          vvinit(j) = sgn*state_cvec(i)
       !       enddo
       !       do i=1,sectorI%Dim
       !          call apply_op_C(i,j,sgn,iorb,ialfa,2,sectorI,sectorJ) !xi*c_a,dw
       !          if(sgn==0d0.OR.j==0)cycle
       !          vvinit(j) = vvinit(j) + xi*sgn*state_cvec(i)
       !       enddo
       !       call delete_sector(sectorJ)
       !    else
       !       allocate(vvinit(1));vvinit=zero
       !    endif
       !    !
       !    call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
       !    call add_to_lanczos_gf_superc(-xi*norm2,state_e,alfa_,beta_,1,3,iorb,iorb,7,istate)
       !    deallocate(alfa_,beta_)
       !    if(allocated(vvinit))deallocate(vvinit)
       ! else
       !  call allocate_GFmatrix(impGmatrix(1,1,iorb,iorb),istate,7,Nexc=0)
       ! endif
       ! !
       ! !EVALUATE [c_{up,iorb} - xi*c^+_{dw,iorb}]|gs> --> -xi*B^<_{iorb,iorb}
       ! isz = getsz(isector)
       ! if(isz>-Ns)then
       !    jsector = getsector(isz-1,1)
       !    if(MpiMaster)then
       !       call build_sector(jsector,sectorJ)
       !       if(ed_verbose>=3)write(LOGfile,"(A23,I3)")'apply c_a,up - xi*c^+_a,dw:',sectorJ%Sz
       !       allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
       !       do i=1,sectorI%Dim
       !          call apply_op_C(i,j,sgn,iorb,ialfa,1,sectorI,sectorJ)!c_a,up
       !          if(sgn==0d0.OR.j==0)cycle
       !          vvinit(j) = sgn*state_cvec(i)
       !       enddo
       !       do i=1,sectorI%Dim
       !          call apply_op_CDG(i,j,sgn,iorb,ialfa,2,sectorI,sectorJ) !-xi*c^+_a,dw
       !          if(sgn==0d0.OR.j==0)cycle
       !          vvinit(j) = vvinit(j) - xi*sgn*state_cvec(i)
       !       enddo
       !       call delete_sector(sectorJ)
       !    else
       !       allocate(vvinit(1));vvinit=zero
       !    endif
       !    !
       !    call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
       !    call add_to_lanczos_gf_superc(-xi*norm2,state_e,alfa_,beta_,-1,3,iorb,iorb,8,istate)
       !    deallocate(alfa_,beta_)
       !    if(allocated(vvinit))deallocate(vvinit)
       ! else
       !  call allocate_GFmatrix(impGmatrix(1,1,iorb,iorb),istate,8,Nexc=0)
       ! endif
       ! !>ANOMAL
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
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")&
         "DEBUG lanc_build_gf SUPERC DIAG: build mixed GF l"//str(iorb)//",m"//str(jorb)
#endif
    !
    ialfa = 1
    !
    !
    do istate=1,state_list%size
       !
       call allocate_GFmatrix(impGmatrix(1,1,iorb,jorb),istate,Nchan=4) !2*[aux1,aux2]
       !
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
       !EVALUATE [c^+_{up,iorb} + c_{dw,jorb}]|gs> --> A^>_{iorb,jorb}
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
          call add_to_lanczos_gf_superc(one*norm2,state_e,alfa_,beta_,1,4,iorb,jorb,1,istate)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,iorb,jorb),istate,1,Nexc=0)
       endif
       !
       !EVALUATE [c_{up,iorb} + c^+_{dw,jorb}]|gs>  --> A^<_{iorb,jorb}
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
          call add_to_lanczos_gf_superc(one*norm2,state_e,alfa_,beta_,-1,4,iorb,jorb,2,istate)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,iorb,jorb),istate,2,Nexc=0)
       endif
       !
       !
       !
       !EVALUATE [c^+_{up,iorb} + xi*c_{dw,jorb}]|gs> --> -xi*B^>_{iorb,jorb}
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
          call add_to_lanczos_gf_superc(-xi*norm2,state_e,alfa_,beta_,1,4,iorb,jorb,3,istate)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,iorb,jorb),istate,3,Nexc=0)
       endif
       !
       !EVALUATE [c_{up,iorb} - xi*c^+_{dw,jorb}]|gs> --> -xi*B^<_{iorb,jorb}
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
          call add_to_lanczos_gf_superc(-xi*norm2,state_e,alfa_,beta_,-1,4,iorb,jorb,4,istate)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,iorb,jorb),istate,4,Nexc=0)
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







  subroutine add_to_lanczos_gf_superc(vnorm2,Ei,alanc,blanc,isign,ichan,iorb,jorb,ic,istate)
    complex(8)                                 :: vnorm2,pesoBZ,peso
    real(8)                                    :: Ei,Egs,de
    integer                                    :: nlanc,itype
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    integer                                    :: isign,ichan,iorb,jorb,istate,ic
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
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
    call allocate_GFmatrix(impGmatrix(1,1,iorb,jorb),istate,ic,Nlanc)
    !
    do j=1,nlanc
       de = diag(j)-Ei
       peso = pesoBZ*Z(1,j)*Z(1,j)
       !
       impGmatrix(1,1,iorb,jorb)%state(istate)%channel(ic)%weight(j) = peso
       impGmatrix(1,1,iorb,jorb)%state(istate)%channel(ic)%poles(j)  = isign*de
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




  subroutine rebuild_gf_superc_diag(iorb)
    integer,intent(in) :: iorb
    integer            :: Nstates,istate
    integer            :: Nchannels,ic,ichan
    integer            :: Nexcs,iexc
    real(8)            :: peso,de
    !
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")&
         "DEBUG rebuild_gf SUPERC: reconstruct diagonal impurity GFs"
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
          case(5:8);ichan=3
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
  end subroutine rebuild_gf_superc_diag

  subroutine rebuild_gf_superc_mix(iorb,jorb)
    integer,intent(in) :: iorb,jorb
    integer            :: Nstates,istate
    integer            :: Nchannels,ic,ichan
    integer            :: Nexcs,iexc
    real(8)            :: peso,de
    !
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")&
         "DEBUG rebuild_gf SUPERC: reconstruct mixed impurity GFs"
#endif
    !
    !
    if(.not.allocated(impGmatrix(1,1,iorb,jorb)%state)) then
       print*, "ED_GF_SUPERC WARNING: impGmatrix%state not allocated. Nothing to do"
       return
    endif
    !
    ichan = 4
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
  end subroutine rebuild_gf_superc_mix






  !################################################################
  !################################################################
  !################################################################
  !################################################################



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
    case ("hybrid")
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





END MODULE ED_GF_SUPERC
