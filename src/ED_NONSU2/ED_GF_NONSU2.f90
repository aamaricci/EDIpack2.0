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
  USE ED_BATH_FUNCTIONS
  USE ED_SETUP
  USE ED_SECTOR
  USE ED_HAMILTONIAN_NONSU2
  implicit none
  private


  public :: build_gf_nonsu2
  public :: build_sigma_nonsu2


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
  complex(8),dimension(:),pointer       :: state_cvec
  real(8)                               :: state_e

  !AUX GF
  !=========================================================
  complex(8),allocatable,dimension(:,:) :: auxGmats,auxGreal


contains



  !+------------------------------------------------------------------+
  !                            NONSU2
  !+------------------------------------------------------------------+
  subroutine build_gf_nonsu2()
    integer                                     :: iorb,jorb,ispin,jspin,i,io,jo
    logical                                     :: MaskBool
    logical(8),dimension(Nspin,Nspin,Norb,Norb) :: Hmask
    !
    if(.not.allocated(impGmats))stop "build_gf_nonsu2: Gmats not allocated"
    if(.not.allocated(impGreal))stop "build_gf_nonsu2: Greal not allocated"
    impGmats=zero
    impGreal=zero
    !
    !Here we evaluate the same orbital, same spin GF: G_{aa}^{ss}(z)
    do ispin=1,Nspin
       do iorb=1,Norb
          write(LOGfile,"(A)")""
          write(LOGfile,"(A)")"Get G_l"//str(iorb)//str(iorb)//"_s"//str(ispin)//str(ispin)
          if(MPIMASTER)call start_timer()
          call lanc_build_gf_nonsu2_diagOrb_diagSpin(iorb,ispin)
          if(MPIMASTER)call stop_timer(unit=logfile)
       enddo
    enddo
    !
    Hmask=mask_hloc(impHloc,wdiag=.true.,uplo=.true.)
    !
    !same orbital, different spin GF: G_{aa}^{ss'}(z)
    do ispin=1,Nspin
       do jspin=1,Nspin
          if(ispin==jspin)cycle
          do iorb=1,Norb
             MaskBool=.true.   
             if(bath_type=="replica")MaskBool=Hmask(ispin,jspin,iorb,iorb)
             if(.not.MaskBool)cycle
             !
             write(LOGfile,"(A)")""
             write(LOGfile,"(A)")"Get G_l"//str(iorb)//str(iorb)//"_s"//str(ispin)//str(jspin)
             if(MPIMASTER)call start_timer()
             call lanc_build_gf_nonsu2_mixOrb_mixSpin(iorb,iorb,ispin,jspin)
             if(MPIMASTER)call stop_timer(unit=logfile)
          enddo
       enddo
    enddo
    !Here we put the symmetry manipulation
    do ispin=1,Nspin
       do jspin=1,Nspin
          if(ispin==jspin)cycle
          do iorb=1,Norb
             !
             impGmats(ispin,jspin,iorb,iorb,:) = 0.5d0*(impGmats(ispin,jspin,iorb,iorb,:) &
                  - (one+xi)*impGmats(ispin,ispin,iorb,iorb,:) &
                  - (one+xi)*impGmats(jspin,jspin,iorb,iorb,:))
             !
             impGreal(ispin,jspin,iorb,iorb,:) = 0.5d0*(impGreal(ispin,jspin,iorb,iorb,:) &
                  - (one+xi)*impGreal(ispin,ispin,iorb,iorb,:) &
                  - (one+xi)*impGreal(jspin,jspin,iorb,iorb,:))
             !
          enddo
       enddo
    enddo
    !
    !
    select case(bath_type)
    case default;
    case("hybrid","replica")
       !Here we evaluate the different orbital, same spin GF: G_{ab}^{ss}(z)
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                if(iorb==jorb)cycle
                MaskBool=.true.   
                if(bath_type=="replica")MaskBool=Hmask(ispin,ispin,iorb,jorb)
                if(.not.MaskBool)cycle
                !
                write(LOGfile,"(A)")""
                write(LOGfile,"(A)")"Get G_l"//str(iorb)//str(jorb)//"_s"//str(ispin)//str(ispin)
                if(MPIMASTER)call start_timer()
                call lanc_build_gf_nonsu2_mixOrb_mixSpin(iorb,jorb,ispin,ispin)
                if(MPIMASTER)call stop_timer(unit=logfile)
             enddo
          enddo
       enddo
       !Here we put the symmetry manipulation
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                if(iorb==jorb)cycle
                impGmats(ispin,ispin,iorb,jorb,:) = 0.5d0*(impGmats(ispin,ispin,iorb,jorb,:) &
                     - (one+xi)*impGmats(ispin,ispin,iorb,iorb,:) &
                     - (one+xi)*impGmats(ispin,ispin,jorb,jorb,:))
                !
                impGreal(ispin,ispin,iorb,jorb,:) = 0.5d0*(impGreal(ispin,ispin,iorb,jorb,:) &
                     - (one+xi)*impGreal(ispin,ispin,iorb,iorb,:) &
                     - (one+xi)*impGreal(ispin,ispin,jorb,jorb,:))
                !
             enddo
          enddo
       enddo
       !
       !Here we evaluate the different orbital, different spin GF: G_{ab}^{ss'}(z)
       do ispin=1,Nspin
          do jspin=1,Nspin
             if(ispin==jspin)cycle
             do iorb=1,Norb
                do jorb=1,Norb
                   if(iorb==jorb)cycle
                   MaskBool=.true.   
                   if(bath_type=="replica")MaskBool=Hmask(ispin,jspin,iorb,jorb)
                   if(.not.MaskBool)cycle
                   !
                   write(LOGfile,"(A)")""
                   write(LOGfile,"(A)")"Get G_l"//str(iorb)//str(jorb)//"_s"//str(ispin)//str(jspin)
                   if(MPIMASTER)call start_timer()
                   call lanc_build_gf_nonsu2_mixOrb_mixSpin(iorb,jorb,ispin,jspin)
                   if(MPIMASTER)call stop_timer(unit=logfile)
                enddo
             enddo
          enddo
       enddo
       !Here we put the symmetry manipulation
       do ispin=1,Nspin
          do jspin=1,Nspin
             if(ispin==jspin)cycle
             do iorb=1,Norb
                do jorb=1,Norb
                   if(iorb==jorb)cycle
                   impGmats(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGmats(ispin,jspin,iorb,jorb,:) &
                        - (one+xi)*impGmats(ispin,ispin,iorb,iorb,:) &
                        - (one+xi)*impGmats(jspin,jspin,jorb,jorb,:))
                   !
                   impGreal(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGreal(ispin,jspin,iorb,jorb,:) &
                        - (one+xi)*impGreal(ispin,ispin,iorb,iorb,:) &
                        - (one+xi)*impGreal(jspin,jspin,jorb,jorb,:))
                enddo
             enddo
          enddo
       enddo
    end select
    !
    !
    !
    !
    ! if(ed_mode=="replica".AND.ed_para)then
    !    call SOC_jz_symmetrize(impGmats,dmft_bath%mask)
    !    call SOC_jz_symmetrize(impGreal,dmft_bath%mask)
    ! endif
    !
  end subroutine build_gf_nonsu2








  !PURPOSE: Evaluate the same orbital IORB, same spin ISPIN impurity GF.
  subroutine lanc_build_gf_nonsu2_diagOrb_diagSpin(iorb,ispin)
    integer      :: iorb,ispin,isite
    type(sector) :: sectorI,sectorJ
    integer :: ib(2*Ns)
    !
    isite = iorb+(ispin-1)*Ns
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
          if(ed_verbose>=3)then
             if(Jz_basis)then
                write(LOGfile,"(3(A,1F5.1,1X))")&
                     "add Jz:",Lzdiag(iorb)+Szdiag(ispin)/2.,&
                     "  from:",gettwoJz(isector)/2.,"  to:",gettwoJz(jsector)/2.
             else
                write(LOGfile,"(A15,2I6)")' From sector:',isector,sectorI%Ntot
             endif
          endif
       endif
       !
       !ADD ONE PARTICLE with IORB,ISPIN:
       if(Jz_basis)then
          jsector = getCDGsector_Jz(iorb,ispin,isector)
       else
          jsector = getCDGsector(1,ispin,isector)
       endif
       if(getN(isector)/=Nlevels.and.jsector>=0)then
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(Jz_basis)then
                if(ed_verbose>=3)write(LOGfile,"(A26,I3)")'apply c^+_a,s:',sectorJ%twoJz/2.
             else
                if(ed_verbose>=3)write(LOGfile,"(A26,I3)")'apply c^+_a,s:',sectorJ%Ntot
             endif
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             do i=1,sectorI%Dim
                call apply_op_CDG(i,j,sgn,iorb,1,ispin,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*state_cvec(i)
             enddo
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=zero
          endif
          !
          call tridiag_Hv_sector_nonsu2(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_nonsu2(one*norm2,state_e,alfa_,beta_,1,iorb,iorb,ispin,ispin)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
       endif
       !
       !REMOVE ONE PARTICLE with ISPIN:
       !
       if(Jz_basis)then
          jsector = getCsector_Jz(iorb,ispin,isector)
       else
          jsector = getCsector(1,ispin,isector)
       endif
       !
       if(getN(isector)/=0.and.jsector>=0)then
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(Jz_basis)then
                if(ed_verbose>=3)write(LOGfile,"(A26,I3)")'apply c_a,s:',sectorJ%twoJz/2.
             else
                if(ed_verbose>=3)write(LOGfile,"(A26,I3)")'apply c_a,s:',sectorJ%Ntot
             endif
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             do i=1,sectorI%Dim
                call apply_op_C(i,j,sgn,iorb,1,ispin,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*state_cvec(i)
             enddo
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=zero
          endif
          !
          call tridiag_Hv_sector_nonsu2(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_nonsu2(one*norm2,state_e,alfa_,beta_,-1,iorb,iorb,ispin,ispin)
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
    return
  end subroutine lanc_build_gf_nonsu2_diagOrb_diagSpin




  !PURPOSE: Evaluate the  different orbital IORB,JORB, different spin ISPIN,JSPIN impurity GF.
  subroutine lanc_build_gf_nonsu2_mixOrb_mixSpin(iorb,jorb,ispin,jspin)
    integer      :: iorb,jorb,ispin,jspin
    type(sector) :: sectorI,sectorJ
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
          if(ed_verbose>=3)then
             if(Jz_basis)then
                write(LOGfile,"(3(A,1F5.1,1X))")&
                     "add Jz:",Lzdiag(iorb)+Szdiag(ispin)/2.,&
                     "  from:",gettwoJz(isector)/2.,"  to:",gettwoJz(jsector)/2.
             else
                write(LOGfile,"(A15,2I6)")' From sector:',isector,sectorI%Ntot
             endif
          endif
       endif
       !
       !
       !APPLY (c^+_{jorb,jspin} + c^+_{iorb,ispin})|gs>
       if(Jz_basis)then
          jsector = getCDGsector_Jz(iorb,ispin,isector)
          ksector = getCDGsector_Jz(jorb,jspin,isector)
          if(getDim(jsector)/=getDim(ksector))stop "lanczos builgf dimensional error"
       else
          jsector = getCDGsector(1,ispin,isector)
       endif
       if(getN(isector)/=Nlevels.and.jsector>=0)then
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(Jz_basis)then
                if(ed_verbose>=3)write(LOGfile,"(A26,I3)")'apply c^+_b,s + c^+_a,s:',sectorJ%twoJz/2.
             else
                if(ed_verbose>=3)write(LOGfile,"(A26,I3)")'apply c^+_b,s + c^+_a,s:',sectorJ%Ntot
             endif
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             do i=1,sectorI%Dim
                call apply_op_CDG(i,j,sgn,iorb,1,ispin,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*state_cvec(i)
             enddo
             do i=1,sectorI%Dim
                call apply_op_CDG(i,j,sgn,jorb,1,jspin,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = vvinit(j) + sgn*state_cvec(i)
             enddo
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=zero
          endif
          !
          call tridiag_Hv_sector_nonsu2(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_nonsu2(one*norm2,state_e,alfa_,beta_,1,iorb,jorb,ispin,jspin)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
       endif
       !
       !
       !APPLY (c_{jorb,jspin} + c_{iorb,ispin})|gs>
       if(Jz_basis)then
          jsector = getCsector_Jz(iorb,ispin,isector)
          ksector = getCsector_Jz(jorb,jspin,isector)
          if(getDim(jsector)/=getDim(ksector))stop "lanczos builgf dimensional error"
       else
          jsector = getCsector(1,ispin,isector)
       endif
       if(getN(isector)/=0.and.jsector>=0)then
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(Jz_basis)then
                if(ed_verbose>=3)write(LOGfile,"(A26,I3)")'apply c_b,s + c_a,s:',sectorJ%twoJz/2.
             else
                if(ed_verbose>=3)write(LOGfile,"(A26,I3)")'apply c_b,s + c_a,s:',sectorJ%Ntot
             endif
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             do i=1,sectorI%Dim
                call apply_op_C(i,j,sgn,iorb,1,ispin,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*state_cvec(i)
             enddo
             do i=1,sectorI%Dim
                call apply_op_C(i,j,sgn,jorb,1,jspin,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = vvinit(j) + sgn*state_cvec(i)
             enddo
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=zero
          endif
          !
          call tridiag_Hv_sector_nonsu2(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_nonsu2(one*norm2,state_e,alfa_,beta_,-1,iorb,jorb,ispin,jspin)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
       endif
       !
       !
       !APPLY (+i*c^+_{jorb,jspin} + c^+_{iorb,ispin})|gs>
       if(Jz_basis)then
          jsector = getCDGsector_Jz(iorb,ispin,isector)
          ksector = getCDGsector_Jz(jorb,jspin,isector)
          if(getDim(jsector)/=getDim(ksector))stop "lanczos builgf dimensional error"
       else
          jsector = getCDGsector(1,ispin,isector)
       endif
       !
       if(getN(isector)/=Nlevels.and.jsector>=0)then
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(Jz_basis)then
                if(ed_verbose>=3)write(LOGfile,"(A26,I3)")'apply i*c^+_b,s + c^+_a,s:',sectorJ%twoJz/2.
             else
                if(ed_verbose>=3)write(LOGfile,"(A26,I3)")'apply i*c^+_b,s + c^+_a,s:',sectorJ%Ntot
             endif
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             do i=1,sectorI%Dim
                call apply_op_CDG(i,j,sgn,iorb,1,ispin,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*state_cvec(i)
             enddo
             do i=1,sectorI%Dim
                call apply_op_CDG(i,j,sgn,jorb,1,jspin,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = vvinit(j) + xi*sgn*state_cvec(i)
             enddo
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=zero
          endif
          !
          call tridiag_Hv_sector_nonsu2(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_nonsu2(xi*norm2,state_e,alfa_,beta_,1,iorb,jorb,ispin,jspin)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
       endif
       !
       !
       !APPLY (-xi*c_{jorb,jspin} + c_{iorb,ispin})|gs>
       if(Jz_basis)then
          jsector = getCsector_Jz(iorb,ispin,isector)
          ksector = getCsector_Jz(jorb,jspin,isector)
          if(getDim(jsector)/=getDim(ksector))stop "lanczos builgf dimensional error"
       else
          jsector = getCsector(1,ispin,isector)
       endif
       if(getN(isector)/=0.and.jsector>=0)then
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(Jz_basis)then
                if(ed_verbose>=3)write(LOGfile,"(A26,I3)")'apply -i*c_b,s + c_a,s:',sectorJ%twoJz/2.
             else
                if(ed_verbose>=3)write(LOGfile,"(A26,I3)")'apply -i*c_b,s + c_a,s:',sectorJ%Ntot
             endif
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             do i=1,sectorI%Dim
                call apply_op_C(i,j,sgn,iorb,1,ispin,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*state_cvec(i)
             enddo
             do i=1,sectorI%Dim
                call apply_op_C(i,j,sgn,jorb,1,jspin,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = vvinit(j) - xi*sgn*state_cvec(i)
             enddo
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=zero
          endif
          !
          call tridiag_Hv_sector_nonsu2(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_nonsu2(xi*norm2,state_e,alfa_,beta_,-1,iorb,jorb,ispin,jspin)
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
    return
  end subroutine lanc_build_gf_nonsu2_mixOrb_mixSpin




  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################






  subroutine add_to_lanczos_gf_nonsu2(vnorm2,Ei,alanc,blanc,isign,iorb,jorb,ispin,jspin)
    complex(8)                                 :: vnorm2,pesoBZ,peso
    real(8)                                    :: Ei,Egs,de
    integer                                    :: nlanc,itype
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    integer                                    :: isign,iorb,jorb,ispin,jspin
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
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,alanc)
       call Bcast_MPI(MpiComm,blanc)
    endif
#endif
    !pesoBZ = vnorm2/zeta_function
    !if(finiteT)pesoBZ = vnorm2*exp(-beta*(Ei-Egs))/zeta_function
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
          impGmats(ispin,jspin,iorb,jorb,i)=impGmats(ispin,jspin,iorb,jorb,i) + peso/(iw-isign*de)
       enddo
       do i=1,Lreal
          iw=dcmplx(wr(i),eps)
          impGreal(ispin,jspin,iorb,jorb,i)=impGreal(ispin,jspin,iorb,jorb,i) + peso/(iw-isign*de)
       enddo
    enddo
  end subroutine add_to_lanczos_gf_nonsu2








  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################









  !+------------------------------------------------------------------+
  !PURPOSE  : Build the Self-energy functions, NONSU2 case
  !+------------------------------------------------------------------+
  subroutine build_sigma_nonsu2
    integer                                           :: i,j,isign,unit(7),iorb,jorb,ispin,jspin,io,jo
    complex(8)                                        :: fg0
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: invG0mats,invimpG0mats,impDeltamats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal) :: invG0real,invimpG0real,impDeltareal
    complex(8),dimension(Nspin*Norb,Nspin*Norb)       :: invGimp,Foo
    character(len=20)                                 :: suffix
    !
    !
    impG0mats=zero
    impG0real=zero
    invG0mats = zero
    invG0real = zero
    !
    !Get G0^-1
    invG0mats(:,:,:,:,:)=invg0_bath_function(dcmplx(0d0,wm(:)),dmft_bath)
    invG0real(:,:,:,:,:)=invg0_bath_function(dcmplx(wr(:),eps),dmft_bath)
    !Get impDelta_anderson
    impDeltamats(:,:,:,:,:)=delta_bath_function(dcmplx(0d0,wm(:)),dmft_bath)
    impDeltareal(:,:,:,:,:)=delta_bath_function(dcmplx(wr(:),eps),dmft_bath)
    !Get inverse functions
    invimpG0mats=invG0mats
    invimpG0real=invG0real
    !
    select case(bath_type)
       !
    case ("normal")
       !
       !Get Gimp^-1 - Matsubara freq.
       do i=1,Lmats
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      if (iorb.eq.jorb) then
                         io = iorb + (ispin-1)*Norb
                         jo = jorb + (jspin-1)*Norb
                         invGimp(io,jo) = impGmats(ispin,jspin,iorb,jorb,i)
                      endif
                   enddo
                enddo
             enddo
          enddo
          call inv(invGimp)!<--- get [G_{imp}]^-1
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      if (iorb.eq.jorb) then
                         io = iorb + (ispin-1)*Norb
                         jo = jorb + (jspin-1)*Norb
                         impSmats(ispin,jspin,iorb,jorb,i) = invG0mats(ispin,jspin,iorb,jorb,i) - invGimp(io,jo) !<-- calG0_imp^-1 - Gimp^-1
                      endif
                   enddo
                enddo
             enddo
          enddo
       enddo
       !Get Gimp^-1 - Real freq.
       do i=1,Lreal
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      if (iorb.eq.jorb) then
                         io = iorb + (ispin-1)*Norb
                         jo = jorb + (jspin-1)*Norb
                         invGimp(io,jo) = impGreal(ispin,jspin,iorb,jorb,i)
                      endif
                   enddo
                enddo
             enddo
          enddo
          call inv(invGimp)!<--- get [G_{imp}]^-1
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      if (iorb.eq.jorb) then
                         io = iorb + (ispin-1)*Norb
                         jo = jorb + (jspin-1)*Norb
                         impSreal(ispin,jspin,iorb,jorb,i) = invG0real(ispin,jspin,iorb,jorb,i) - invGimp(io,jo) !<-- calG0_imp^-1 - Gimp^-1
                      endif
                   enddo
                enddo
             enddo
          enddo
       enddo
       !
    case ("hybrid","replica")
       !
       !Get Gimp^-1 - Matsubara freq.
       do i=1,Lmats
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb + (ispin-1)*Norb
                      jo = jorb + (jspin-1)*Norb
                      invGimp(io,jo) = impGmats(ispin,jspin,iorb,jorb,i)
                   enddo
                enddo
             enddo
          enddo
          call inv(invGimp)!<--- get [G_{imp}]^-1
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb + (ispin-1)*Norb
                      jo = jorb + (jspin-1)*Norb
                      impSmats(ispin,jspin,iorb,jorb,i) = invG0mats(ispin,jspin,iorb,jorb,i) - invGimp(io,jo) !<-- calG0_imp^-1 - Gimp^-1
                   enddo
                enddo
             enddo
          enddo
       enddo
       !Get Gimp^-1 - Real freq.
       do i=1,Lreal
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb + (ispin-1)*Norb
                      jo = jorb + (jspin-1)*Norb
                      invGimp(io,jo) = impGreal(ispin,jspin,iorb,jorb,i)
                   enddo
                enddo
             enddo
          enddo
          call inv(invGimp)!<--- get [G_{imp}]^-1
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb + (ispin-1)*Norb
                      jo = jorb + (jspin-1)*Norb
                      impSreal(ispin,jspin,iorb,jorb,i) = invG0real(ispin,jspin,iorb,jorb,i) - invGimp(io,jo) !<-- calG0_imp^-1 - Gimp^-1
                   enddo
                enddo
             enddo
          enddo
       enddo
       !
    end select
    !
    !Get G0and:
    impG0mats(:,:,:,:,:) = g0and_bath_function(dcmplx(0d0,wm(:)),dmft_bath)
    impG0real(:,:,:,:,:) = g0and_bath_function(dcmplx(wr(:),eps),dmft_bath)
    !
  end subroutine build_sigma_nonsu2





END MODULE ED_GF_NONSU2
