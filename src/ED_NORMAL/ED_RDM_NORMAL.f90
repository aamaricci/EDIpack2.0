MODULE ED_RDM_NORMAL
  USE SF_CONSTANTS, only:zero,pi,xi
  USE SF_IOTOOLS, only:free_unit,reg,txtfy
  USE SF_ARRAYS, only: arange
  USE SF_LINALG
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_EIGENSPACE
  USE ED_SETUP
  USE ED_SECTOR
  USE ED_BATH
  USE ED_HAMILTONIAN_NORMAL
  implicit none
  private
  !
  public                           :: imp_rdm_normal

  real(8)                          :: Egs ! Ground-state energy
  real(8)                          :: Ei
  integer                          :: iorb,jorb,iorb1,jorb1
  integer                          :: ispin,jspin
  integer                          :: isite,jsite
  integer                          :: ibath
  integer                          :: r,m,k,k1,k2,k3,k4
  integer                          :: iup,idw
  integer                          :: jup,jdw
  integer                          :: mup,mdw
  integer                          :: iph,i_el,isectorDim
  real(8)                          :: sgn,sgn1,sgn2,sg1,sg2,sg3,sg4
  real(8)                          :: gs_weight
  !
  real(8)                          :: peso
  real(8)                          :: norm
  !
  integer                          :: i,j,ii,io,jo
  integer                          :: isector,jsector
  !
  real(8),dimension(:),allocatable :: vvinit
  real(8),dimension(:),allocatable :: state_dvec
  logical                          :: Jcondition
  !
  integer                          :: iImpUp,iImpDw
  integer                          :: jImpUp,jImpDw
  integer                          :: iBathUp,iBathDw
  integer                          :: iBup,iBdw

  integer                          :: LenBathUp,LenBathDw
  integer,allocatable              :: BathUp(:),BathDw(:)

  type(sector)                     :: sectorI,sectorJ
  character(len=128)               :: fmt

contains 

  !+-------------------------------------------------------------------+
  !PURPOSE  : Lanc method
  !+-------------------------------------------------------------------+
  subroutine imp_rdm_normal()
    integer                         :: istate,val
    real(8),dimension(Norb)         :: nup,ndw

    !
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG imp_rdm_normal"
#endif
    Egs     = state_list%emin
    !
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")&
         "DEBUG imp_rdm_normal: get imp_rdm"
#endif
    !
    !IMPURITY DENSITY MATRIX
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")&
         "DEBUG observables_normal: eval impurity density matrix \rho_IMP = Tr_BATH(\rho)"
#endif
    impurity_density_matrix=zero
    do istate=1,state_list%size
#ifdef _DEBUG
       if(ed_verbose>3)write(Logfile,"(A)")&
            "DEBUG observables_normal: get contribution from state:"//str(istate)
#endif

       isector = es_return_sector(state_list,istate)
       Ei      = es_return_energy(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          call es_return_dvector(MpiComm,state_list,istate, state_dvec)
       else
          call es_return_dvector(state_list,istate, state_dvec)
       endif
#else
       call es_return_dvector(state_list,istate, state_dvec)
#endif
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       if(MpiMaster)then
          call build_sector(isector,sectorI,itrace=.true.)
          !
          do IimpUp=0,2**Norb-1
             do JimpUp=0,2**Norb-1
                !
                !Finding the unique bath states connecting IimpUp and JimpUp -> BATHup(:)
                call sp_return_intersection(sectorI%H(1)%sp,IimpUp,JimpUp,BATHup,lenBATHup)
                if(lenBATHup==0)cycle  !there are no bath states intersecting IimpUp,JimpUp
                !
                do IimpDw=0,2**Norb-1
                   do JimpDw=0,2**Norb-1
                      !
                      !Finding the unique bath states connecting IimpDw and JimpDw -> BATHdw(:)
                      call sp_return_intersection(sectorI%H(2)%sp,IimpDw,JimpDw,BATHdw,lenBATHdw)
                      if(lenBATHdw==0)cycle  !there are no bath states intersecting IimpDw,JimpDw
                      !
                      !=== >>> TRACE over bath states <<< =================================================
                      do ibUP=1,lenBATHup
                         IbathUp = BATHup(ibUP)
                         do ibDW=1,lenBATHdw
                            IbathDw = BATHdw(ibDW)
                            !-----------------------------------------------------
                            !Allowed spin Fock space Istates:
                            !Iup = IimpUp +  2^Norb * IbathUp
                            !Idw = IimpDw +  2^Norb * IbathDw
                            !
                            !Corresponding sector indices per spin:
                            iUP= binary_search(sectorI%H(1)%map,IimpUp + 2**Norb*IbathUp)
                            iDW= binary_search(sectorI%H(2)%map,IimpDw + 2**Norb*IbathDw)
                            !
                            !Global sector index:
                            i  = iUP + (iDW-1)*sectorI%DimUp
                            !-----------------------------------------------------
                            !Allowed spin Fock space Jstates:
                            !Jup = JimpUp +  2^Norb * IbathUp
                            !Jdw = JimpDw +  2^Norb * IbathDw
                            !
                            !Corresponding sector jndices per spin:
                            jUP= binary_search(sectorI%H(1)%map,JimpUp + 2**Norb*IbathUp)
                            jDW= binary_search(sectorI%H(2)%map,JimpDw + 2**Norb*IbathDw)
                            !
                            !Global sector jndex:
                            j  = jUP + (jDW-1)*sectorI%DimUp
                            !-----------------------------------------------------
                            !Final [Up,Dw] composition for the impurity:
                            io = (IimpUp + 2**Norb*IimpDw) + 1
                            jo = (JimpUp + 2**Norb*JimpDw) + 1
                            !-----------------------------------------------------------------
                            !(i,j)_th contribution to the (io,jo)_th element of \rho_IMP
                            impurity_density_matrix(io,jo) = impurity_density_matrix(io,jo) + &
                                 state_dvec(i)*state_dvec(j)*peso
                            !-----------------------------------------------------------------
                         enddo
                      enddo !=============================================================================
                      !
                   enddo
                enddo
                !
             enddo
          enddo
          !
          call delete_sector(sectorI)
       endif
       !
       if(allocated(state_dvec))deallocate(state_dvec)
       !
    enddo
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")""
#endif
    !Useful prints to test:
    if(MPIMASTER)then
       if(Norb<=3)then
          write(LOGfile,*)"RDM:"
          write(fmt,"(A1,I0,A2,I0,A1,I0,A5)")"(",4**Norb,"(F5.2,1x))"
          do i=1,4**Norb
             write(LOGfile,fmt)(dreal(impurity_density_matrix(i,j)),j=1,4**Norb)
          enddo
#ifdef _DEBUG
          if(Norb==1)then
             ! Cfr Eq. 4 in Mod Phys Lett B 2013 27:05
             write(LOGfile,*)1-ed_dens_up(1)-ed_dens_dw(1)+ed_docc(1),abs(1-ed_dens_up(1)-ed_dens_dw(1)+ed_docc(1)-impurity_density_matrix(1,1))
             write(LOGfile,*)ed_dens_up(1)-ed_docc(1),abs(ed_dens_up(1)-ed_docc(1)-impurity_density_matrix(2,2))
             write(LOGfile,*)ed_dens_dw(1)-ed_docc(1),abs(ed_dens_dw(1)-ed_docc(1)-impurity_density_matrix(3,3))
             write(LOGfile,*)ed_docc(1),abs(ed_docc(1)-impurity_density_matrix(4,4))
          endif
#endif
       endif
    endif
    !
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")""
#endif
  end subroutine imp_rdm_normal






end MODULE ED_RDM_NORMAL

















