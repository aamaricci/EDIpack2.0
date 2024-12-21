MODULE ED_RDM_SUPERC
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
  USE ED_HAMILTONIAN_SUPERC
  implicit none
  private
  !
  public                           :: imp_rdm_superc


  real(8)                          :: Egs ! Ground-state energy
  real(8)                          :: Ei
  integer                          :: iorb,jorb,iorb1,jorb1
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
  complex(8),dimension(:),allocatable :: state_cvec
  logical                          :: Jcondition
  !
  integer                          :: iImpUp,iImpDw,iImp
  integer                          :: jImpUp,jImpDw,jImp
  integer                          :: ib,iBath

  integer                          :: LenBath
  integer,allocatable              :: Bath(:)

  type(sector)                     :: sectorI,sectorJ
  character(len=128)               :: fmt

contains 

  !+-------------------------------------------------------------------+
  !PURPOSE  : Lanc method
  !+-------------------------------------------------------------------+
  subroutine imp_rdm_superc()
    integer                         :: istate,val
    real(8),dimension(Norb)         :: nup,ndw

    !
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG imp_rdm_superc"
#endif
    Egs     = state_list%emin
    !
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")&
         "DEBUG imp_rdm_superc: get imp_rdm"
#endif
    !
    !IMPURITY DENSITY MATRIX
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")&
         "DEBUG observables_superc: eval impurity density matrix \rho_IMP = Tr_BATH(\rho)"
#endif
    impurity_density_matrix=zero
    do istate=1,state_list%size
#ifdef _DEBUG
       if(ed_verbose>3)write(Logfile,"(A)")&
            "DEBUG observables_superc: get contribution from state:"//str(istate)
#endif

       isector = es_return_sector(state_list,istate)
       Ei      = es_return_energy(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          call es_return_cvector(MpiComm,state_list,istate,state_cvec) 
       else
          call es_return_cvector(state_list,istate,state_cvec) 
       endif
#else
       call es_return_cvector(state_list,istate,state_cvec) 
#endif
       !
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       if(MpiMaster)then
          call build_sector(isector,sectorI,itrace=.true.)
          !
          do IimpUp=0,2**Norb-1
             do IimpDw=0,2**Norb-1
                do JimpUp=0,2**Norb-1
                   do JimpDw=0,2**Norb-1
                      !
                      !Build indices of the RDM in 1:4**Norb
                      iImp = iImpUp + iImpDw*2**Norb
                      jImp = jImpUp + jImpDw*2**Norb
                      call sp_return_intersection(sectorI%H(1)%sp,iImp,jImp,Bath,lenBath)
                      if(lenBATH==0)cycle
                      !
                      !=== >>> TRACE over bath states <<< =================================================
                      do ib=1,lenBath
                         iBath = Bath(ib)
                         !get state I
                         i = binary_search(sectorI%H(1)%map,iImpUp + iImpDw*2**Ns + iBath*2**Norb)
                         !
                         !get state J
                         j = binary_search(sectorI%H(1)%map,jImpUp + jImpDw*2**Ns + iBath*2**Norb)
                         !
                         !Build (i,j)_th contribution to the (io,jo)_th element of \rho_IMP
                         io = (iImpUp+1) + 2**Norb*iImpDw
                         jo = (jImpUp+1) + 2**Norb*jImpDw
                         impurity_density_matrix(io,jo) = impurity_density_matrix(io,jo) + &
                              state_cvec(i)*conjg(state_cvec(j))*peso
                         !-----------------------------------------------------------------
                      enddo
                   enddo !=============================================================================
                   !
                enddo
             enddo
             !
          enddo
          !
          call delete_sector(sectorI)
       endif
       !
       if(allocated(state_cvec))deallocate(state_cvec)
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
  end subroutine imp_rdm_superc






end MODULE ED_RDM_SUPERC

















