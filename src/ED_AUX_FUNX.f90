MODULE ED_AUX_FUNX
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE SF_TIMER
  USE SF_LINALG, only: eye
  USE SF_PAULI
  USE SF_MISC, only: assert_shape
  USE SF_IOTOOLS, only:free_unit,reg,txtfy
  USE SF_ARRAYS,  only: arange,linspace
  implicit none
  private


  interface lso2nnn_reshape
     module procedure d_nlso2nnn
     module procedure c_nlso2nnn
  end interface lso2nnn_reshape

  interface so2nn_reshape
     module procedure d_nso2nn
     module procedure c_nso2nn
  end interface so2nn_reshape

  interface nnn2lso_reshape
     module procedure d_nnn2nlso
     module procedure c_nnn2nlso
  end interface nnn2lso_reshape

  interface nn2so_reshape
     module procedure d_nn2nso
     module procedure c_nn2nso
  end interface nn2so_reshape

  interface print_state_vector
     module procedure print_state_vector_ivec
     module procedure print_state_vector_ivec_ud
     module procedure print_state_vector_int
  end interface print_state_vector


  interface print_Hloc
     module procedure print_Hloc_so_c
     module procedure print_Hloc_nn_c
  end interface print_Hloc


  interface ed_set_suffix
     module procedure :: ed_set_suffix_i
     module procedure :: ed_set_suffix_d
     module procedure :: ed_set_suffix_c
  end interface ed_set_suffix

  !   !This is to overload read/write procedure to GFmatrix derived types
  ! #if __GNUC__ > 6
  !   interface read(formatted)
  !      procedure read_formatted
  !   end interface read(formatted)

  !   interface write(formatted)
  !      procedure write_formatted
  !   end interface write(formatted)
  ! #endif


  !FERMIONIC OPERATORS IN BITWISE OPERATIONS
  public :: c,cdg
  !AUX BIT OPERATIONS
  public :: bdecomp
  public :: breorder
  public :: bjoin
  !BINARY SEARCH
  public :: binary_search
  !AUX RESHAPE FUNCTIONS (internal use)
  public :: index_stride_so
  public :: lso2nnn_reshape
  public :: so2nn_reshape
  public :: nnn2lso_reshape
  public :: nn2so_reshape
  public :: so2os_reshape
  public :: os2so_reshape
  !SEARCH CHEMICAL POTENTIAL, this should go into DMFT_TOOLS I GUESS
  public :: ed_search_variable
  public :: search_chemical_potential
  !SOC RELATED STUFF
  public :: SOC_jz_symmetrize
  public :: atomic_SOC
  public :: atomic_SOC_rotation
  public :: orbital_Lz_rotation_NorbNspin
  public :: orbital_Lz_rotation_Norb
  public :: atomic_j
  !ALLOCATE/DEALLOCATE GRIDS
  public :: allocate_grids
  public :: deallocate_grids
  !PRINT STATE VECTORS
  public :: print_state_vector
  !PRINT LOCAL HAMILTONIAN
  public :: print_hloc
  !SET/RESET GLOBAL FILE SUFFIX
  public :: ed_set_suffix
  public :: ed_reset_suffix
  !MPI PROCEDURES
#ifdef _MPI
  interface scatter_vector_MPI
     module procedure :: d_scatter_vector_MPI
     module procedure :: c_scatter_vector_MPI
  end interface scatter_vector_MPI

  interface scatter_basis_MPI
     module procedure :: d_scatter_basis_MPI
     module procedure :: c_scatter_basis_MPI
  end interface scatter_basis_MPI

  interface gather_vector_MPI
     module procedure :: d_gather_vector_MPI
     module procedure :: c_gather_vector_MPI
  end interface gather_vector_MPI

  interface allgather_vector_MPI
     module procedure :: d_allgather_vector_MPI
     module procedure :: c_allgather_vector_MPI
  end interface allgather_vector_MPI
  public :: scatter_vector_MPI
  public :: scatter_basis_MPI
  public :: gather_vector_MPI
  public :: allgather_vector_MPI
#endif


contains


  subroutine ed_reset_suffix()
    ed_file_suffix=''
  end subroutine ed_reset_suffix

  subroutine ed_set_suffix_i(indx)
    integer :: indx
    ed_file_suffix=reg(ineq_site_suffix)//str(indx,site_indx_padding)
  end subroutine ed_set_suffix_i
  subroutine ed_set_suffix_d(indx)
    real(8) :: indx
    ed_file_suffix=reg(ineq_site_suffix)//str(indx)
  end subroutine ed_set_suffix_d
  subroutine ed_set_suffix_c(indx)
    character(len=*) :: indx
    ed_file_suffix=reg(ineq_site_suffix)//str(indx)
  end subroutine ed_set_suffix_c






  !##################################################################
  !##################################################################
  !               CREATION / DESTRUCTION OPERATORS
  !##################################################################
  !##################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE: input state |in> of the basis and calculates 
  !   |out>=C_pos|in>  OR  |out>=C^+_pos|in> ; 
  !   the sign of |out> has the phase convention, pos labels the sites
  !+-------------------------------------------------------------------+
  subroutine c(pos,in,out,fsgn)
    integer,intent(in)    :: pos
    integer,intent(in)    :: in
    integer,intent(inout) :: out
    real(8),intent(inout) :: fsgn    
    integer               :: l
    if(.not.btest(in,pos-1))stop "C error: C_i|...0_i...>"
    fsgn=1d0
    do l=1,pos-1
       if(btest(in,l-1))fsgn=-fsgn
    enddo
    out = ibclr(in,pos-1)
  end subroutine c

  subroutine cdg(pos,in,out,fsgn)
    integer,intent(in)    :: pos
    integer,intent(in)    :: in
    integer,intent(inout) :: out
    real(8),intent(inout) :: fsgn    
    integer               :: l
    if(btest(in,pos-1))stop "C^+ error: C^+_i|...1_i...>"
    fsgn=1d0
    do l=1,pos-1
       if(btest(in,l-1))fsgn=-fsgn
    enddo
    out = ibset(in,pos-1)
  end subroutine cdg



  !##################################################################
  !##################################################################
  !               BITWISE OPERATIONS
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE  : input a state |i> and output a vector ivec(Nlevels)
  !with its binary decomposition
  !(corresponds to the decomposition of the number i-1)
  !+------------------------------------------------------------------+
  function bdecomp(i,Ntot) result(ivec)
    integer :: Ntot,ivec(Ntot),l,i
    logical :: busy
    !this is the configuration vector |1,..,Ns,Ns+1,...,Ntot>
    !obtained from binary decomposition of the state/number i\in 2^Ntot
    do l=0,Ntot-1
       busy=btest(i,l)
       ivec(l+1)=0
       if(busy)ivec(l+1)=1
    enddo
  end function bdecomp


  !+------------------------------------------------------------------+
  ! Reorder a binary decomposition so to have a state of the form:
  ! default: |(1:Norb),([1:Nbath]_1, [1:Nbath]_2, ... ,[1:Nbath]_Norb)>_spin
  ! hybrid:  |(1:Norb),([1:Nbath])_spin
  ! replica: |(1:Norb),([1:Norb]_1, [1:Norb]_2, ...  , [1:Norb]_Nbath)>_spin
  !
  !> case (ed_total_ud):
  !   (T): Ns_Ud=1, Ns_Orb=Ns.
  !        bdecomp is already of the form above [1:Ns]
  !   (F): Ns_Ud=Norb, Ns_Orb=Ns/Norb==1+Nbath
  !        bdecomp is
  !        |( [1:1+Nbath]_1,...,[1:1+Nbath]_Norb)>_spin
  !+------------------------------------------------------------------+
  function breorder(Nins) result(Ivec)
    integer,intent(in),dimension(Ns_Ud,Ns_Orb) :: Nins ![1,Ns] - [Norb,1+Nbath]
    integer,dimension(Ns)                      :: Ivec ![Ns]
    integer                                    :: iud,ibath,indx
    select case (ed_total_ud)
    case (.true.)
       Ivec = Nins(1,:)
    case (.false.)
       do iud=1,Ns_Ud           ![1:Norb]
          Ivec(iud) = Nins(iud,1)
          do ibath=1,Nbath
             indx = getBathStride(iud,ibath) !take care of normal/
             Ivec(indx) = Nins(iud,1+ibath)
          enddo
       enddo
    end select
  end function breorder


  !+------------------------------------------------------------------+
  !PURPOSE  : input a vector ib(Nlevels) with the binary sequence 
  ! and output the corresponding state |i>
  !(corresponds to the recomposition of the number i-1)
  !+------------------------------------------------------------------+
  function bjoin(ib,Ntot) result(i)
    integer                 :: Ntot
    integer,dimension(Ntot) :: ib
    integer                 :: i,j
    i=0
    do j=0,Ntot-1
       i=i+ib(j+1)*2**j
    enddo
  end function bjoin


  !+------------------------------------------------------------------+
  !PURPOSE : binary search of a value in an array
  !+------------------------------------------------------------------+
  recursive function binary_search(a,value) result(bsresult)
    integer,intent(in) :: a(:), value
    integer            :: bsresult, mid
    mid = size(a)/2 + 1
    if (size(a) == 0) then
       bsresult = 0        ! not found
       !stop "binary_search error: value not found"
    else if (a(mid) > value) then
       bsresult= binary_search(a(:mid-1), value)
    else if (a(mid) < value) then
       bsresult = binary_search(a(mid+1:), value)
       if (bsresult /= 0) then
          bsresult = mid + bsresult
       end if
    else
       bsresult = mid      ! SUCCESS!!
    end if
  end function binary_search



















  !##################################################################
  !##################################################################
  !AUXILIARY COMPUTATIONAL ROUTINES ARE HERE BELOW:
  !##################################################################
  !##################################################################
#ifdef _MPI
  !! Scatter V into the arrays Vloc on each thread: sum_threads(size(Vloc)) must be equal to size(v)
  subroutine d_scatter_vector_MPI(MpiComm,v,vloc)
    integer                          :: MpiComm
    real(8),dimension(:)             :: v    !size[N]
    real(8),dimension(:)             :: vloc !size[Nloc]
    integer                          :: i,iph,irank,Nloc,N
    integer                          :: v_start,v_end,vloc_start,vloc_end
    integer,dimension(:),allocatable :: Counts,Offset
    integer                          :: MpiSize,MpiIerr
    logical                          :: MpiMaster
    !
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A)")"DEBUG d_scatter_vector_MPI: scatter v into vloc"
#endif
    !
    if( MpiComm == MPI_UNDEFINED .OR. MpiComm == Mpi_Comm_Null )return
    ! stop "scatter_vector_MPI error: MpiComm == MPI_UNDEFINED"
    !
    MpiSize   = get_size_MPI(MpiComm)
    MpiMaster = get_master_MPI(MpiComm)
    !
    Nloc = size(Vloc)
    N = 0
    call AllReduce_MPI(MpiComm,Nloc,N)
    if(MpiMaster.AND.N /= size(V)) stop "scatter_vector_MPI error: size(V) != Mpi_Allreduce(Nloc)"
    !
    allocate(Counts(0:MpiSize-1)) ; Counts=0
    allocate(Offset(0:MpiSize-1)) ; Offset=0
    !
    !Get Counts;
    call MPI_AllGather(Nloc/DimPh,1,MPI_INTEGER,Counts,1,MPI_INTEGER,MpiComm,MpiIerr)
    !
    !Get Offset:
    Offset(0)=0
    do i=1,MpiSize-1
       Offset(i) = Offset(i-1) + Counts(i-1)
    enddo
    !
    Vloc=0d0
    do iph=1,Dimph
       if(MpiMaster)then
          v_start = 1 + (iph-1)*(N/Dimph)
          v_end = iph*(N/Dimph)
       else
          v_start = 1
          v_end = 1
       endif
       vloc_start = 1 + (iph-1)*(Nloc/Dimph)
       vloc_end = iph*(Nloc/Dimph)
       call MPI_Scatterv(V(v_start:v_end),Counts,Offset,MPI_DOUBLE_PRECISION,&
            Vloc(vloc_start:vloc_end),Nloc/DimPh,MPI_DOUBLE_PRECISION,0,MpiComm,MpiIerr)
    enddo
    !
    return
  end subroutine d_scatter_vector_MPI

  subroutine c_scatter_vector_MPI(MpiComm,v,vloc)
    integer                          :: MpiComm
    complex(8),dimension(:)          :: v    !size[N]
    complex(8),dimension(:)          :: vloc !size[Nloc]
    integer                          :: i,irank,Nloc,N
    integer,dimension(:),allocatable :: Counts,Offset
    integer                          :: MpiSize,MpiIerr
    logical                          :: MpiMaster
    !
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A)")"DEBUG c_scatter_vector_MPI: scatter v into vloc"
#endif
    !
    if( MpiComm == MPI_UNDEFINED ) stop "scatter_vector_MPI error: MpiComm == MPI_UNDEFINED"
    !
    MpiSize   = get_size_MPI(MpiComm)
    MpiMaster = get_master_MPI(MpiComm)
    !
    Nloc = size(Vloc)
    N = 0
    call AllReduce_MPI(MpiComm,Nloc,N)
    if(MpiMaster.AND.N /= size(V)) stop "scatter_vector_MPI error: size(V) != Mpi_Allreduce(Nloc)"
    !
    allocate(Counts(0:MpiSize-1)) ; Counts=0
    allocate(Offset(0:MpiSize-1)) ; Offset=0
    !
    !Get Counts;
    call MPI_AllGather(Nloc,1,MPI_INTEGER,Counts,1,MPI_INTEGER,MpiComm,MpiIerr)
    !
    !Get Offset:
    Offset(0)=0
    do i=1,MpiSize-1
       Offset(i) = Offset(i-1) + Counts(i-1)
    enddo
    !
    Vloc=0
    call MPI_Scatterv(V,Counts,Offset,MPI_DOUBLE_COMPLEX,Vloc,Nloc,MPI_DOUBLE_COMPLEX,0,MpiComm,MpiIerr)
    !
    return
  end subroutine c_scatter_vector_MPI


  subroutine d_scatter_basis_MPI(MpiComm,v,vloc)
    integer                :: MpiComm
    real(8),dimension(:,:) :: v    !size[N,N]
    real(8),dimension(:,:) :: vloc !size[Nloc,Neigen]
    integer                :: N,Nloc,Neigen,i
    !
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A)")"DEBUG d_scatter_basis_MPI: scatter many v"
#endif
    N      = size(v,1)
    Nloc   = size(vloc,1)
    Neigen = size(vloc,2)
    if( size(v,2) < Neigen ) stop "error scatter_basis_MPI: size(v,2) < Neigen"
    !
    do i=1,Neigen
       call scatter_vector_MPI(MpiComm,v(:,i),vloc(:,i))
    end do
    !
    return
  end subroutine d_scatter_basis_MPI

  subroutine c_scatter_basis_MPI(MpiComm,v,vloc)
    integer                   :: MpiComm
    complex(8),dimension(:,:) :: v    !size[N,N]
    complex(8),dimension(:,:) :: vloc !size[Nloc,Neigen]
    integer                   :: N,Nloc,Neigen,i
    !
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A)")"DEBUG c_scatter_basis_MPI: scatter many v"
#endif
    !
    N      = size(v,1)
    Nloc   = size(vloc,1)
    Neigen = size(vloc,2)
    if( size(v,2) < Neigen ) stop "error scatter_basis_MPI: size(v,2) < Neigen"
    !
    do i=1,Neigen
       call scatter_vector_MPI(MpiComm,v(:,i),vloc(:,i))
    end do
    !
    return
  end subroutine c_scatter_basis_MPI




  !! AllGather Vloc on each thread into the array V: sum_threads(size(Vloc)) must be equal to size(v)
  subroutine d_gather_vector_MPI(MpiComm,vloc,v)
    integer                          :: MpiComm
    real(8),dimension(:)             :: vloc !size[Nloc]
    real(8),dimension(:)             :: v    !size[N]
    integer                          :: i,iph,irank,Nloc,N
    integer                          :: v_start,v_end,vloc_start,vloc_end
    integer,dimension(:),allocatable :: Counts,Offset
    integer                          :: MpiSize,MpiIerr
    logical                          :: MpiMaster
    !
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A)")"DEBUG d_gather_basis_MPI: gather  v"
#endif
    !
    if(  MpiComm == MPI_UNDEFINED .OR. MpiComm == Mpi_Comm_Null ) return
    !stop "gather_vector_MPI error: MpiComm == MPI_UNDEFINED"
    !
    MpiSize   = get_size_MPI(MpiComm)
    MpiMaster = get_master_MPI(MpiComm)
    !
    Nloc = size(Vloc)
    N = 0
    call AllReduce_MPI(MpiComm,Nloc,N)
    if(MpiMaster.AND.N /= size(V)) stop "gather_vector_MPI error: size(V) != Mpi_Allreduce(Nloc)"
    !
    allocate(Counts(0:MpiSize-1)) ; Counts=0
    allocate(Offset(0:MpiSize-1)) ; Offset=0
    !
    !Get Counts;
    call MPI_AllGather(Nloc/Dimph,1,MPI_INTEGER,Counts,1,MPI_INTEGER,MpiComm,MpiIerr)
    !
    !Get Offset:
    Offset(0)=0
    do i=1,MpiSize-1
       Offset(i) = Offset(i-1) + Counts(i-1)
    enddo
    !
    do iph=1,Dimph
       if(MpiMaster)then
          v_start = 1 + (iph-1)*(N/Dimph)
          v_end = iph*(N/Dimph)
       else
          v_start = 1
          v_end = 1
       endif
       vloc_start = 1 + (iph-1)*(Nloc/Dimph)
       vloc_end = iph*(Nloc/Dimph)
       !
       call MPI_Gatherv(Vloc(vloc_start:vloc_end),Nloc/DimPh,MPI_DOUBLE_PRECISION,&
            V(v_start:v_end),Counts,Offset,MPI_DOUBLE_PRECISION,0,MpiComm,MpiIerr)
    enddo
    !
    return
  end subroutine d_gather_vector_MPI

  subroutine c_gather_vector_MPI(MpiComm,vloc,v)
    integer                          :: MpiComm
    complex(8),dimension(:)          :: vloc !size[Nloc]
    complex(8),dimension(:)          :: v    !size[N]
    integer                          :: i,irank,Nloc,N
    integer,dimension(:),allocatable :: Counts,Offset
    integer                          :: MpiSize,MpiIerr
    logical                          :: MpiMaster
    !
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A)")"DEBUG c_gather_basis_MPI: gather  v"
#endif
    !
    if( MpiComm == MPI_UNDEFINED ) stop "gather_vector_MPI error: MpiComm == MPI_UNDEFINED"
    !
    MpiSize   = get_size_MPI(MpiComm)
    MpiMaster = get_master_MPI(MpiComm)
    !
    Nloc = size(Vloc)
    N = 0
    call AllReduce_MPI(MpiComm,Nloc,N)
    if(MpiMaster.AND.N /= size(V)) stop "gather_vector_MPI error: size(V) != Mpi_Allreduce(Nloc)"
    !
    allocate(Counts(0:MpiSize-1)) ; Counts=0
    allocate(Offset(0:MpiSize-1)) ; Offset=0
    !
    !Get Counts;
    call MPI_AllGather(Nloc,1,MPI_INTEGER,Counts,1,MPI_INTEGER,MpiComm,MpiIerr)
    !
    !Get Offset:
    Offset(0)=0
    do i=1,MpiSize-1
       Offset(i) = Offset(i-1) + Counts(i-1)
    enddo
    !
    call MPI_Gatherv(Vloc,Nloc,MPI_DOUBLE_COMPLEX,V,Counts,Offset,MPI_DOUBLE_COMPLEX,0,MpiComm,MpiIerr)
    !
    return
  end subroutine c_gather_vector_MPI



  !! AllGather Vloc on each thread into the array V: sum_threads(size(Vloc)) must be equal to size(v)
  subroutine d_allgather_vector_MPI(MpiComm,vloc,v)
    integer                          :: MpiComm
    real(8),dimension(:)             :: vloc !size[Nloc]
    real(8),dimension(:)             :: v    !size[N]
    integer                          :: i,iph,irank,Nloc,N
    integer                          :: v_start,v_end,vloc_start,vloc_end
    integer,dimension(:),allocatable :: Counts,Offset
    integer                          :: MpiSize,MpiIerr
    logical                          :: MpiMaster
    !
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A)")"DEBUG d_allgather_basis_MPI: allgather v"
#endif
    !
    if(  MpiComm == MPI_UNDEFINED .OR. MpiComm == Mpi_Comm_Null ) return
    ! stop "gather_vector_MPI error: MpiComm == MPI_UNDEFINED"
    !
    MpiSize   = get_size_MPI(MpiComm)
    MpiMaster = get_master_MPI(MpiComm)
    !
    Nloc = size(Vloc)
    N    = 0
    call AllReduce_MPI(MpiComm,Nloc,N)
    if(MpiMaster.AND.N /= size(V)) stop "allgather_vector_MPI error: size(V) != Mpi_Allreduce(Nloc)"
    !
    allocate(Counts(0:MpiSize-1)) ; Counts=0
    allocate(Offset(0:MpiSize-1)) ; Offset=0
    !
    !Get Counts;
    call MPI_AllGather(Nloc/Dimph,1,MPI_INTEGER,Counts,1,MPI_INTEGER,MpiComm,MpiIerr)
    !
    !Get Offset:
    Offset(0)=0
    do i=1,MpiSize-1
       Offset(i) = Offset(i-1) + Counts(i-1)
    enddo
    !
    V = 0d0
    do iph=1,Dimph
       v_start = 1 + (iph-1)*(N/Dimph)
       v_end = iph*(N/Dimph)
       vloc_start = 1 + (iph-1)*(Nloc/Dimph)
       vloc_end = iph*(Nloc/Dimph)
       call MPI_AllGatherv(Vloc(vloc_start:vloc_end),Nloc/DimPh,MPI_DOUBLE_PRECISION,&
            V(v_start:v_end),Counts,Offset,MPI_DOUBLE_PRECISION,MpiComm,MpiIerr)
    enddo
    !
    return
  end subroutine d_allgather_vector_MPI

  subroutine c_allgather_vector_MPI(MpiComm,vloc,v)
    integer                          :: MpiComm
    complex(8),dimension(:)          :: vloc !size[Nloc]
    complex(8),dimension(:)          :: v    !size[N]
    integer                          :: i,irank,Nloc,N
    integer,dimension(:),allocatable :: Counts,Offset
    integer                          :: MpiSize,MpiIerr
    logical                          :: MpiMaster
    !
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A)")"DEBUG c_allgather_basis_MPI: allgather v"
#endif
    !
    if( MpiComm == MPI_UNDEFINED ) stop "gather_vector_MPI error: MpiComm == MPI_UNDEFINED"
    !
    MpiSize   = get_size_MPI(MpiComm)
    MpiMaster = get_master_MPI(MpiComm)
    !
    Nloc = size(Vloc)
    N = 0
    call AllReduce_MPI(MpiComm,Nloc,N)
    if(MpiMaster.AND.N /= size(V)) stop "gather_vector_MPI error: size(V) != Mpi_Allreduce(Nloc)"
    !
    allocate(Counts(0:MpiSize-1)) ; Counts=0
    allocate(Offset(0:MpiSize-1)) ; Offset=0
    !
    !Get Counts;
    call MPI_AllGather(Nloc,1,MPI_INTEGER,Counts,1,MPI_INTEGER,MpiComm,MpiIerr)
    !
    !Get Offset:
    Offset(0)=0
    do i=1,MpiSize-1
       Offset(i) = Offset(i-1) + Counts(i-1)
    enddo
    !
    call MPI_AllGatherv(Vloc,Nloc,MPI_DOUBLE_COMPLEX,V,Counts,Offset,MPI_DOUBLE_COMPLEX,MpiComm,MpiIerr)
    !
    return
  end subroutine c_allgather_vector_MPI
#endif



  !##################################################################
  !                   RESHAPE ROUTINES
  !##################################################################

  !> Get stride position in the one-particle many-body space 
  function index_stride_so(ispin,iorb) result(indx)
    integer :: iorb
    integer :: ispin
    integer :: indx
    indx = iorb  + (ispin-1)*Norb
  end function index_stride_so



  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  ! reshape a matrix from the [Nlso][Nlso] shape
  ! from/to the [Nlat][Nspin][Nspin][Norb][Norb] shape.
  ! _nlso2nnn : from [Nlso][Nlso] to [Nlat][Nspin][Nspin][Norb][Norb]  !
  ! _nso2nn   : from [Nso][Nso]   to [Nspin][Nspin][Norb][Norb]
  !+-----------------------------------------------------------------------------+!
  function d_nlso2nnn(Hlso,Nlat,Nspin,Norb) result(Hnnn)
    integer                                            :: Nlat,Nspin,Norb
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                            :: iorb,ispin,ilat,is
    integer                                            :: jorb,jspin,js
    Hnnn=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hnnn(ilat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nlso2nnn
  function c_nlso2nnn(Hlso,Nlat,Nspin,Norb) result(Hnnn)
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,js
    Hnnn=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hnnn(ilat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_nlso2nnn

  function d_nso2nn(Hso,Nspin,Norb) result(Hnn)
    integer                                  :: Nspin,Norb
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                  :: iorb,ispin,is
    integer                                  :: jorb,jspin,js
    Hnn=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hnn(ispin,jspin,iorb,jorb) = Hso(is,js)
             enddo
          enddo
       enddo
    enddo
  end function d_nso2nn
  function c_nso2nn(Hso,Nspin,Norb) result(Hnn)
    integer                                     :: Nspin,Norb
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                     :: iorb,ispin,is
    integer                                     :: jorb,jspin,js
    Hnn=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hnn(ispin,jspin,iorb,jorb) = Hso(is,js)
             enddo
          enddo
       enddo
    enddo
  end function c_nso2nn




  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  ! reshape a matrix from the [Nlat][Nspin][Nspin][Norb][Norb] shape
  ! from/to the [Nlso][Nlso] shape.
  ! _nnn2nlso : from [Nlat][Nspin][Nspin][Norb][Norb] to [Nlso][Nlso]
  ! _nn2nso   : from [Nspin][Nspin][Norb][Norb]       to [Nso][Nso]
  !+-----------------------------------------------------------------------------+!
  function d_nnn2nlso(Hnnn,Nlat,Nspin,Norb) result(Hlso)
    integer                                            :: Nlat,Nspin,Norb
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                            :: iorb,ispin,ilat,is
    integer                                            :: jorb,jspin,js
    Hlso=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hlso(is,js) = Hnnn(ilat,ispin,jspin,iorb,jorb)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nnn2nlso

  function c_nnn2nlso(Hnnn,Nlat,Nspin,Norb) result(Hlso)
    integer                                            :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,js
    Hlso=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hlso(is,js) = Hnnn(ilat,ispin,jspin,iorb,jorb)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_nnn2nlso

  function d_nn2nso(Hnn,Nspin,Norb) result(Hso)
    integer                                  :: Nspin,Norb
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                  :: iorb,ispin,is
    integer                                  :: jorb,jspin,js
    Hso=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hso(is,js) = Hnn(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function d_nn2nso

  function c_nn2nso(Hnn,Nspin,Norb) result(Hso)
    integer                                  :: Nspin,Norb
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                     :: iorb,ispin,is
    integer                                     :: jorb,jspin,js
    Hso=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hso(is,js) = Hnn(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function c_nn2nso




  function so2os_reshape(fg,Nspin,Norb) result(g)
    integer                                     :: Nspin,Norb
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: fg
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: g
    integer                                     :: iorb,jorb,ispin,jspin
    integer                                     :: io1,jo1,io2,jo2
    g = zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                !O-index
                io1 = iorb + (ispin-1)*Norb
                jo1 = jorb + (jspin-1)*Norb
                !I-index
                io2 = ispin + (iorb-1)*Nspin
                jo2 = jspin + (jorb-1)*Nspin
                !switch
                g(io1,jo1)  = fg(io2,jo2)
                !
             enddo
          enddo
       enddo
    enddo
  end function so2os_reshape

  function os2so_reshape(fg,Nspin,Norb) result(g)
    integer                                     :: Nspin,Norb
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: fg
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: g
    integer                                     :: iorb,jorb,ispin,jspin
    integer                                     :: io1,jo1,io2,jo2
    g = zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                !O-index
                io1 = ispin + (iorb-1)*Nspin
                jo1 = jspin + (jorb-1)*Nspin
                !I-index
                io2 = iorb + (ispin-1)*Norb
                jo2 = jorb + (jspin-1)*Norb
                !switch
                g(io1,jo1)  = fg(io2,jo2)
                !
             enddo
          enddo
       enddo
    enddo
  end function os2so_reshape













  !+------------------------------------------------------------------+
  !PURPOSE  : Allocate arrays and setup frequencies and times
  !+------------------------------------------------------------------+
  subroutine allocate_grids
    integer :: i
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(vm))allocate(vm(0:Lmats))          !bosonic frequencies
    if(.not.allocated(wr))allocate(wr(Lreal))
    if(.not.allocated(vr))allocate(vr(Lreal))
    if(.not.allocated(tau))allocate(tau(0:Ltau))
    wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    do i=0,Lmats
       vm(i) = pi/beta*2*i
    enddo
    wr     = linspace(wini,wfin,Lreal)
    vr     = linspace(0d0,wfin,Lreal)
    tau(0:)= linspace(0d0,beta,Ltau+1)
  end subroutine allocate_grids


  subroutine deallocate_grids
    if(allocated(wm))deallocate(wm)
    if(allocated(vm))deallocate(vm)
    if(allocated(tau))deallocate(tau)
    if(allocated(wr))deallocate(wr)
    if(allocated(vr))deallocate(vr)
  end subroutine deallocate_grids










  !##################################################################
  !##################################################################
  ! ROUTINES TO SEARCH CHEMICAL POTENTIAL UP TO SOME ACCURACY
  ! can be used to fix any other *var so that  *ntmp == nread
  !##################################################################
  !##################################################################

  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine ed_search_variable(var,ntmp,converged)
    real(8),intent(inout) :: var
    real(8),intent(in)    :: ntmp
    logical,intent(inout) :: converged
    logical               :: bool
    real(8),save          :: chich
    real(8),save          :: nold
    real(8),save          :: var_new
    real(8),save          :: var_old
    real(8)               :: var_sign
    real(8)               :: delta_n,delta_v,chi_shift
    !
    real(8)               :: ndiff
    integer,save          :: count=0,totcount=0,i
    integer               :: unit
    logical :: master
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG ed_search_variable: adjust var"
#endif
    !
    if(nread==0d0)return
    master=.true.
#ifdef _MPI    
    if(check_MPI())master  = get_master_MPI()
#endif
    !
    if(master)then
       !check actual value of the density *ntmp* with respect to goal value *nread*
       count=count+1
       totcount=totcount+1
       !  
       if(count==1)then
          chich = ndelta        !~0.2
          inquire(file="var_compressibility.restart",EXIST=bool)
          if(bool)then
             write(LOGfile,"(A)")"Reading compressibility from file"
             open(free_unit(unit),file="var_compressibility.restart")
             read(unit,*)chich
             close(unit)
          endif
          var_old = var
       endif
       !
       ndiff=ntmp-nread
       !
       open(free_unit(unit),file="var_compressibility.used")
       write(unit,*)chich
       close(unit)
       !
       ! if(abs(ndiff)>nerr)then
       !Get 'charge compressibility"
       delta_n = ntmp-nold
       delta_v = var-var_old
       if(count>1)chich = delta_v/(delta_n+1d-10) !1d-4*nerr)  !((ntmp-nold)/(var-var_old))**-1
       !
       !Add here controls on chich: not to be too small....
       if(chich>10d0)chich=2d0*chich/abs(chich) !do nothing?
       !
       chi_shift = ndiff*chich
       !
       !update chemical potential
       var_new = var - chi_shift
       !
       !
       !re-define variables:
       nold    = ntmp
       var_old = var
       var     = var_new
       !
       !Print information
       write(LOGfile,"(A11,F16.9,A,F15.9)")  "n      = ",ntmp,"| instead of",nread
       write(LOGfile,"(A11,ES16.9,A,ES16.9)")"n_diff = ",ndiff,"/",nerr
       write(LOGfile,"(A11,ES16.9,A,ES16.9)")"dv     = ",delta_v
       write(LOGfile,"(A11,ES16.9,A,ES16.9)")"dn     = ",delta_n
       write(LOGfile,"(A11,F16.9,A,F15.9)")  "dv/dn  = ",chich
       var_sign = (var-var_old)/abs(var-var_old)
       if(var_sign>0d0)then
          write(LOGfile,"(A11,ES16.9,A4)")"shift    = ",chi_shift," ==>"
       else
          write(LOGfile,"(A11,ES16.9,A4)")"shift    = ",chi_shift," <=="
       end if
       write(LOGfile,"(A11,F16.9)")"var     = ",var
       !
       ! else
       !    count=0
       ! endif
       !Save info about search variable iteration:
       open(free_unit(unit),file="search_variable_iteration_info"//reg(ed_file_suffix)//".ed",position="append")
       ! if(count==1)write(unit,*)"#var,ntmp,ndiff"
       write(unit,*)totcount,var,ntmp,ndiff
       close(unit)
       !
       !If density is not converged set convergence to .false.
       if(abs(ndiff)>nerr)converged=.false.
       !
       write(LOGfile,"(A18,I5)")"Search var count= ",count
       write(LOGfile,"(A19,L2)")"Converged       = ",converged
       print*,""
       !
       open(free_unit(unit),file="var_compressibility.restart")
       write(unit,*)chich
       close(unit)
       !
    endif
#ifdef _MPI
    if(check_MPI())then
       call Bcast_MPI(MPI_COMM_WORLD,converged)
       call Bcast_MPI(MPI_COMM_WORLD,var)
    endif
#endif
  end subroutine ed_search_variable




  subroutine search_chemical_potential(var,ntmp,converged)
    real(8),intent(inout) :: var
    real(8),intent(in)    :: ntmp
    logical,intent(inout) :: converged
    logical               :: bool
    real(8)               :: ndiff
    integer,save          :: count=0,totcount=0,i
    integer,save          :: nindex=0
    integer               :: nindex_old(3)
    real(8)               :: ndelta_old,nratio,kcompr=0d0
    integer,save          :: nth_magnitude=-2,nth_magnitude_old=-2
    real(8),save          :: nth=1.d-2,var_old,ntmp_old
    logical,save          :: ireduce=.true.
    integer               :: unit
    logical :: master
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG search_chemical_potential: adjust mu"
#endif
    !
    if(nread==0d0)return
    master=.true.
#ifdef _MPI    
    if(check_MPI())master  = get_master_MPI()
#endif
    !
    if(master)then
       !
       if(count==0)then
          inquire(file="var.restart",EXIST=bool)
          if(bool)then
             open(free_unit(unit),file="var.restart")
             read(unit,*)var,ndelta
             ndelta=abs(ndelta)*ncoeff
             close(unit)
          endif
       endif
       !
       ndiff=ntmp-nread
       nratio = 0.5d0;!nratio = 1.d0/(6.d0/11.d0*pi)
       !
       !check actual value of the density *ntmp* with respect to goal value *nread*
       count=count+1
       totcount=totcount+1
       if(count>2)then
          do i=1,2
             nindex_old(i+1)=nindex_old(i)
          enddo
       endif
       nindex_old(1)=nindex
       !
       if(ndiff >= nth)then
          nindex=-1
       elseif(ndiff <= -nth)then
          nindex=1
       else
          nindex=0
       endif
       !
       !
       ndelta_old=ndelta
       bool=nindex/=0.AND.( (nindex+nindex_old(1)==0).OR.(nindex+sum(nindex_old(:))==0) )
       !if(nindex_old(1)+nindex==0.AND.nindex/=0)then !avoid loop forth and back
       if(bool)then
          ndelta=ndelta_old*nratio !decreasing the step
       else
          ndelta=ndelta_old
       endif
       !
       if(abs(ndelta_old)<1.d-9)then
          ndelta_old=0.d0
          nindex=0
       endif
       !
       !update chemical potential
       var=var+dble(nindex)*ndelta
       !
       ! if(count>0)kcompr = (ntmp - ntmp_old)/(var - var_old)       
       !
       !Print information
       write(LOGfile,"(A,f16.9,A,f15.9)")"n    = ",ntmp," /",nread
       if(nindex>0)then
          write(LOGfile,"(A,I0,A1,es16.9,A)")"shift= ",nindex,"*",ndelta," ==>"
       elseif(nindex<0)then
          write(LOGfile,"(A,I0,A1,es16.9,A)")"shift= ",nindex,"*",ndelta," <=="
       else
          write(LOGfile,"(A,I0,A1,es16.9,A)")"shift= ",nindex,"*",ndelta," ==="
       endif
       write(LOGfile,"(A,ES16.9,A,ES16.9)")"dn   = ",ndiff,"/",nth
       ! write(LOGfile,"(A,10F16.9)")"k    = ",kcompr,1d0/kcompr,ntmp,ntmp_old,ntmp-ntmp_old,var,var_old,var-var_old
       write(LOGfile,"(A,f15.9)")"var  = ",var
       !
       open(free_unit(unit),file="search_mu_iteration"//reg(ed_file_suffix)//".ed",position="append")
       write(unit,*)var,ntmp,ndiff
       close(unit)
       !
       !check convergence within actual threshold
       !if reduce is activetd
       !if density is in the actual threshold
       !if DMFT is converged
       !if threshold is larger than nerror (i.e. this is not last loop)
       bool=ireduce.AND.(abs(ndiff)<nth).AND.converged.AND.(nth>nerr)
       if(bool)then
          nth_magnitude_old=nth_magnitude        !save old threshold magnitude
          nth_magnitude=nth_magnitude_old-1      !decrease threshold magnitude || floor(log10(abs(ntmp-nread)))
          nth=max(nerr,10.d0**(nth_magnitude))   !set the new threshold 
          count=0                                !reset the counter
          converged=.false.                      !reset convergence
          ndelta=ndelta_old*nratio                  !reduce the delta step
          !
       endif
       !
       !if density is not converged set convergence to .false.
       if(abs(ntmp-nread)>nth)converged=.false.
       !
       !check convergence for this threshold
       !!---if smallest threshold-- NO MORE
       !if reduce is active (you reduced the treshold at least once)
       !if # iterations > max number
       !if not yet converged
       !set threshold back to the previous larger one.
       !bool=(nth==nerr).AND.ireduce.AND.(count>niter).AND.(.not.converged)
       bool=ireduce.AND.(count>niter).AND.(.not.converged)
       if(bool)then
          ireduce=.false.
          nth=10.d0**(nth_magnitude_old)
       endif
       !
       write(LOGfile,"(A,I5)")"count= ",count
       write(LOGfile,"(A,L2)")"Converged=",converged
       print*,""
       !
       open(free_unit(unit),file="var.restart")
       write(unit,*)var,ndelta
       close(unit)
       !
       ntmp_old = ntmp
       var_old  = var
       !
    endif
#ifdef _MPI
    if(check_MPI())then
       call Bcast_MPI(MPI_COMM_WORLD,converged)
       call Bcast_MPI(MPI_COMM_WORLD,var)
    endif
#endif
  end subroutine search_chemical_potential










  subroutine SOC_jz_symmetrize(funct,mask)
    !passed
    complex(8),allocatable,intent(inout)         ::  funct(:,:,:,:,:)
    logical(8),allocatable,intent(in)            ::  mask(:,:,:,:,:)
    complex(8),allocatable                       ::  funct_in(:,:,:),funct_out(:,:,:)
    complex(8),allocatable                       ::  a_funct(:),b_funct(:)
    integer                                      ::  ispin,io
    integer                                      ::  ifreq,Lfreq
    logical(8)                                   ::  boolmask
    complex(8),allocatable                       ::  U(:,:),Udag(:,:)
    if(size(funct,dim=1)/=Nspin)stop "wrong size 1 in SOC symmetrize input f"
    if(size(funct,dim=2)/=Nspin)stop "wrong size 2 in SOC symmetrize input f"
    if(size(funct,dim=3)/=Norb) stop "wrong size 3 in SOC symmetrize input f"
    if(size(funct,dim=4)/=Norb) stop "wrong size 4 in SOC symmetrize input f"
    Lfreq=size(funct,dim=5)
    allocate(funct_in(Nspin*Norb,Nspin*Norb,Lfreq)); funct_in=zero
    allocate(funct_out(Nspin*Norb,Nspin*Norb,Lfreq));funct_out=zero
    allocate(U(Nspin*Norb,Nspin*Norb));U=zero
    allocate(Udag(Nspin*Norb,Nspin*Norb));Udag=zero
    allocate(a_funct(Lfreq));a_funct=zero
    allocate(b_funct(Lfreq));b_funct=zero
    !
    !function intake
    do ifreq=1,Lfreq
       funct_in(:,:,ifreq)=nn2so_reshape(funct(:,:,:,:,ifreq),Nspin,Norb)
    enddo
    !
    !function diagonalization
    if(Jz_basis)then
       U=matmul(transpose(conjg(orbital_Lz_rotation_NorbNspin())),atomic_SOC_rotation())
       Udag=transpose(conjg(U))
    else
       U=atomic_SOC_rotation()
       Udag=transpose(conjg(U))
    endif
    !
    do ifreq=1,Lfreq
       funct_out(:,:,ifreq)=matmul(Udag,matmul(funct_in(:,:,ifreq),U))
    enddo
    !
    !function symmetrization in the rotated basis
    do io=1,2
       a_funct(:)=a_funct(:)+funct_out(io,io,:)
    enddo
    a_funct = a_funct/2.d0
    do io=3,6
       b_funct(:)=b_funct(:)+funct_out(io,io,:)
    enddo
    b_funct = b_funct/4.d0
    !
    boolmask = .false.
    if(Jz_basis)then
       boolmask = (.not.mask(1,2,3,2,1)).and.(.not.mask(1,2,3,2,2))
    else
       boolmask = (.not.mask(1,2,3,1,1)).and.(.not.mask(1,2,3,1,2))
    endif
    if(boolmask)then
       a_funct = ( a_funct + b_funct ) / 2.d0
       b_funct = a_funct
    endif
    !
    funct_out=zero
    do io=1,2
       funct_out(io,io,:)=a_funct(:)
    enddo
    do io=3,6
       funct_out(io,io,:)=b_funct(:)
    enddo
    !
    !function rotation in the non-diagonal basis
    funct_in=zero
    do ifreq=1,Lfreq
       funct_in(:,:,ifreq)=matmul(U,matmul(funct_out(:,:,ifreq),Udag))
    enddo
    !
    !founction out
    funct=zero
    do ifreq=1,Lfreq
       funct(:,:,:,:,ifreq)=so2nn_reshape(funct_in(:,:,ifreq),Nspin,Norb)
    enddo
  end subroutine SOC_jz_symmetrize




  !+-------------------------------------------------------------------+
  !PURPOSE  : Atomic SOC and j vector components
  !+-------------------------------------------------------------------+
  function atomic_SOC() result (LS)
    complex(8),dimension(Nspin*Norb,Nspin*Norb)  :: LS,LS_
    integer                                      :: i,j
    LS_=zero;LS=zero
    LS_(1:2,3:4) = +Xi * pauli_z / 2.
    LS_(1:2,5:6) = -Xi * pauli_y / 2.
    LS_(3:4,5:6) = +Xi * pauli_x / 2.
    !hermiticity
    do i=1,Nspin*Norb
       do j=1,Nspin*Norb
          LS_(j,i)=conjg(LS_(i,j))
       enddo
    enddo
    LS=so2os_reshape(LS_,Nspin,Norb)
  end function atomic_SOC

  function atomic_SOC_rotation() result (LS_rot)
    complex(8),dimension(Nspin*Norb,Nspin*Norb)  :: LS_rot,LS_rot_
    integer                                      :: i,j
    LS_rot_=zero;LS_rot=zero
    !
    ! {a,Sz}-->{J}
    !
    ![Norb*Norb]*Nspin notation
    !J=1/2 jz=-1/2
    LS_rot_(1,1)=+1.d0
    LS_rot_(2,1)=-Xi
    LS_rot_(6,1)=-1.d0
    LS_rot_(:,1)=LS_rot_(:,1)/sqrt(3.)
    !J=1/2 jz=+1/2
    LS_rot_(4,2)=+1.d0
    LS_rot_(5,2)=+Xi
    LS_rot_(3,2)=+1.d0
    LS_rot_(:,2)=LS_rot_(:,2)/sqrt(3.)
    !J=3/2 jz=-3/2
    LS_rot_(4,3)=+1.d0
    LS_rot_(5,3)=-Xi
    LS_rot_(:,3)=LS_rot_(:,3)/sqrt(2.)
    !J=3/2 jz=+3/2
    LS_rot_(1,4)=-1.d0
    LS_rot_(2,4)=-Xi
    LS_rot_(:,4)=LS_rot_(:,4)/sqrt(2.)
    !J=3/2 jz=-1/2
    LS_rot_(1,5)=+1.d0
    LS_rot_(2,5)=-Xi
    LS_rot_(6,5)=+2.d0
    LS_rot_(:,5)=LS_rot_(:,5)/sqrt(6.)
    !J=3/2 jz=+1/2
    LS_rot_(4,6)=-1.d0
    LS_rot_(5,6)=-Xi
    LS_rot_(3,6)=+2.d0
    LS_rot_(:,6)=LS_rot_(:,6)/sqrt(6.)
    !
    LS_rot=LS_rot_
    !
  end function atomic_SOC_rotation

  function orbital_Lz_rotation_Norb() result (U_rot)
    complex(8),dimension(Norb,Norb)              :: U_rot,U_rot_
    integer                                      :: i,j
    U_rot=zero;U_rot_=zero
    !
    ! {a}-->{Lz}
    !
    ![Norb*Norb] notation
    U_rot_(1,1)=-Xi/sqrt(2.)
    U_rot_(2,2)=+1.d0/sqrt(2.)
    U_rot_(3,3)=+Xi
    U_rot_(1,2)=-Xi/sqrt(2.)
    U_rot_(2,1)=-1.d0/sqrt(2.)
    !
    U_rot=U_rot_
    !
  end function orbital_Lz_rotation_Norb

  function orbital_Lz_rotation_NorbNspin() result (U_rot)
    complex(8),dimension(Norb,Norb)              :: U_rot_
    complex(8),dimension(Norb*Nspin,Norb*Nspin)  :: U_rot
    integer                                      :: i,j
    U_rot=zero;U_rot_=zero
    !
    ! {a,Sz}-->{Lz,Sz}
    !
    ![Norb*Norb]*Nspin notation
    U_rot_(1,1)=-Xi/sqrt(2.)
    U_rot_(2,2)=+1.d0/sqrt(2.)
    U_rot_(3,3)=+Xi
    U_rot_(1,2)=-Xi/sqrt(2.)
    U_rot_(2,1)=-1.d0/sqrt(2.)
    !
    U_rot(1:Norb,1:Norb)=U_rot_
    U_rot(1+Norb:2*Norb,1+Norb:2*Norb)=U_rot_
    !
  end function orbital_Lz_rotation_NorbNspin

  function atomic_j(component) result (ja)
    complex(8),dimension(Nspin*Norb,Nspin*Norb)  :: ja,ja_
    character(len=1)                             :: component
    integer                                      :: i,j
    ja_=zero;ja=zero
    if    (component=="x")then
       ja_(1:2,1:2) = pauli_x / 2.
       ja_(3:4,3:4) = pauli_x / 2.
       ja_(5:6,5:6) = pauli_x / 2.
       ja_(3:4,5:6) = -Xi * eye(2)
    elseif(component=="y")then
       ja_(1:2,1:2) = pauli_y / 2.
       ja_(3:4,3:4) = pauli_y / 2.
       ja_(5:6,5:6) = pauli_y / 2.
       ja_(1:2,5:6) = +Xi * eye(2)
    elseif(component=="z")then
       ja_(1:2,1:2) = pauli_z / 2.
       ja_(3:4,3:4) = pauli_z / 2.
       ja_(5:6,5:6) = pauli_z / 2.
       ja_(1:2,3:4) = -Xi * eye(2)
    endif
    !hermiticity
    do i=1,Nspin*Norb
       do j=1,Nspin*Norb
          ja_(j,i)=conjg(ja_(i,j))
       enddo
    enddo
    ja=so2os_reshape(ja_,Nspin,Norb)
  end function atomic_j








  !+------------------------------------------------------------------+
  !PURPOSE  : Print Hloc
  !+------------------------------------------------------------------+
  subroutine print_Hloc_nn_c(hloc,file)
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: hloc
    character(len=*),optional                   :: file
    integer                                     :: iorb,jorb,ispin,jspin,Nso,unit
    character(len=32)                           :: fmt
    !
    unit=LOGfile
    !
    if(present(file))then
       open(free_unit(unit),file=reg(file))
       write(LOGfile,"(A)")"print_Hloc on file :"//reg(file)
    endif
    !
    do ispin=1,Nspin
       do iorb=1,Norb
          write(unit,"(100(A1,F8.4,A1,F8.4,A1,2x))")&
               (&
               (&
               '(',dreal(Hloc(ispin,jspin,iorb,jorb)),',',dimag(Hloc(ispin,jspin,iorb,jorb)),')',&
               jorb =1,Norb),&
               jspin=1,Nspin)
       enddo
    enddo
    write(unit,*)""
    !
    if(present(file))close(unit)
  end subroutine print_Hloc_nn_c

  subroutine print_Hloc_so_c(hloc,file)
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: hloc
    character(len=*),optional                   :: file
    integer                                     :: is,js,Nso,unit
    character(len=32)                           :: fmt
    !
    unit=LOGfile
    !
    if(present(file))then
       open(free_unit(unit),file=reg(file))
       write(LOGfile,"(A)")"print_Hloc on file :"//reg(file)
    endif
    !
    Nso = Nspin*Norb
    do is=1,Nso
       write(unit,"(20(A1,F8.4,A1,F8.4,A1,2x))")&
            ('(',dreal(Hloc(is,js)),',',dimag(Hloc(is,js)),')',js =1,Nso)
    enddo
    write(unit,*)""
    !
    if(present(file))close(unit)
  end subroutine print_Hloc_so_c
  !





  !+------------------------------------------------------------------+
  !PURPOSE  : print a state vector |{up}>|{dw}>
  !+------------------------------------------------------------------+
  subroutine print_state_vector_ivec(ivec,unit)
    integer,intent(in) :: ivec(:)
    integer,optional   :: unit
    integer            :: unit_
    integer            :: i,j,Ntot
    character(len=2)   :: fbt
    character(len=16)  :: fmt
    unit_=6;if(present(unit))unit_=unit
    Ntot = size(ivec)
    write(fbt,'(I2.2)')Ntot
    fmt="(B"//adjustl(trim(fbt))//"."//adjustl(trim(fbt))//")"
    i= bjoin(ivec,Ntot)
    write(unit_,"(I9,1x,A1)",advance="no")i,"|"
    write(unit_,"(10I1)",advance="no")(ivec(j),j=1,Ntot)
    write(unit_,"(A4)",advance="no")"> - "
    write(unit_,fmt,advance="yes")i
  end subroutine print_state_vector_ivec
  !
  subroutine  print_state_vector_ivec_ud(ivec,jvec,unit)
    integer,intent(in) :: ivec(:),jvec(size(ivec))
    integer,optional   :: unit
    integer            :: unit_
    integer            :: i,j,iup,idw,Ntot
    character(len=2)   :: fbt
    character(len=20)  :: fmt
    unit_=6;if(present(unit))unit_=unit
    Ntot = size(ivec)
    write(fbt,'(I2.2)')Ntot
    fmt="(B"//adjustl(trim(fbt))//"."//adjustl(trim(fbt))//",1x,B"//adjustl(trim(fbt))//"."//adjustl(trim(fbt))//")"
    iup = bjoin(ivec,Ntot)
    idw = bjoin(jvec,Ntot)
    i = bjoin([ivec,jvec],2*Ntot)
    write(unit_,"(I9,1x,I4,1x,A1)",advance="no")i,iup,"|"
    write(unit_,"(10I1)",advance="no")(ivec(j),j=1,Ntot)
    write(unit_,"(A1,I4,A2)",advance="no")">",idw," |"
    write(unit_,"(10I1)",advance="no")(jvec(j),j=1,Ntot)
    write(unit_,"(A4)",advance="no")"> - "
    write(unit_,fmt,advance="yes")ibits(i,0,Ntot),ibits(i,Ntot,2*Ntot)
  end subroutine print_state_vector_ivec_ud
  !
  subroutine print_state_vector_int(i,Ntot,unit)
    integer,intent(in) :: i
    integer,intent(in) :: Ntot
    integer,optional   :: unit
    integer            :: unit_
    integer            :: j
    integer            :: ivec(Ntot)
    character(len=2)   :: fbt
    character(len=16)  :: fmt
    unit_=6;if(present(unit))unit_=unit
    write(fbt,'(I2.2)')Ntot
    fmt="(B"//adjustl(trim(fbt))//"."//adjustl(trim(fbt))//")"
    ivec = bdecomp(i,Ntot)
    write(unit_,"(I9,1x,A1)",advance="no")i,"|"
    write(unit_,"(10I1)",advance="no")(ivec(j),j=1,Ntot)
    write(unit_,"(A4)",advance="no")"> - "
    write(unit_,fmt,advance="yes")i
  end subroutine print_state_vector_int


END MODULE ED_AUX_FUNX
