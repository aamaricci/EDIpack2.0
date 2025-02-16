MODULE ED_VARS_GLOBAL
  !
  !Contains all variables, arrays and derived types instances shared throughout the code.
  !Specifically, it contains definitions of the :f:var:`effective_bath`, the :f:var:`gfmatrix` and the :f:var:`sector` data structures. 
  !
  USE SF_CONSTANTS
  USE SF_IOTOOLS, only:free_unit,reg,str
  USE ED_SPARSE_MATRIX
  USE ED_SPARSE_MAP
  USE ED_GFMATRIX
  USE ED_INPUT_VARS, only:ed_verbose,logfile  
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none


  type H_operator
     ! The matrix storing in the basis [:f:var:`nspin` , :f:var:`nspin` , :f:var:`norb` , :f:var:`norb` ] each element of the Matrix basis decomposing the replica/general bath Hamiltonian :math:`H_p=\sum_{i=1}^{N_{basis}} \lambda_i(p) O_i`, where :math:`N_{basis}` is the dimension of the user defined basis.  
     complex(8),dimension(:,:,:,:),allocatable               :: O !Replica/General hamiltonian
  end type H_operator

  type effective_bath_component
     ! Effective bath component for the replica/general bath. Each istance of this type defines the parameters :math:`\vec{\lambda}` and the amplitudes :math:`\vec{V}`. The first is used to decompose the Hamiltonian of each element of the bath :math:`H_p=\sum_{i=1}^{N_{basis}} \lambda_i(p) O_i`, the latter describes the hopping from/to the impurity.
     real(8)                                                 :: v
     real(8),dimension(:),allocatable                        :: vg
     real(8),dimension(:),allocatable                        :: lambda ![:f:var:`nsym`]
  end type effective_bath_component

  type effective_bath
     ! This structure describes the (effective) discretized bath used in the contruct the Hamiltonian of the quantum impurity system. Each element of this structure is allocated and used according the value of :f:var:`ed_mode` = :code:`normal,superc,nonsu2` and :f:var:`bath_type` = :code:`normal,hybrid,replica,general`.  
     real(8),dimension(:,:,:),allocatable                    :: e !local energies [ :f:var:`nspin` ][ :f:var:`norb` ][ :f:var:`bath` ]/[ :f:var:`nspin` ][ :code:`1` ][ :f:var:`nspin` ]
     real(8),dimension(:,:,:),allocatable                    :: v !spin-keep hyb. [ :f:var:`nspin` ][ :f:var:`norb` ][ :f:var:`nbath` ]
     real(8),dimension(:,:,:),allocatable                    :: d !SC amplitues   [ :f:var:`nspin` ][ :f:var:`norb` ][ :f:var:`nbath` ]/[ :f:var:`norb` ][ :code:`1` ][ :f:var:`norb` ]
     real(8),dimension(:,:,:),allocatable                    :: u !spin-flip hyb. [ :f:var:`nspin` ][ :f:var:`norb` ][ :f:var:`nbath` ] for :f:var:`ed_mode` = :code:`nonsu2`
     integer                                                 :: Nbasis  !The replica/general Matrix basis dimension     
     type(effective_bath_component),dimension(:),allocatable :: item    ![ :f:var:`nbath` ] Replica/General bath components, V included
     logical                                                 :: status=.false.
  end type effective_bath





  !-------------------- CUSTOM OBSERVABLE STRUCTURE ----------------------!
  type observable
     complex(8),dimension(:,:,:),allocatable   :: sij ![Nlso][Nlso][Nk]
     character(len=32)                         :: o_name
     real(8)                                   :: o_value
  end type observable

  type custom_observables
     type(observable),dimension(:),allocatable :: item     ![:]
     complex(8),dimension(:,:,:),allocatable   :: Hk       ![Nlso][Nlso][Nk]
     integer                                   :: N_asked
     integer                                   :: N_filled
     logical                                   :: init=.false.
  end type custom_observables




  !---------------- SECTOR-TO-FOCK SPACE STRUCTURE -------------------!
  type sector_map
     integer,dimension(:),allocatable          :: map
     type(sparse_map)                          :: sp
     logical                                   :: status=.false.
  end type sector_map

  type sector
     integer                                   :: index       !
     type(sector_map),dimension(:),allocatable :: H
     integer,dimension(:),allocatable          :: DimUps
     integer,dimension(:),allocatable          :: DimDws
     integer                                   :: DimUp
     integer                                   :: DimDw
     integer                                   :: DimEl
     integer                                   :: DimPh
     integer                                   :: Dim
     integer,dimension(:),allocatable          :: Nups
     integer,dimension(:),allocatable          :: Ndws
     integer                                   :: Nup
     integer                                   :: Ndw
     integer                                   :: Sz
     integer                                   :: Ntot,twoJz
     integer                                   :: Nlanc
     logical                                   :: status=.false.
  end type sector




  !------------------ ABTRACT INTERFACES PROCEDURES ------------------!
  !SPARSE MATRIX-VECTOR PRODUCTS USED IN ED_MATVEC
  !dbleMat*dbleVec
  abstract interface
     subroutine dd_sparse_HxV(Nloc,v,Hv)
       integer                 :: Nloc
       real(8),dimension(Nloc) :: v
       real(8),dimension(Nloc) :: Hv
     end subroutine dd_sparse_HxV
  end interface


  !cmplxMat*cmplxVec
  abstract interface
     subroutine cc_sparse_HxV(Nloc,v,Hv)
       integer                    :: Nloc
       complex(8),dimension(Nloc) :: v
       complex(8),dimension(Nloc) :: Hv
     end subroutine cc_sparse_HxV
  end interface




  !-------------------------- ED  VARIABLES --------------------------!

  !SIZE OF THE PROBLEM
  !=========================================================
  integer,save                                       :: Ns       !Number of levels per spin
  integer                                            :: Nlevels
  integer,save                                       :: Nsectors !Number of sectors
  integer,save                                       :: Ns_orb
  integer,save                                       :: Ns_ud
  integer                                            :: Nhel
  integer,save                                       :: DimPh    !Number of phonon states

  !Some maps between sectors and full Hilbert space (pointers)
  !PRIVATE:
  !=========================================================
  integer,allocatable,dimension(:)                   :: getDim             ! [Nsectors]
  integer,allocatable,dimension(:,:,:)               :: getCsector         ! [1/Norb,2,NSectors]
  integer,allocatable,dimension(:,:,:)               :: getCDGsector       ! [1/Norb,2,NSectors]
  integer,allocatable,dimension(:,:)                 :: getBathStride
  integer,allocatable,dimension(:,:)                 :: impIndex
  logical,allocatable,dimension(:)                   :: twin_mask
  logical,allocatable,dimension(:)                   :: sectors_mask
  integer,allocatable,dimension(:,:)                 :: getSector
  integer,allocatable,dimension(:)                   :: getSz
  integer,allocatable,dimension(:)                   :: getN
  !
  integer,allocatable,dimension(:,:,:)               :: getCsector_Jz
  integer,allocatable,dimension(:,:,:)               :: getCDGsector_Jz
  integer,allocatable,dimension(:)                   :: gettwoJz
  integer,allocatable,dimension(:)                   :: getmaxtwoJz


  !Effective Bath used in the ED code (this is opaque to user)
  !=========================================================
  type(effective_bath)                               :: dmft_bath !instance of :f:var:`effective_bath` used to store the quantum impurity effective bath in the rest of the code 


  !Global Nambu factor for SC calculations (Nspin=1 but this index is 2 to
  !correctly allocate  Nambu arrays of dim 2*Norb) 
  !=========================================================
  integer                                            :: Nnambu=1
  !Replica/General bath basis set
  !=========================================================
  type(H_operator),dimension(:),allocatable          :: Hreplica_basis   ![Nsym]
  real(8),dimension(:,:),allocatable                 :: Hreplica_lambda  ![Nbath,Nsym]
  logical                                            :: Hreplica_status=.false.
  !
  type(H_operator),dimension(:),allocatable          :: Hgeneral_basis   ![Nsym]
  real(8),dimension(:,:),allocatable                 :: Hgeneral_lambda  ![Nbath,Nsym]
  logical                                            :: Hgeneral_status=.false.

  !local part of the Hamiltonian
  !=========================================================
  complex(8),dimension(:,:,:,:),allocatable          :: impHloc           !local hamiltonian



  !Variables for DIAGONALIZATION
  !PRIVATE
  !=========================================================  
  type(sparse_matrix_csr)                            :: spH0d !diagonal part
  type(sparse_matrix_csr)                            :: spH0nd !non-diagonal part
  type(sparse_matrix_csr),dimension(:),allocatable   :: spH0ups,spH0dws !reduced UP and DW parts
  type(sparse_matrix_csr)                            :: spH0_ph !Hamiltonian for phonons
  type(sparse_matrix_csr)                            :: spH0e_eph, spH0ph_eph !electron-phonon interaction
  type(sparse_matrix_csr)                            :: spH0
  procedure(dd_sparse_HxV),pointer                   :: spHtimesV_p=>null()
  procedure(cc_sparse_HxV),pointer                   :: spHtimesV_cc=>null()




  !Variables for DIAGONALIZATION
  !PRIVATE
  !=========================================================  
  integer,allocatable,dimension(:)                   :: neigen_sector
  logical                                            :: trim_state_list=.false.

  !Partition function
  !PRIVATE
  !=========================================================
  real(8)                                            :: zeta_function
  real(8)                                            :: gs_energy



  
  !Green's functions
  type(GFmatrix),allocatable,dimension(:,:,:,:)      :: impGmatrix
  type(GFmatrix)                                     :: impDmatrix
  !Spin Susceptibilities
  type(GFmatrix),allocatable,dimension(:,:)          :: spinChimatrix
  !Diagonal/Off-diagonal charge-charge Susceptibilities
  type(GFmatrix),allocatable,dimension(:,:)          :: densChimatrix
  !Pair-Pair Susceptibilities
  type(GFmatrix),allocatable,dimension(:,:)          :: pairChimatrix
  !Exciton Susceptibilities
  type(GFmatrix),allocatable,dimension(:,:,:)        :: exctChimatrix




  !Density and double occupancy
  !Local energies and generalized double occupancies
  !PRIVATE (now public but accessible thru routines)
  !=========================================================
  real(8),dimension(:),allocatable                   :: ed_dens
  real(8),dimension(:),allocatable                   :: ed_dens_up,ed_dens_dw
  real(8),dimension(:),allocatable                   :: ed_docc
  real(8),dimension(:,:),allocatable                 :: ed_phisc
  real(8),dimension(:,:),allocatable                 :: ed_mag
  real(8),dimension(:),allocatable                   :: ed_imp_info
  real(8)                                            :: ed_Ekin
  real(8)                                            :: ed_Epot
  real(8)                                            :: ed_Eint
  real(8)                                            :: ed_Ehartree
  real(8)                                            :: ed_Eknot
  real(8)                                            :: ed_Dust,ed_Dund,ed_Dse,ed_Dph


  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable                   :: wm,tau,wr,vm,vr


  !Impurity operators
  !PRIVATE (now public but accessible thru routine)
  !=========================================================
  complex(8),allocatable,dimension(:,:,:,:)          :: single_particle_density_matrix ![Nspin,Nspin,Norb,Norb]
  complex(8),allocatable,dimension(:,:)              :: impurity_density_matrix        ![2**2Norb,2**2Norb]
  integer,dimension(3)                               :: Lzdiag = [-1,+1,0]
  integer,dimension(2)                               :: Szdiag = [1,-1]
  real(8),dimension(:,:),allocatable                 :: spin_field ![Norb,3=x,y,z]


  !--------------- LATTICE WRAP VARIABLES -----------------!
  complex(8),dimension(:,:,:,:,:),allocatable        :: Hloc_ineq
  complex(8),dimension(:,:),allocatable,save         :: Dmats_ph_ineq,Dreal_ph_ineq
  complex(8),dimension(:,:,:,:,:),allocatable,save   :: single_particle_density_matrix_ineq
  complex(8),dimension(:,:,:),allocatable,save       :: impurity_density_matrix_ineq
  real(8),dimension(:,:),allocatable,save            :: dens_ineq 
  real(8),dimension(:,:),allocatable,save            :: docc_ineq
  real(8),dimension(:,:,:),allocatable,save          :: mag_ineq
  real(8),dimension(:,:,:),allocatable,save          :: phisc_ineq
  real(8),dimension(:,:),allocatable,save            :: dd_ineq,e_ineq
  integer,allocatable,dimension(:,:)                 :: neigen_sector_ineq
  integer,allocatable,dimension(:)                   :: neigen_total_ineq
  real(8),dimension(:,:,:),allocatable               :: Hreplica_lambda_ineq ![Nineq,Nbath,Nsym]
  real(8),dimension(:,:,:),allocatable               :: Hgeneral_lambda_ineq ![Nineq,Nbath,Nsym]



  !File suffixes for printing fine tuning.
  !=========================================================
  character(len=32)                                  :: ed_file_suffix=""       !suffix string attached to the output files.
  character(len=10)                                  :: ineq_site_suffix="_ineq"
  integer                                            :: site_indx_padding=4
  !logical                                           :: Jhflag              !spin-exchange and pair-hopping flag.
  logical                                            :: offdiag_gf_flag=.false.
  ! character(len=200)                               :: ed_input_file=""


  !This is the internal Mpi Communicator and variables.
  !=========================================================
#ifdef _MPI
  integer                                            :: MpiComm_Global=MPI_COMM_NULL
  integer                                            :: MpiComm=MPI_COMM_NULL
#endif
  integer                                            :: MpiGroup_Global=MPI_GROUP_NULL
  integer                                            :: MpiGroup=MPI_GROUP_NULL
  logical                                            :: MpiStatus=.false.
  logical                                            :: MpiMaster=.true.
  integer                                            :: MpiRank=0
  integer                                            :: MpiSize=1
  integer,allocatable,dimension(:)                   :: MpiMembers
  integer                                            :: mpiQup=0
  integer                                            :: mpiRup=0
  integer                                            :: mpiQdw=0
  integer                                            :: mpiRdw=0
  integer                                            :: mpiQ=0
  integer                                            :: mpiR=0
  integer                                            :: mpiIstart
  integer                                            :: mpiIend
  integer                                            :: mpiIshift
  logical                                            :: mpiAllThreads=.true.
  !

contains






  !=========================================================
  subroutine ed_set_MpiComm()
#ifdef _MPI
    integer :: ierr
    ! call MPI_Comm_dup(Comm,MpiComm_Global,ierr)
    ! call MPI_Comm_dup(Comm,MpiComm,ierr)
    MpiComm_Global = MPI_COMM_WORLD
    MpiComm        = MPI_COMM_WORLD
    call Mpi_Comm_group(MpiComm_Global,MpiGroup_Global,ierr)
    MpiStatus      = .true.
    MpiSize        = get_Size_MPI(MpiComm_Global)
    MpiRank        = get_Rank_MPI(MpiComm_Global)
    MpiMaster      = get_Master_MPI(MpiComm_Global)
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG ed_set_MpiComm: setting MPI comm"
#endif
#endif
  end subroutine ed_set_MpiComm

  subroutine ed_del_MpiComm()
#ifdef _MPI    
    MpiComm_Global = MPI_UNDEFINED
    MpiComm        = MPI_UNDEFINED
    MpiGroup_Global= MPI_GROUP_NULL
    MpiStatus      = .false.
    MpiSize        = 1
    MpiRank        = 0
    MpiMaster      = .true.
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG ed_del_MpiComm: deleting MPI comm"
#endif
#endif
  end subroutine ed_del_MpiComm
  !=========================================================




END MODULE ED_VARS_GLOBAL
