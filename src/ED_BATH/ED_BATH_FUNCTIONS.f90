MODULE ED_BATH_FUNCTIONS
  !A comprehensive set of procedures to evaluate the non-interacting impurity Green's functions :math:`\hat{G}^{\rm And}` and hybridizations :math:`\hat{F}^{\rm And}` in the complex frequency domain given the :f:var:`effective_bath` instance.
  !
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,txtfy,str
  USE SF_LINALG, only: eye,inv,diag,zeye,inv_her,kron,det
  USE SF_SPIN, only: pauli_sigma_z
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_BATH_AUX
  USE ED_BATH_DIM
  USE ED_BATH_USER
  USE ED_BATH_DMFT
  USE ED_BATH_REPLICA
  implicit none

  private


  !
  !\DELTA, THE HYBRIDIZATION FUNCTION
  interface delta_bath_function
     !
     ! Evaluates the normal hybridization function :math:`\Delta(x)`.
     !
     ! Output:
     !   * :f:var:`delta` : complex rank-5 array with dimension [ |Nspin| , |Nspin| , |Norb| , |Norb| , :code:`size(x)` ]
     !
     module procedure delta_bath_array
  end interface delta_bath_function

  interface Fdelta_bath_function
     !
     ! Evaluates the anomalous hybridization function :math:`\Theta(x)`.
     !
     ! Output:
     !   * :f:var:`fdelta` : complex rank-5 array with dimension [ |Nspin| , |Nspin| , |Norb| , |Norb| , :code:`size(x)` ]
     !
     module procedure Fdelta_bath_array
  end interface Fdelta_bath_function



  !NON-INTERACTING GREEN'S FUNCTION 
  interface g0and_bath_function
     !
     ! Evaluates the normal non-interacting Green's function :math:`G_0(x)`.
     !
     ! Output:
     !   * :f:var:`g0and` : complex rank-5 array with dimension [ |Nspin| , |Nspin| , |Norb| , |Norb| , :code:`size(x)` ]
     !
     module procedure g0and_bath_array
  end interface g0and_bath_function

  interface f0and_bath_function
     !
     ! Evaluates the anomalous non-interacting Green's function :math:`F_0(x)`.
     !
     ! Output:
     !   * :f:var:`f0and` : complex rank-5 array with dimension [ |Nspin| , |Nspin| , |Norb| , |Norb| , :code:`size(x)` ]
     !
     module procedure f0and_bath_array
  end interface f0and_bath_function



  !INVERSE NON-INTERACTING GREEN'S FUNCTION 
  interface invg0_bath_function
     !
     ! Evaluates the inverse of the normal non-interacting Green's function :math:`G^{-1}_0(x)`.
     !
     ! Output:
     !   * :f:var:`g0and` : complex rank-5 array with dimension [ |Nspin| , |Nspin| , |Norb| , |Norb| , :code:`size(x)` ]
     !
     module procedure invg0_bath_array
  end interface invg0_bath_function

  interface invf0_bath_function
     !
     ! Evaluates the inverse of the anomalous non-interacting Green's function :math:`F^{-1}_0(x)`.
     !
     ! Output:
     !   * :f:var:`f0and` : complex rank-5 array with dimension [ |Nspin| , |Nspin| , |Norb| , |Norb| , :code:`size(x)` ]
     !
     module procedure invf0_bath_array
  end interface invf0_bath_function




  public :: delta_bath_function
  public :: fdelta_bath_function
  public :: g0and_bath_function
  public :: f0and_bath_function
  public :: invg0_bath_function
  public :: invf0_bath_function






contains

  include "delta_functions/delta_normal.f90"
  include "delta_functions/delta_hybrid.f90"
  include "delta_functions/delta_replica.f90"
  include "delta_functions/delta_general.f90"
  function delta_bath_array(x,dmft_bath_,axis) result(Delta)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    complex(8),dimension(:),intent(in)                                :: x          !complex  array for the frequency
    type(effective_bath)                                              :: dmft_bath_ !the current :f:var:`effective_bath` instance
    character(len=*),optional                                         :: axis       !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x))               :: Delta
    character(len=4)                                                  :: axis_
    !
    axis_="mats";if(present(axis))axis_=str(axis)
    !
    Delta=zero
    !
    select case(bath_type)
    case default;stop "delta_bath_array error: bath_type not supported"
    case("normal")  ;Delta = delta_bath_array_normal(x,dmft_bath_,axis_)
    case("hybrid")  ;Delta = delta_bath_array_hybrid(x,dmft_bath_,axis_)
    case("replica") ;Delta = delta_bath_array_replica(x,dmft_bath_,axis_)
    case("general") ;Delta = delta_bath_array_general(x,dmft_bath_,axis_)
    end select
    !
  end function delta_bath_array



  include "delta_functions/fdelta_normal.f90"
  include "delta_functions/fdelta_hybrid.f90"
  include "delta_functions/fdelta_replica.f90"
  include "delta_functions/fdelta_general.f90"
  function fdelta_bath_array(x,dmft_bath_,axis) result(Fdelta)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    complex(8),dimension(:),intent(in)                                :: x          !complex  array for the frequency
    type(effective_bath)                                              :: dmft_bath_ !the current :f:var:`effective_bath` instance
    character(len=*),optional                                         :: axis       !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x))               :: Fdelta
    character(len=4)                                                  :: axis_
    !
    axis_="mats";if(present(axis))axis_=str(axis)
    !
    Fdelta=zero
    !    
    if(ed_mode/="superc")then
       print*,"fdelta_bath_array WARNING: called with ed_mode != superc. Return zero"
       call sleep(2000)
       return
    endif
    !
    select case(bath_type)
    case default;stop "fdelta_bath_array error: bath_type not supported"
    case("normal")  ;FDelta = fdelta_bath_array_normal(x,dmft_bath_,axis_)
    case("hybrid")  ;FDelta = fdelta_bath_array_hybrid(x,dmft_bath_,axis_)
    case("replica") ;FDelta = fdelta_bath_array_replica(x,dmft_bath_,axis_)
    case("general") ;FDelta = fdelta_bath_array_general(x,dmft_bath_,axis_)
    end select
    !
  end function fdelta_bath_array






  include "g0and_functions/g0and_normal.f90"
  include "g0and_functions/g0and_hyrege.f90"
  function g0and_bath_array(x,dmft_bath_,axis) result(G0and)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    complex(8),dimension(:),intent(in)                  :: x          !complex  array for the frequency
    type(effective_bath)                                :: dmft_bath_ !the current :f:var:`effective_bath` instance
    character(len=*),optional                           :: axis       !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis    
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
    character(len=4)                                    :: axis_
    !
    axis_="mats";if(present(axis))axis_=str(axis)
    !
    G0and = zero
    !
    select case(bath_type)
    case default;stop "G0and_bath_array error: bath_type not supported"
    case("normal")
       G0and = G0and_bath_array_normal(x,dmft_bath_,axis_)
    case("hybrid","replica","general")
       G0and = G0and_bath_array_hyrege(x,dmft_bath_,axis_)
    end select
    !
  end function g0and_bath_array


  include "g0and_functions/f0and_normal.f90"
  include "g0and_functions/f0and_hyrege.f90"
  function f0and_bath_array(x,dmft_bath_,axis) result(F0and)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    complex(8),dimension(:),intent(in)                  :: x !complex  array for the frequency
    type(effective_bath)                                :: dmft_bath_ !the current :f:var:`effective_bath` instance
    character(len=*),optional                           :: axis!string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis    
    character(len=4)                                    :: axis_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: F0and
    !
    axis_="mats";if(present(axis))axis_=str(axis)
    !
    F0and=zero
    !
    if(ed_mode/="superc")then
       print*,"f0and_bath_array WARNING: called with ed_mode != superc. Return zero"
       call sleep(2000)
       return
    endif
    !
    select case(bath_type)
    case default;stop "f0and_bath_array error: bath_type not supported"
    case("normal")
       F0and = f0and_bath_array_normal(x,dmft_bath_,axis_)
    case("hybrid","replica","general")
       F0and = f0and_bath_array_hyrege(x,dmft_bath_,axis_)
    end select
    !
  end function f0and_bath_array




  include "invg0_functions/invg0_normal.f90"
  include "invg0_functions/invg0_hyrege.f90"
  function invg0_bath_array(x,dmft_bath_,axis) result(G0and)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    complex(8),dimension(:),intent(in)                  :: x          !complex  array for the frequency
    type(effective_bath)                                :: dmft_bath_ !the current :f:var:`effective_bath` instance
    character(len=*),optional                           :: axis       !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis    
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and    
    character(len=4)                                    :: axis_
    !
    axis_="mats";if(present(axis))axis_=str(axis)
    !
    G0and = zero
    select case(bath_type)
    case default;stop "invg0_bath_array error: bath_type not supported"
    case("normal")
       G0and = invg0_bath_array_normal(x,dmft_bath_,axis_)
    case("hybrid","replica","general")
       G0and = invg0_bath_array_hyrege(x,dmft_bath_,axis_)
    end select
  end function invg0_bath_array



  include "invg0_functions/invf0_normal.f90"
  include "invg0_functions/invf0_hyrege.f90"
  function invf0_bath_array(x,dmft_bath_,axis) result(F0and)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    complex(8),dimension(:),intent(in)                  :: x          !complex  array for the frequency
    type(effective_bath)                                :: dmft_bath_ !the current :f:var:`effective_bath` instance
    character(len=*),optional                           :: axis       !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis    
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: F0and
    character(len=4)                                    :: axis_
    !
    axis_="mats";if(present(axis))axis_=str(axis)
    !
    F0and=zero
    !
    if(ed_mode/="superc")then
       print*,"f0and_bath_array WARNING: called with ed_mode != superc. Return zero"
       call sleep(2000)
       return
    endif
    !
    select case(bath_type)
    case default;stop "f0and_bath_array error: bath_type not supported"
    case("normal")
       F0and = invf0_bath_array_normal(x,dmft_bath_,axis_)
    case("hybrid","replica","general")
       F0and = invf0_bath_array_hyrege(x,dmft_bath_,axis_)
    end select
    !
  end function invf0_bath_array











  function zeta_superc(x,mu,axis) result(zeta)
    complex(8),dimension(:)                        :: x
    real(8)                                        :: mu
    character(len=*)                               :: axis
    complex(8),dimension(2*Nspin*Norb,size(x))     :: zeta
    integer                                        :: iorb,N,L
    N = Nspin*Norb
    L = size(x)
    do concurrent(iorb=1:N)
       select case(axis)
       case default
          zeta(iorb ,1:L)   = x(1:L) + mu
          zeta(N+iorb,1:L)  = x(1:L) - mu
       case ('real')
          zeta(iorb ,1:L)   = x(1:L) + mu
          zeta(N+iorb,1:L) = -conjg(x(L:1:-1) + mu)
       end select
    enddo
  end function zeta_superc


END MODULE ED_BATH_FUNCTIONS
