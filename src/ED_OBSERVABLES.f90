MODULE ED_OBSERVABLES
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: str,free_unit,reg,free_units,txtfy
  USE SF_LINALG,  only: inv,eigh,eye
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_IO                     !< this contains the routine to print GF,Sigma and G0
  USE ED_EIGENSPACE
  USE ED_BATH
  USE ED_BATH_FUNCTIONS
  USE ED_SETUP
  USE ED_HAMILTONIAN
  USE ED_AUX_FUNX
  !
  USE ED_OBSERVABLES_NORMAL
  USE ED_OBSERVABLES_SUPERC
  USE ED_OBSERVABLES_NONSU2
  !
  implicit none
  private 

  public :: observables_impurity
  public :: local_energy_impurity

contains


  subroutine observables_impurity()
    write(LOGfile,"(A)")"Get observables:"
    select case(ed_mode)
    case default  ;call observables_normal()
    case("superc");call observables_superc()
    case("nonsu2");call observables_nonsu2()
    end select
  end subroutine observables_impurity


  subroutine local_energy_impurity()
    write(LOGfile,"(A)")"Get local energy:"
    select case(ed_mode)
    case default  ;call local_energy_normal()
    case("superc");call local_energy_superc()
    case("nonsu2");call local_energy_nonsu2()
    end select
  end subroutine local_energy_impurity


end MODULE ED_OBSERVABLES
