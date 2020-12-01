MODULE ED_DIAG
  USE ED_INPUT_VARS
  USE ED_DIAG_NORMAL
  USE ED_DIAG_SUPERC
  USE ED_DIAG_NONSU2
  !
  implicit none
  private

  !>Diag hamiltonian
  public  :: diagonalize_impurity

contains

  subroutine  diagonalize_impurity()
    !
    write(LOGfile,"(A)")"Diagonalize impurity problem:"
    select case(ed_mode)
    case default  ;call diagonalize_impurity_normal()
    case("superc");call diagonalize_impurity_superc()
    case("nonsu2");call diagonalize_impurity_nonsu2()
    end select
  end subroutine diagonalize_impurity

end MODULE ED_DIAG
