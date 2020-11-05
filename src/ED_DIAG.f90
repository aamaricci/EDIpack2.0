MODULE ED_DIAG
  USE ED_DIAG_NORMAL
  USE ED_DIAG_SUPERC
  USE ED_DIAG_NONSU2
  !
  implicit none
  private

  !>Diag hamiltonian
  public  :: diagonalize_impurity_normal
  public  :: diagonalize_impurity_superc
  public  :: diagonalize_impurity_nonsu2


end MODULE ED_DIAG
