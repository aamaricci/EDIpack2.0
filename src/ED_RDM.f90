MODULE ED_RDM
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE
  USE ED_IO
  USE ED_BATH
  USE ED_SETUP
  USE ED_HAMILTONIAN
  USE ED_AUX_FUNX
  !
  USE ED_RDM_NORMAL
  USE ED_RDM_SUPERC
  USE ED_RDM_NONSU2
  !
  implicit none
  private 

  public :: rdm_impurity

contains


  subroutine rdm_impurity()
    ! 
    ! Calculate the impurity RDM calling the correct procedure according to the value of :f:var:`ed_mode` .
    ! Write the resulting matrix to a plain-text file.
    !
    ! * :code:`normal` : :f:func:`imp_rdm_normal`
    ! * :code:`superc` : :f:func:`imp_rdm_superc`
    ! * :code:`nonsu2` : :f:func:`imp_rdm_nonsu2`
    !
    if(.not.rdm_flag)return
    write(LOGfile,"(A)")"Get RDM:"
    select case(ed_mode)
    case default  ;call imp_rdm_normal()
    case("superc");call imp_rdm_superc()
    case("nonsu2");call imp_rdm_nonsu2()
    end select
    !
    call ed_print_dm(impurity_density_matrix,size(impurity_density_matrix,1))
  end subroutine rdm_impurity



end MODULE ED_RDM
