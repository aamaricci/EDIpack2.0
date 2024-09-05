module edi2py_bindings
    use edipack2
    use scifor
    use iso_c_binding
    implicit none
    
    real(c_double),dimension(5),bind(c, name="Uloc_cbind")                   :: Uloc_cbind                !local interactions
    
contains
    
    !integer to logical
    function i2l(var_integer) result (var_logical)
      integer    :: var_integer
      logical    :: var_logical   
    
      if (var_integer == 1) then
        var_logical = .true.
      else
        var_logical = .false.
      endif
    end function i2l
    
    !logical to integer
    function l2i(var_logical) result (var_integer)
      integer    :: var_integer
      logical    :: var_logical   
    
      if (var_logical) then
        var_integer = 1
      else
        var_integer = 0
      endif
    end function l2i
  
    !c string to fortran string
    subroutine c2f(c_str)
        character(kind=c_char), dimension(*),intent(IN) :: c_str
        character(len=120), allocatable                 :: f_str
        integer                                         :: length
        integer                                         :: i
        
        length=0
        f_str=" "
        do
           if (c_str(length+1) == C_NULL_CHAR) exit
           length = length + 1
        end do
        do i = 1, length
          f_str(i:i) = c_str(i)
        enddo
        f_str=trim(f_str)
    end subroutine c2f
    
    subroutine update_array(c_str1,c_str2) bind(c, name="update_array")
        character(kind=c_char), dimension(*),intent(IN) :: c_str1, c_str2
        character(len=120), allocatable                 :: f_str1, f_str2
        integer                                         :: length
        integer                                         :: i
        
        f_str1=" "
        f_str2=" "
        
        length=0
        do
           if (c_str1(length+1) == C_NULL_CHAR) exit
           length = length + 1
        end do
        do i = 1, length
          f_str1(i:i) = c_str1(i)
        enddo
        f_str1=trim(f_str1)
        
        length=0
        do
           if (c_str2(length+1) == C_NULL_CHAR) exit
           length = length + 1
        end do
        do i = 1, length
          f_str2(i:i) = c_str2(i)
        enddo
        f_str1=trim(f_str1)
        
        select case (f_str1)
          case default
            STOP "This variable is not updatable"
          case("Uloc")
            if(allocated(Uloc))deallocate(Uloc)
            allocate(Uloc(Norb))
            if(f_str2=="set") Uloc(1:Norb) = Uloc_cbind(1:Norb)
            if(f_str2=="get") Uloc_cbind(1:Norb) = Uloc(1:Norb)
        end select
    end subroutine update_array
    
    !include library functions
    include "edi2py_read_input.f90"
    include "edi2py_main.f90"
    include "edi2py_bath.f90"
    include "edi2py_io.f90"
    include "edi2py_bath_fit.f90"
    include "edi2py_aux_funx.f90"

end module edi2py_bindings
