module edi2py_bindings
    use edipack2
    use scifor
    implicit none
    
    
contains
    
    !integer to logical
    function i2l(var_integer) result (var_logical)
      use, intrinsic :: iso_c_binding
      integer        :: var_integer
      logical        :: var_logical   
    
      if (var_integer == 1) then
        var_logical = .true.
      else
        var_logical = .false.
      endif
    end function i2l
    
    !logical to integer
    function l2i(var_logical) result (var_integer)
      use, intrinsic :: iso_c_binding
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
      use, intrinsic :: iso_c_binding
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
    
    !Get Bath Type
    integer(c_int) function get_bath_type_c() result(bt) bind(c, name='get_bath_type')
      use, intrinsic :: iso_c_binding
      select case(bath_type)
        case("normal");  bt = 1
        case("hybrid");  bt = 2
        case("replica"); bt = 3
        case("general"); bt = 4
      end select
    end function get_bath_type_c
    
    !Get Bath Type
    integer(c_int) function get_ed_mode_c() result(edm) bind(c, name='get_ed_mode')
      use, intrinsic :: iso_c_binding
      select case(ed_mode)
        case("normal");  edm = 1
        case("superc");  edm = 2
        case("nonsu2");  edm = 3
      end select
    end function get_ed_mode_c
    
    !include library functions
    include "edi2py_read_input.f90"
    include "edi2py_main.f90"
    include "edi2py_bath.f90"
    include "edi2py_io.f90"
    include "edi2py_bath_fit.f90"
    include "edi2py_aux_funx.f90"

end module edi2py_bindings
