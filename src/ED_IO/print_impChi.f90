! SPIN-SPIN
subroutine print_chi_spin(Self,axis)
  complex(8),dimension(:,:,:) :: Self
  character(len=1)            :: axis
  integer                     :: i,j,iorb,jorb
  integer                     :: L,i,ispin,isign
  character(len=20)           :: suffix
  call allocate_grids
  do iorb=1,Norb
     do jorb=1,Norb
        suffix="_l"//str(iorb)//str(jorb)
        select case(to_lower(axis))
        case default;stop "print_chi_spib error: axis not supported"
        case("m");call splot("spinChi"//str(suffix)//"_iv"//reg(ed_file_suffix)//".ed",vm,Self(iorb,jorb,:))
        case("r");call splot("spinChi"//str(suffix)//"_realw"//reg(ed_file_suffix)//".ed",vr,Self(iorb,jorb,:))
        case("t");call splot("spinChi"//str(suffix)//"_tau"//reg(ed_file_suffix)//".ed",tau,Self(iorb,jorb,:))
        end select
     enddo
  enddo
  call deallocate_grids
end subroutine print_chi_spin

! DENSITY-DENSITY
subroutine print_chi_dens(Self,axis)
  complex(8),dimension(:,:,:,:,:) :: Self
  character(len=1)                :: axis
  integer                         :: i,j,iorb,jorb
  integer                         :: L,i,ispin,isign
  character(len=20)               :: suffix
  call allocate_grids
  do iorb=1,Norb
     do jorb=1,Norb
        suffix="_l"//str(iorb)//str(jorb)
        select case(to_lower(axis))
        case default;stop "print_chi_spib error: axis not supported"
        case("m");call splot("spinChi"//str(suffix)//"_iv"//reg(ed_file_suffix)//".ed",vm,Self(iorb,jorb,:))
        case("r");call splot("spinChi"//str(suffix)//"_realw"//reg(ed_file_suffix)//".ed",vr,Self(iorb,jorb,:))
        case("t");call splot("spinChi"//str(suffix)//"_tau"//reg(ed_file_suffix)//".ed",tau,Self(iorb,jorb,:))
        end select
     enddo
  enddo
  call deallocate_grids
end subroutine print_chi_dens

! PAIR-PAIR
subroutine print_chi_pair(Self,axis)
  complex(8),dimension(:,:,:,:,:) :: Self
  character(len=1)                :: axis
  integer                         :: i,j,iorb,jorb
  integer                         :: L,i,ispin,isign
  character(len=20)               :: suffix
  call allocate_grids
  do iorb=1,Norb
     do jorb=1,Norb
        suffix="_l"//str(iorb)//str(jorb)
        select case(to_lower(axis))
        case default;stop "print_chi_spib error: axis not supported"
        case("m");call splot("spinChi"//str(suffix)//"_iv"//reg(ed_file_suffix)//".ed",vm,Self(iorb,jorb,:))
        case("r");call splot("spinChi"//str(suffix)//"_realw"//reg(ed_file_suffix)//".ed",vr,Self(iorb,jorb,:))
        case("t");call splot("spinChi"//str(suffix)//"_tau"//reg(ed_file_suffix)//".ed",tau,Self(iorb,jorb,:))
        end select
     enddo
  enddo
  call deallocate_grids
end subroutine print_chi_pair

! EXCITON
subroutine print_chi_exct(Self,axis)
  complex(8),dimension(:,:,:,:,:) :: Self
  character(len=1)                :: axis
  integer                         :: i,j,iorb,jorb
  integer                         :: L,i,ispin,isign
  character(len=20)               :: suffix
  call allocate_grids
  do iorb=1,Norb
     do jorb=iorb+1,Norb
        suffix="_l"//str(iorb)//str(jorb)
        select case(to_lower(axis))
        case default;stop "print_chi_spib error: axis not supported"
        case("m")
           call splot("exctChi_singlet"//str(suffix)//"_iv"//reg(ed_file_suffix)//".ed",vm,Self(1,orb,jorb,:))
           call splot("exctChi_tripletXY"//str(suffix)//"_iv"//reg(ed_file_suffix)//".ed",vm,Self(1,iorb,jorb,:))
           call splot("exctChi_tripletZ"//str(suffix)//"_iv"//reg(ed_file_suffix)//".ed",vm,Self(1,iorb,jorb,:))
        case("r")
           call splot("exctChi_singlet"//str(suffix)//"_realw"//reg(ed_file_suffix)//".ed",vr,Self(2,iorb,jorb,:))
           call splot("exctChi_tripletXY"//str(suffix)//"_realw"//reg(ed_file_suffix)//".ed",vr,Self(2,iorb,jorb,:))
           call splot("exctChi_tripletZ"//str(suffix)//"_realw"//reg(ed_file_suffix)//".ed",vr,Self(2,iorb,jorb,:))
        case("t")
           call splot("exctChi_singlet"//str(suffix)//"_tau"//reg(ed_file_suffix)//".ed",tau,Self(3,iorb,jorb,:))
           call splot("exctChi_tripletXY"//str(suffix)//"_tau"//reg(ed_file_suffix)//".ed",tau,Self(3,iorb,jorb,:))
           call splot("exctChi_tripletZ"//str(suffix)//"_tau"//reg(ed_file_suffix)//".ed",tau,Self(3,iorb,jorb,:))
        end select
     enddo
  enddo
  call deallocate_grids
end subroutine print_chi_exct

