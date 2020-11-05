!                         SPIN-SPIN
subroutine print_chi_spin
  integer                               :: i,j,iorb,jorb
  do iorb=1,Norb
     do jorb=1,Norb
        call splot("spinChi_l"//str(iorb)//str(jorb)//"_tau"//reg(ed_file_suffix)//".ed",tau,spinChi_tau(iorb,jorb,0:))
        call splot("spinChi_l"//str(iorb)//str(jorb)//"_realw"//reg(ed_file_suffix)//".ed",vr,spinChi_w(iorb,jorb,:))
        call splot("spinChi_l"//str(iorb)//str(jorb)//"_iw"//reg(ed_file_suffix)//".ed",vm,spinChi_iv(iorb,jorb,:))
     enddo
  enddo
end subroutine print_chi_spin
!                     DENSITY-DENSITY
subroutine print_chi_dens
  integer                               :: i,j,iorb,jorb
  do iorb=1,Norb
     do jorb=1,Norb
        call splot("densChi_l"//str(iorb)//str(jorb)//"_tau"//reg(ed_file_suffix)//".ed",tau,densChi_tau(iorb,jorb,0:))
        call splot("densChi_l"//str(iorb)//str(jorb)//"_realw"//reg(ed_file_suffix)//".ed",vr,densChi_w(iorb,jorb,:))
        call splot("densChi_l"//str(iorb)//str(jorb)//"_iw"//reg(ed_file_suffix)//".ed",vm,densChi_iv(iorb,jorb,:))
     enddo
  enddo
end subroutine print_chi_dens
!                     PAIR-PAIR
subroutine print_chi_pair
  integer                               :: i,j,iorb,jorb
  do iorb=1,Norb
     do jorb=1,Norb
        call splot("pairChi_l"//str(iorb)//str(jorb)//"_tau"//reg(ed_file_suffix)//".ed",tau,pairChi_tau(iorb,jorb,0:))
        call splot("pairChi_l"//str(iorb)//str(jorb)//"_realw"//reg(ed_file_suffix)//".ed",vr,pairChi_w(iorb,jorb,:))
        call splot("pairChi_l"//str(iorb)//str(jorb)//"_iw"//reg(ed_file_suffix)//".ed",vm,pairChi_iv(iorb,jorb,:))
     enddo
  enddo
end subroutine print_chi_pair
!                     EXCITON
subroutine print_chi_exct
  integer                               :: i,j,iorb,jorb
  do iorb=1,Norb
     do jorb=iorb+1,Norb
        call splot("exctChi_singlet_l"//str(iorb)//str(jorb)//"_tau"//reg(ed_file_suffix)//".ed",tau,exctChi_tau(0,iorb,jorb,0:))
        call splot("exctChi_singlet_l"//str(iorb)//str(jorb)//"_realw"//reg(ed_file_suffix)//".ed",vr,exctChi_w(0,iorb,jorb,:))
        call splot("exctChi_singlet_l"//str(iorb)//str(jorb)//"_iw"//reg(ed_file_suffix)//".ed",vm,exctChi_iv(0,iorb,jorb,:))
     enddo
  enddo
  !
  do iorb=1,Norb
     do jorb=iorb+1,Norb
        call splot("exctChi_tripletXY_l"//str(iorb)//str(jorb)//"_tau"//reg(ed_file_suffix)//".ed",tau,exctChi_tau(1,iorb,jorb,0:))
        call splot("exctChi_tripletXY_l"//str(iorb)//str(jorb)//"_realw"//reg(ed_file_suffix)//".ed",vr,exctChi_w(1,iorb,jorb,:))
        call splot("exctChi_tripletXY_l"//str(iorb)//str(jorb)//"_iw"//reg(ed_file_suffix)//".ed",vm,exctChi_iv(1,iorb,jorb,:))
     enddo
  enddo
  !
  do iorb=1,Norb
     do jorb=iorb+1,Norb
        call splot("exctChi_tripletZ_l"//str(iorb)//str(jorb)//"_tau"//reg(ed_file_suffix)//".ed",tau,exctChi_tau(2,iorb,jorb,0:))
        call splot("exctChi_tripletZ_l"//str(iorb)//str(jorb)//"_realw"//reg(ed_file_suffix)//".ed",vr,exctChi_w(2,iorb,jorb,:))
        call splot("exctChi_tripletZ_l"//str(iorb)//str(jorb)//"_iw"//reg(ed_file_suffix)//".ed",vm,exctChi_iv(2,iorb,jorb,:))
     enddo
  enddo
  !
end subroutine print_chi_exct

