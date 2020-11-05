  !diagonal bath hamiltonian: +energy of the bath=\sum_a=1,Norb\sum_{l=1,Nbath}\e^a_l n^a_l
  select case(bath_type)
  case default 
     htmp=zero
     do iorb=1,size(dmft_bath%e,2)
        do kp=1,Nbath
           ialfa=getBathStride(iorb,kp)
           htmp =htmp + dmft_bath%e(1,iorb,kp)*ib(ialfa)        !UP
           htmp =htmp + dmft_bath%e(Nspin,iorb,kp)*ib(ialfa+Ns) !DW
        enddo
     enddo
     !
     i = j
     hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(i)
     !
  case ("replica")
     htmp=zero
     do kp=1,Nbath
        do iorb=1,Norb
           ialfa = getBathStride(iorb,kp)
           htmp = htmp + bath_diag(1    ,iorb,kp)*ib(ialfa)    !UP
           htmp = htmp + bath_diag(Nspin,iorb,kp)*ib(ialfa+Ns) !DW
        enddo
     enddo
     !
     i = j
     hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(i)
     !
     !off-diagonal elements
     do kp=1,Nbath
        do iorb=1,Norb
           do jorb=1,Norb
              !UP
              ialfa = getBathStride(iorb,kp)
              ibeta = getBathStride(jorb,kp)
              Jcondition = &
                   (hbath_tmp(1,1,iorb,jorb,kp)/=zero) .AND. &
                   (ib(ibeta)==1) .AND. (ib(ialfa)==0)
              if (Jcondition)then
                 call c(ibeta,m,k1,sg1)
                 call cdg(ialfa,k1,k2,sg2)
                 i = binary_search(Hsector%H(1)%map,k2)
                 htmp = conjg(hbath_tmp(1,1,iorb,jorb,kp))*sg1*sg2
                 !
                 hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
                 !
              endif
              !DW
              ialfa = getBathStride(iorb,kp) + Ns
              ibeta = getBathStride(jorb,kp) + Ns
              Jcondition = &
                   (hbath_tmp(Nspin,Nspin,iorb,jorb,kp)/=zero) .AND. &
                   (ib(ibeta)==1)  .AND. (ib(ialfa)==0)
              if (Jcondition)then
                 call c(ibeta,m,k1,sg1)
                 call cdg(ialfa,k1,k2,sg2)
                 i = binary_search(Hsector%H(1)%map,k2)
                 htmp = conjg(hbath_tmp(Nspin,Nspin,iorb,jorb,kp))*sg1*sg2
                 !
                 hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
                 !
              endif
           enddo
        enddo
     enddo
     !     
  end select



  !anomalous pair-creation/destruction
  do iorb=1,size(dmft_bath%e,2)
     do kp=1,Nbath
        ms=getBathStride(iorb,kp)
        !\Delta_l c_{\up,ms} c_{\dw,ms}
        if( (dmft_bath%d(1,iorb,kp)/=0d0) .AND. (ib(ms)==1) .AND. (ib(ms+Ns)==1) )then
           call c(ms,m,k1,sg1)
           call c(ms+Ns,k1,k2,sg2)
           i = binary_search(Hsector%H(1)%map,k2)
           htmp=one*dmft_bath%d(1,iorb,kp)*sg1*sg2
           !
           hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
           !
        endif
        !\Delta_l cdg_{\up,ms} cdg_{\dw,ms}
        if( (dmft_bath%d(1,iorb,kp)/=0d0) .AND. (ib(ms)==0) .AND. (ib(ms+Ns)==0) )then
           call cdg(ms+Ns,m,k1,sg1)
           call cdg(ms,k1,k2,sg2)
           i=binary_search(Hsector%H(1)%map,k2)
           htmp=one*dmft_bath%d(1,iorb,kp)*sg1*sg2 !
           !
           hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
           !
        endif
     enddo
  enddo


