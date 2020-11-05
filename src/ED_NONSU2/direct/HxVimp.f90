  !Diagonal Elements, i.e. local part
  htmp = zero
  htmp = htmp - xmu*(sum(nup)+sum(ndw))
  !
  do iorb=1,Norb
     htmp = htmp + impHloc(1,1,iorb,iorb)*nup(iorb)
     htmp = htmp + impHloc(Nspin,Nspin,iorb,iorb)*ndw(iorb)
  enddo
  !
  i = j
  hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(i)
  !
  !Off-diagonal elements, i.e. non-local part
  !1. same spin:
  do iorb=1,Norb
     do jorb=1,Norb
        !UP
        Jcondition = &
             (impHloc(1,1,iorb,jorb)/=zero) .AND. &
             (ib(jorb)==1) .AND. (ib(iorb)==0)
        if (Jcondition) then
           call c(jorb,m,k1,sg1)
           call cdg(iorb,k1,k2,sg2)
           i = binary_search(Hsector%H(1)%map,k2)
           htmp = impHloc(1,1,iorb,jorb)*sg1*sg2
           !
           hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
           !
        endif
        !DW
        Jcondition = &
             (impHloc(Nspin,Nspin,iorb,jorb)/=zero) .AND. &
             (ib(jorb+Ns)==1) .AND. (ib(iorb+Ns)==0)
        if (Jcondition) then
           call c(jorb+Ns,m,k1,sg1)
           call cdg(iorb+Ns,k1,k2,sg2)
           i = binary_search(Hsector%H(1)%map,k2)
           htmp = impHloc(Nspin,Nspin,iorb,jorb)*sg1*sg2
           !
           hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
           !
        endif
     enddo
  enddo
  !
  !2. spin-flip part (only for the nonSU2 channel!)
  do ispin=1,Nspin
     jspin = 3-ispin !ispin=1,jspin=2, ispin=2,jspin=1
     do iorb=1,Norb
        do jorb=1,Norb           
           ialfa = iorb + (ispin-1)*Ns
           ibeta = jorb + (jspin-1)*Ns
           Jcondition=&
                (impHloc(ispin,jspin,iorb,jorb)/=zero) .AND. &
                (ib(ibeta)==1) .AND. (ib(ialfa)==0)
           if(Jcondition)then
              call c(ibeta,m,k1,sg1)
              call cdg(ialfa,k1,k2,sg2)
              i = binary_search(Hsector%H(1)%map,k2)
              htmp = impHloc(ispin,jspin,iorb,jorb)*sg1*sg2
              !
              hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
              !
           endif
           !
        enddo
     enddo
  enddo

