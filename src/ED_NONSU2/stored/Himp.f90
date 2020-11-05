  do i=MpiIstart,MpiIend
     m  = Hsector%H(1)%map(i)
     ib = bdecomp(m,2*Ns)
     !
     do iorb=1,Norb
        nup(iorb)=dble(ib(iorb))
        ndw(iorb)=dble(ib(iorb+Ns))
     enddo

     !> H_Imp: Diagonal Elements, i.e. local part
     htmp = zero
     do iorb=1,Norb
        htmp = htmp + impHloc(1,1,iorb,iorb)*nup(iorb)
        htmp = htmp + impHloc(Nspin,Nspin,iorb,iorb)*ndw(iorb)
        htmp = htmp - xmu*(nup(iorb)+ndw(iorb))
     enddo
     !
     select case(MpiStatus)
     case (.true.)
        call sp_insert_element(MpiComm,spH0,htmp,i,i)
     case (.false.)
        call sp_insert_element(spH0,htmp,i,i)
     end select
     !

     !Off-diagonal elements, i.e. non-local part
     !1. same spin:
     do iorb=1,Norb
        do jorb=1,Norb
           !this loop considers only the orbital off-diagonal terms
           !because iorb=jorb can not have simultaneously
           !occupation 0 and 1, as required by this if Jcondition:
           !UP
           Jcondition = &
                (impHloc(1,1,iorb,jorb)/=zero) .AND. &
                (ib(jorb)==1) .AND. (ib(iorb)==0)
           if (Jcondition) then
              call c(jorb,m,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              j = binary_search(Hsector%H(1)%map,k2)
              htmp = conjg(impHloc(1,1,iorb,jorb))*sg1*sg2
              !
              select case(MpiStatus)
              case (.true.)
                 call sp_insert_element(MpiComm,spH0,htmp,i,j)
              case (.false.)
                 call sp_insert_element(spH0,htmp,i,j)
              end select
              !
           endif
           !DW
           Jcondition = &
                (impHloc(Nspin,Nspin,iorb,jorb)/=zero) .AND. &
                (ib(jorb+Ns)==1) .AND. (ib(iorb+Ns)==0)
           if (Jcondition) then
              call c(jorb+Ns,m,k1,sg1)
              call cdg(iorb+Ns,k1,k2,sg2)
              j = binary_search(Hsector%H(1)%map,k2)
              htmp = conjg(impHloc(Nspin,Nspin,iorb,jorb))*sg1*sg2
              !
              select case(MpiStatus)
              case (.true.)
                 call sp_insert_element(MpiComm,spH0,htmp,i,j)
              case (.false.)
                 call sp_insert_element(spH0,htmp,i,j)
              end select
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
                 j = binary_search(Hsector%H(1)%map,k2)
                 htmp = conjg(impHloc(ispin,jspin,iorb,jorb))*sg1*sg2
                 !
                 select case(MpiStatus)
                 case (.true.)
                    call sp_insert_element(MpiComm,spH0,htmp,i,j)
                 case (.false.)
                    call sp_insert_element(spH0,htmp,i,j)
                 end select
                 !
              endif
              !
           enddo
        enddo
     enddo

  enddo
