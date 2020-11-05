
  subroutine ed_get_quantum_SOC_operators_single()
    implicit none
    complex(8),allocatable            ::  Simp(:,:,:)
    complex(8),allocatable            ::  Limp(:,:,:)
    complex(8),allocatable            ::  Jimp(:),Jimp_sq(:)
    complex(8)                        ::  LSimp
    complex(8),allocatable            ::  rho_so(:,:),U(:,:),Udag(:,:)
    complex(8),allocatable            ::  rho_nn(:,:,:,:)
    integer                           ::  unit_
    integer                           ::  iorb,ispin,jorb,jspin,io,ibath
    !
    if(Norb/=3)  stop "SOC_operators implemented only for 3 orbitals"
    if(Nspin/=2) stop "SOC_operators implemented only for 2 spins"
    !
    if(allocated(Simp))   deallocate(Simp)   ;allocate(Simp(3,Norb,Norb))            ;Simp=zero
    if(allocated(Limp))   deallocate(Limp)   ;allocate(Limp(3,Nspin,Nspin))          ;Limp=zero
    if(allocated(Jimp))   deallocate(Jimp)   ;allocate(Jimp(3))                      ;Jimp=zero
    if(allocated(Jimp_sq))deallocate(Jimp_sq);allocate(Jimp_sq(3))                   ;Jimp_sq=zero
    !
    !rotations definition
    if(allocated(U))      deallocate(U)      ;allocate(U(Nspin*Norb,Nspin*Norb))     ;U=zero
    if(allocated(Udag))   deallocate(Udag)   ;allocate(Udag(Nspin*Norb,Nspin*Norb))  ;Udag=zero
    if(allocated(rho_so)) deallocate(rho_so) ;allocate(rho_so(Nspin*Norb,Nspin*Norb));rho_so=zero
    if(allocated(rho_nn)) deallocate(rho_nn) ;allocate(rho_nn(Nspin,Nspin,Norb,Norb));rho_nn=zero
    !
    if(bath_type=="replica".and.(.not.Jz_basis))then
       !
       !impurity dm in {t2g,Sz}. Rotation U: 1
       U=eye(Nspin*Norb)
       Udag=transpose(conjg(U))
       !
    elseif(bath_type=="replica".and.Jz_basis)then
       !
       !impurity dm in {Lz,Sz}. Rotation U: {Lz,Sz}-->{t2g,Sz}
       U=transpose(conjg(orbital_Lz_rotation_NorbNspin()))
       Udag=transpose(conjg(U))
       !
    elseif(bath_type=="normal")then
       !
       !impurity dm in {J,jz}. Rotation U: {J,jz}-->{t2g,Sz}
       U=transpose(conjg(atomic_SOC_rotation()))
       Udag=transpose(conjg(U))
       !
    endif
    !
    rho_so=zero;rho_nn=zero
    rho_so=nn2so_reshape(imp_density_matrix,Nspin,Norb)
    rho_so=matmul(Udag,matmul(rho_so,U))
    rho_nn=so2nn_reshape(rho_so,Nspin,Norb)
    !
    !#####################################################
    !#                  < S(iorb,jorb) >                 #
    !#####################################################
    !
    ! Sx =    [ <c+_up,c_dw> + <c+_dw,c_up> ]_(iorb,jorb)
    ! Sy = xi*[ <c+_dw,c_up> - <c+_up,c_dw> ]_(iorb,jorb)
    ! Sz =    [ <c+_up,c_up> - <c+_dw,c_dw> ]_(iorb,jorb)
    !
    do iorb=1,Norb
       do jorb=1,Norb
          Simp(1,iorb,jorb) = 0.5d0*( rho_nn(1,2,iorb,jorb) + rho_nn(2,1,iorb,jorb) )
          Simp(2,iorb,jorb) = 0.5d0*( rho_nn(2,1,iorb,jorb) - rho_nn(1,2,iorb,jorb) )*xi
          Simp(3,iorb,jorb) = 0.5d0*( rho_nn(1,1,iorb,jorb) - rho_nn(2,2,iorb,jorb) )
       enddo
    enddo
    !
    !#####################################################
    !#                 < L(ispin,jspin) >                #
    !#####################################################
    ! 1=yz 2=zx 3=xy
    ! Lx = xi*[ <c+_3,c_2> - <c+_2,c_3> ]_(ispin,jspin)
    ! Ly = xi*[ <c+_1,c_3> - <c+_3,c_1> ]_(ispin,jspin)
    ! Lz = xi*[ <c+_2,c_1> - <c+_1,c_2> ]_(ispin,jspin)
    !
    do ispin=1,Nspin
       do jspin=1,Nspin
          Limp(1,ispin,jspin) = ( rho_nn(ispin,jspin,3,2) - rho_nn(ispin,jspin,2,3) )*xi
          Limp(2,ispin,jspin) = ( rho_nn(ispin,jspin,1,3) - rho_nn(ispin,jspin,3,1) )*xi
          Limp(3,ispin,jspin) = ( rho_nn(ispin,jspin,2,1) - rho_nn(ispin,jspin,1,2) )*xi
       enddo
    enddo
    !
    !#####################################################
    !#                    < j_{x,y,z} >                  #
    !#####################################################
    !
    jimp(1) = trace(matmul(rho_so,atomic_j("x")))
    jimp(2) = trace(matmul(rho_so,atomic_j("y")))
    jimp(3) = trace(matmul(rho_so,atomic_j("z")))
    !
    !#####################################################
    !#                   < j^2_{x,y,z} >                 #
    !#####################################################
    !
    jimp_sq(1) = trace(matmul(rho_so,matmul(atomic_j("x"),atomic_j("x"))))
    jimp_sq(2) = trace(matmul(rho_so,matmul(atomic_j("y"),atomic_j("y"))))
    jimp_sq(3) = trace(matmul(rho_so,matmul(atomic_j("z"),atomic_j("z"))))
    !
    !#####################################################
    !#                        < LS >                     #
    !#####################################################
    !
    LSimp = trace(matmul(rho_so,atomic_SOC()))
    !
    write(LOGfile,"(A,10f18.12,A)") " Ji   = ",(real(jimp(io)),io=1,3)
    !
    call print_operators(Simp,Limp,Jimp,Jimp_sq,LSimp,1)
    !
  end subroutine ed_get_quantum_SOC_operators_single




  subroutine ed_get_quantum_SOC_operators_lattice()
    implicit none
    complex(8),allocatable            ::  Simp(:,:,:,:)
    complex(8),allocatable            ::  Limp(:,:,:,:)
    complex(8),allocatable            ::  Jimp(:,:),Jimp_sq(:,:)
    complex(8),allocatable            ::  LSimp(:)
    complex(8),allocatable            ::  rho_so(:,:,:),U(:,:),Udag(:,:)
    complex(8),allocatable            ::  rho_nn(:,:,:,:,:)
    integer                           ::  unit_
    integer                           ::  iorb,ispin,jorb,jspin,io,ibath
    integer                           ::  ilat,Nlat
    complex(8),allocatable            ::  Simp_tmp(:,:,:)
    complex(8),allocatable            ::  Limp_tmp(:,:,:)
    complex(8),allocatable            ::  Jimp_tmp(:),Jimp_sq_tmp(:)
    complex(8)                        ::  LSimp_tmp
    !
    if(Norb/=3)  stop "SOC_operators implemented only for 3 orbitals"
    if(Nspin/=2) stop "SOC_operators implemented only for 2 spins"
    !
    Nlat=size(imp_density_matrix_ineq,1)
    !
    if(allocated(Simp))       deallocate(Simp)       ;allocate(Simp(Nlat,3,Norb,Norb))            ;Simp=zero
    if(allocated(Limp))       deallocate(Limp)       ;allocate(Limp(Nlat,3,Nspin,Nspin))          ;Limp=zero
    if(allocated(Jimp))       deallocate(Jimp)       ;allocate(Jimp(Nlat,3))                      ;Jimp=zero
    if(allocated(Jimp_sq))    deallocate(Jimp_sq)    ;allocate(Jimp_sq(Nlat,3))                   ;Jimp_sq=zero
    if(allocated(LSimp))      deallocate(LSimp)      ;allocate(LSimp(Nlat))                       ;LSimp=zero
    !
    if(allocated(Simp_tmp))   deallocate(Simp_tmp)   ;allocate(Simp_tmp(3,Norb,Norb))             ;Simp_tmp=zero
    if(allocated(Limp_tmp))   deallocate(Limp_tmp)   ;allocate(Limp_tmp(3,Nspin,Nspin))           ;Limp_tmp=zero
    if(allocated(Jimp_tmp))   deallocate(Jimp_tmp)   ;allocate(Jimp_tmp(3))                       ;Jimp_tmp=zero
    if(allocated(Jimp_sq_tmp))deallocate(Jimp_sq_tmp);allocate(Jimp_sq_tmp(3))                    ;Jimp_sq_tmp=zero
    !
    !rotations definition
    if(allocated(U))          deallocate(U)          ;allocate(U(Nspin*Norb,Nspin*Norb))          ;U=zero
    if(allocated(Udag))       deallocate(Udag)       ;allocate(Udag(Nspin*Norb,Nspin*Norb))       ;Udag=zero
    if(allocated(rho_so))     deallocate(rho_so)     ;allocate(rho_so(Nlat,Nspin*Norb,Nspin*Norb));rho_so=zero
    if(allocated(rho_nn))     deallocate(rho_nn)     ;allocate(rho_nn(Nlat,Nspin,Nspin,Norb,Norb));rho_nn=zero
    !
    if(bath_type=="replica".and.(.not.Jz_basis))then
       !
       !impurity dm in {t2g,Sz}. Rotation U: 1
       U=eye(Nspin*Norb)
       Udag=transpose(conjg(U))
       !
    elseif(bath_type=="replica".and.Jz_basis)then
       !
       !impurity dm in {Lz,Sz}. Rotation U: {Lz,Sz}-->{t2g,Sz}
       U=transpose(conjg(orbital_Lz_rotation_NorbNspin()))
       Udag=transpose(conjg(U))
       !
    elseif(bath_type=="normal")then
       !
       !impurity dm in {J,jz}. Rotation U: {J,jz}-->{t2g,Sz}
       U=transpose(conjg(atomic_SOC_rotation()))
       Udag=transpose(conjg(U))
       !
    endif
    !
    do ilat=1,Nlat
       !
       rho_so(ilat,:,:)=nn2so_reshape(imp_density_matrix_ineq(ilat,:,:,:,:),Nspin,Norb)
       rho_so(ilat,:,:)=matmul(Udag,matmul(rho_so(ilat,:,:),U))
       rho_nn(ilat,:,:,:,:)=so2nn_reshape(rho_so(ilat,:,:),Nspin,Norb)
       !
       !#####################################################
       !#                  < S(iorb,jorb) >                 #
       !#####################################################
       !
       ! Sx =    [ <c+_up,c_dw> + <c+_dw,c_up> ]_(iorb,jorb)
       ! Sy = xi*[ <c+_dw,c_up> - <c+_up,c_dw> ]_(iorb,jorb)
       ! Sz =    [ <c+_up,c_up> - <c+_dw,c_dw> ]_(iorb,jorb)
       !
       do iorb=1,Norb
          do jorb=1,Norb
             Simp(ilat,1,iorb,jorb) = 0.5d0*( rho_nn(ilat,1,2,iorb,jorb) + rho_nn(ilat,2,1,iorb,jorb) )
             Simp(ilat,2,iorb,jorb) = 0.5d0*( rho_nn(ilat,2,1,iorb,jorb) - rho_nn(ilat,1,2,iorb,jorb) )*xi
             Simp(ilat,3,iorb,jorb) = 0.5d0*( rho_nn(ilat,1,1,iorb,jorb) - rho_nn(ilat,2,2,iorb,jorb) )
          enddo
       enddo
       !
       !#####################################################
       !#                 < L(ispin,jspin) >                #
       !#####################################################
       ! 1=yz 2=zx 3=xy
       ! Lx = xi*[ <c+_3,c_2> - <c+_2,c_3> ]_(ispin,jspin)
       ! Ly = xi*[ <c+_1,c_3> - <c+_3,c_1> ]_(ispin,jspin)
       ! Lz = xi*[ <c+_2,c_1> - <c+_1,c_2> ]_(ispin,jspin)
       !
       do ispin=1,Nspin
          do jspin=1,Nspin
             Limp(ilat,1,ispin,jspin) = ( rho_nn(ilat,ispin,jspin,3,2) - rho_nn(ilat,ispin,jspin,2,3) )*xi
             Limp(ilat,2,ispin,jspin) = ( rho_nn(ilat,ispin,jspin,1,3) - rho_nn(ilat,ispin,jspin,3,1) )*xi
             Limp(ilat,3,ispin,jspin) = ( rho_nn(ilat,ispin,jspin,2,1) - rho_nn(ilat,ispin,jspin,1,2) )*xi
          enddo
       enddo
       !
       !#####################################################
       !#                    < j_{x,y,z} >                  #
       !#####################################################
       !
       jimp(ilat,1) = trace(matmul(rho_so(ilat,:,:),atomic_j("x")))
       jimp(ilat,2) = trace(matmul(rho_so(ilat,:,:),atomic_j("y")))
       jimp(ilat,3) = trace(matmul(rho_so(ilat,:,:),atomic_j("z")))
       !
       !#####################################################
       !#                   < j^2_{x,y,z} >                 #
       !#####################################################
       !
       jimp_sq(ilat,1) = trace(matmul(rho_so(ilat,:,:),matmul(atomic_j("x"),atomic_j("x"))))
       jimp_sq(ilat,2) = trace(matmul(rho_so(ilat,:,:),matmul(atomic_j("y"),atomic_j("y"))))
       jimp_sq(ilat,3) = trace(matmul(rho_so(ilat,:,:),matmul(atomic_j("z"),atomic_j("z"))))
       !
       !#####################################################
       !#                        < LS >                     #
       !#####################################################
       !
       LSimp(ilat) = trace(matmul(rho_so(ilat,:,:),atomic_SOC()))
       !
       write(LOGfile,"(A8,I3,A10,10f18.12,A)") "site:",ilat," Ji   = ",(real(jimp(ilat,io)),io=1,3)
       !
       Simp_tmp    = zero ; Simp_tmp    = Simp(ilat,:,:,:)
       Limp_tmp    = zero ; Limp_tmp    = Limp(ilat,:,:,:)
       Jimp_tmp    = zero ; Jimp_tmp    = Jimp(ilat,:)
       Jimp_sq_tmp = zero ; Jimp_sq_tmp = Jimp_sq(ilat,:)
       LSimp_tmp   = zero ; LSimp_tmp   = LSimp(ilat)
       call print_operators(Simp_tmp,Limp_tmp,Jimp_tmp,Jimp_sq_tmp,LSimp_tmp,ilat)
       !
    enddo
    !
  end subroutine ed_get_quantum_SOC_operators_lattice


  subroutine print_operators(S,L,J,Jsq,LS,ndx)
    implicit none
    complex(8),allocatable,intent(in)   :: S(:,:,:)
    complex(8),allocatable,intent(in)   :: L(:,:,:)
    complex(8),allocatable,intent(in)   :: J(:)
    complex(8),allocatable,intent(in)   :: Jsq(:)
    complex(8)            ,intent(in)   :: LS
    integer                             :: ndx
    !internal
    integer                             :: unit_
    character(len=11)                   :: suffix
    integer                             :: iorb,jorb,ispin,jspin,io,jo
    !
    suffix="S_imp_"//reg(str(ndx))//".dat"
    unit_ = free_unit()
    open(unit=unit_,file=suffix,status='unknown',position='rewind',action='write',form='formatted')
    write(unit_,'(A)')"# Re{Sx_(iorb,jorb)}, Im{Sx_(iorb,jorb)}"
    do iorb=1,Norb
       write(unit_,'(30(F20.12,1X))') (real(S(1,iorb,jorb)),jorb=1,Norb),(aimag(S(1,iorb,jorb)),jorb=1,Norb)
    enddo
    write(unit_,*)
    write(unit_,'(A)')"# Re{Sy_(iorb,jorb)}, Im{Sy_(iorb,jorb)}"
    do iorb=1,Norb
       write(unit_,'(30(F20.12,1X))') (real(S(2,iorb,jorb)),jorb=1,Norb),(aimag(S(2,iorb,jorb)),jorb=1,Norb)
    enddo
    write(unit_,*)
    write(unit_,'(A)')"# Re{Sz_(iorb,jorb)}, Im{Sz_(iorb,jorb)}"
    do iorb=1,Norb
       write(unit_,'(30(F20.12,1X))') (real(S(3,iorb,jorb)),jorb=1,Norb),(aimag(S(3,iorb,jorb)),jorb=1,Norb)
    enddo
    write(unit_,*)
    write(unit_,'(30(a20,1X))')"#1-Re{Tr[Sx]}","2-Im{Tr[Sx]}","3-Re{Tr[Sy]}","4-Im{Tr[Sy]}","5-Re{Tr[Sz]}","6-Im{Tr[Sz]}"
    write(unit_,'(30(F20.12,1X))') real(trace(S(1,:,:))),aimag(trace(S(1,:,:))),&
         real(trace(S(2,:,:))),aimag(trace(S(2,:,:))),&
         real(trace(S(3,:,:))),aimag(trace(S(3,:,:)))
    close(unit_)
    !
    suffix="L_imp_"//reg(str(ndx))//".dat"
    unit_ = free_unit()
    open(unit=unit_,file=suffix,status='unknown',position='rewind',action='write',form='formatted')
    write(unit_,'(A)')"# Re{Lx_(ipin,jspin)}, Im{Lx_(ipin,jspin)}"
    do ispin=1,Nspin
       write(unit_,'(30(F20.12,1X))') (real(L(1,ispin,jspin)),jspin=1,Nspin),(aimag(L(1,ispin,jspin)),jspin=1,Nspin)
    enddo
    write(unit_,*)
    write(unit_,'(A)')"# Re{Ly_(ipin,jspin)}, Im{Ly_(ipin,jspin)}"
    do ispin=1,Nspin
       write(unit_,'(30(F20.12,1X))') (real(L(2,ispin,jspin)),jspin=1,Nspin),(aimag(L(2,ispin,jspin)),jspin=1,Nspin)
    enddo
    write(unit_,*)
    write(unit_,'(A)')"# Re{Lz_(ipin,jspin)}, Im{Lz_(ipin,jspin)}"
    do ispin=1,Nspin
       write(unit_,'(30(F20.12,1X))') (real(L(3,ispin,jspin)),jspin=1,Nspin),(aimag(L(3,ispin,jspin)),jspin=1,Nspin)
    enddo
    write(unit_,*)
    write(unit_,'(30(a20,1X))')"#1-Re{Tr[Lx]}","2-Im{Tr[Lx]}","3-Re{Tr[Ly]}","4-Im{Tr[ly]}","5-Re{Tr[Lz]}","6-Im{Tr[Lz]}"
    write(unit_,'(30(F20.12,1X))') real(trace(L(1,:,:))),aimag(trace(L(1,:,:))),&
         real(trace(L(2,:,:))),aimag(trace(L(2,:,:))),&
         real(trace(L(3,:,:))),aimag(trace(L(3,:,:)))
    close(unit_)
    !
    suffix="J_imp_"//reg(str(ndx))//".dat"
    unit_ = free_unit()
    open(unit=unit_,file=suffix,status='unknown',position='rewind',action='write',form='formatted')
    write(unit_,'(30(a20,1X))') "# 1-Re{jx}   "," 2-Im{jx}   "," 3-Re{jy}   "," 4-Im{jy}   "," 5-Re{jz}   "," 6-Im{jz}   ", &
         " 7-Re{jx_sq}"," 8-Im{jx_sq}"," 9-Re{jy_sq}","10-Im{jy_sq}","11-Re{jz_sq}","12-Im{jz_sq}", &
         "13-Re{L.S}","14-Im{L.S}"
    write(unit_,'(30(F20.12,1X))') real(J(1))   ,aimag(J(1))   ,real(J(2))   ,aimag(J(2))   ,real(J(3))   ,aimag(J(3))   , &
         real(Jsq(1)),aimag(Jsq(1)),real(Jsq(2)),aimag(Jsq(2)),real(Jsq(3)),aimag(Jsq(3)), &
         real(LS),aimag(LS)!,(real(LSbth(io)),io=1,Nbath)
    close(unit_)
    !
  end subroutine print_operators

