!+-------------------------------------------------------------------+
!PURPOSE  : Allocate the ED bath
!+-------------------------------------------------------------------+
subroutine allocate_dmft_bath(dmft_bath_)
  type(effective_bath) :: dmft_bath_
  integer              :: Nsym,ibath
  if(dmft_bath_%status)call deallocate_dmft_bath(dmft_bath_)
  !
  select case(bath_type)
  case default
     !
     select case(ed_mode)
     case default                                 !normal [N,Sz]
        allocate(dmft_bath_%e(Nspin,Norb,Nbath))  !local energies of the bath
        allocate(dmft_bath_%v(Nspin,Norb,Nbath))  !same-spin hybridization 
     case ("superc")                              !superc [Sz] 
        allocate(dmft_bath_%e(Nspin,Norb,Nbath))  !local energies of the bath
        allocate(dmft_bath_%d(Nspin,Norb,Nbath))  !local SC order parameters the bath
        allocate(dmft_bath_%v(Nspin,Norb,Nbath))  !same-spin hybridization 
     case ("nonsu2")                              !nonsu2 [N]
        allocate(dmft_bath_%e(Nspin,Norb,Nbath))  !local energies of the bath
        allocate(dmft_bath_%v(Nspin,Norb,Nbath))  !same-spin hybridization 
        allocate(dmft_bath_%u(Nspin,Norb,Nbath))  !spin-flip hybridization
     end select
     !
  case('hybrid')
     !
     select case(ed_mode)
     case default                                 !normal  [N,Sz]
        allocate(dmft_bath_%e(Nspin,1,Nbath))     !local energies of the bath
        allocate(dmft_bath_%v(Nspin,Norb,Nbath))  !same-spin hybridization 
     case ("superc")                              !superc  [Sz]
        allocate(dmft_bath_%e(Nspin,1,Nbath))     !local energies of the bath
        allocate(dmft_bath_%d(Nspin,1,Nbath))     !local SC order parameters the bath
        allocate(dmft_bath_%v(Nspin,Norb,Nbath))  !same-spin hybridization 
     case ("nonsu2")                              !nonsu2 case [N] qn
        allocate(dmft_bath_%e(Nspin,1,Nbath))     !local energies of the bath
        allocate(dmft_bath_%v(Nspin,Norb,Nbath))  !same-spin hybridization 
        allocate(dmft_bath_%u(Nspin,Norb,Nbath))  !spin-flip hybridization
     end select
     !
  case('replica')
     !
     if(.not.allocated(lambda_impHloc))stop "lambda_impHloc not allocated in allocate_dmft_bath" !FIXME
     call deallocate_dmft_bath(dmft_bath_)     !
     allocate(dmft_bath_%item(Nbath))
     !
     Nsym=size(lambda_impHloc)
     !
     !ALLOCATE coefficients vectors
     do ibath=1,Nbath
        dmft_Bath_%item(ibath)%N_dec=Nsym
        allocate(dmft_bath_%item(ibath)%v(Nspin))
        allocate(dmft_bath_%item(ibath)%lambda(Nsym))
     enddo
     !
  end select
  !
  dmft_bath_%status=.true.
  !
end subroutine allocate_dmft_bath


!+-------------------------------------------------------------------+
!PURPOSE  : Deallocate the ED bath
!+-------------------------------------------------------------------+
subroutine deallocate_dmft_bath(dmft_bath_)
  type(effective_bath) :: dmft_bath_
  integer              :: ibath,isym
  if(.not.dmft_bath_%status)return
  if(allocated(dmft_bath_%e))   deallocate(dmft_bath_%e)
  if(allocated(dmft_bath_%d))   deallocate(dmft_bath_%d)
  if(allocated(dmft_bath_%v))   deallocate(dmft_bath_%v)
  if(allocated(dmft_bath_%u))   deallocate(dmft_bath_%u)
  if(bath_type=="replica")then
     do ibath=1,Nbath
        dmft_bath_%item(ibath)%N_dec= 0
        deallocate(dmft_bath_%item(ibath)%v)
        deallocate(dmft_bath_%item(ibath)%lambda)
     enddo
     deallocate(dmft_bath_%item)
  endif
  dmft_bath_%status=.false.
end subroutine deallocate_dmft_bath







!+-------------------------------------------------------------------+
!PURPOSE  : Reconstruct bath matrix from lambda vector
!+-------------------------------------------------------------------+
function bath_from_sym(lambdavec) result (Hbath)
  integer                                  :: Nsym,isym
  real(8),dimension(:)                     :: lambdavec
  real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hbath
  !
  Nsym=size(lambdavec)
  !
  Hbath=zero
  !
  do isym=1,Nsym
     Hbath=Hbath+lambdavec(isym)*H_Basis(isym)%O
  enddo
  !
end function bath_from_sym




!+------------------------------------------------------------------+
!PURPOSE  : Initialize the DMFT loop, builindg H parameters and/or 
!reading previous (converged) solution
!+------------------------------------------------------------------+
subroutine init_dmft_bath(dmft_bath_)
  type(effective_bath) :: dmft_bath_
  integer,dimension(Nbath) :: Nlambdas
  integer              :: i,unit,flen,Nh,isym,Nsym
  integer              :: io,jo,iorb,ispin,jorb,jspin,ibath
  logical              :: IOfile
  real(8)              :: de
  real(8)              :: rescale(Nbath),offset_b(Nbath)
  !
  if(.not.dmft_bath_%status)stop "init_dmft_bath error: bath not allocated"
  !
  if(Nbath>1)then
     rescale=linspace(HWBAND/Nbath,HWBAND,Nbath)
  else
     rescale(1)=0.d0
  endif
  !
  select case(bath_type)
  case default
     !Get energies:
     dmft_bath_%e(:,:,1)    =-hwband
     dmft_bath_%e(:,:,Nbath)= hwband
     Nh=Nbath/2
     if(mod(Nbath,2)==0.and.Nbath>=4)then
        de=hwband/max(Nh-1,1)
        dmft_bath_%e(:,:,Nh)  = -1.d-1
        dmft_bath_%e(:,:,Nh+1)=  1.d-1
        do i=2,Nh-1
           dmft_bath_%e(:,:,i)   =-hwband + (i-1)*de
           dmft_bath_%e(:,:,Nbath-i+1)= hwband - (i-1)*de
        enddo
     elseif(mod(Nbath,2)/=0.and.Nbath>=3)then
        de=hwband/Nh
        dmft_bath_%e(:,:,Nh+1)= 0d0
        do i=2,Nh
           dmft_bath_%e(:,:,i)        =-hwband + (i-1)*de
           dmft_bath_%e(:,:,Nbath-i+1)= hwband - (i-1)*de
        enddo
     endif
     !Get spin-keep yhbridizations
     do i=1,Nbath
        dmft_bath_%v(:,:,i)=max(0.1d0,1d0/sqrt(dble(Nbath)))
     enddo
     !Get SC amplitudes
     if(ed_mode=="superc")dmft_bath_%d(:,:,:)=deltasc
     !Get spin-flip hybridizations
     if(ed_mode=="nonsu2")then
        do i=1,Nbath
           dmft_bath_%u(:,:,i) = dmft_bath_%v(:,:,i)!*ed_vsf_ratio
        enddo
     endif
     !
  case('replica')
     !BATH V INITIALIZATION
     do ibath=1,Nbath
        do ispin=1,Nspin
           dmft_bath%item(ibath)%v(ispin)=max(0.1d0,1d0/sqrt(dble(Nbath)))
        enddo
     enddo
     !
     !BATH LAMBDAS INITIALIZATION
     do ibath=1,Nbath
        Nsym = dmft_bath%item(ibath)%N_dec
        do isym=1,Nsym
           if(is_diagonal(H_basis(isym)%O))then
              dmft_bath%item(ibath)%lambda(isym)=rescale(ibath)*lambda_impHloc(isym)
           else
              dmft_bath%item(ibath)%lambda(isym) =  lambda_impHloc(isym)
           endif
        enddo
     enddo
     !
  end select
  !
  !
  !
  !Read from file if exist:
  !
  inquire(file=trim(Hfile)//trim(ed_file_suffix)//".restart",exist=IOfile)
  if(IOfile)then
     write(LOGfile,"(A)")'Reading bath from file'//trim(Hfile)//trim(ed_file_suffix)//".restart"
     unit = free_unit()
     flen = file_length(trim(Hfile)//trim(ed_file_suffix)//".restart")
     !
     open(unit,file=trim(Hfile)//trim(ed_file_suffix)//".restart")
     !
     select case(bath_type)
     case default
        !
        read(unit,*)
        select case(ed_mode)
        case default
           do i=1,min(flen,Nbath)
              read(unit,*)((&
                   dmft_bath_%e(ispin,iorb,i),&
                   dmft_bath_%v(ispin,iorb,i),&
                   iorb=1,Norb),ispin=1,Nspin)
           enddo
        case ("superc")
           do i=1,min(flen,Nbath)
              read(unit,*)((&
                   dmft_bath_%e(ispin,iorb,i),&
                   dmft_bath_%d(ispin,iorb,i),&
                   dmft_bath_%v(ispin,iorb,i),&
                   iorb=1,Norb),ispin=1,Nspin)
           enddo
        case("nonsu2")
           do i=1,min(flen,Nbath)
              read(unit,*)((&
                   dmft_bath_%e(ispin,iorb,i),&
                   dmft_bath_%v(ispin,iorb,i),&
                   dmft_bath_%u(ispin,iorb,i),&
                   iorb=1,Norb),ispin=1,Nspin)
           enddo
        end select
        !
     case ('hybrid')
        read(unit,*)
        !
        select case(ed_mode)
        case default
           do i=1,min(flen,Nbath)
              read(unit,*)(&
                   dmft_bath_%e(ispin,1,i),&
                   (&
                   dmft_bath_%v(ispin,iorb,i),&
                   iorb=1,Norb),&
                   ispin=1,Nspin)
           enddo
        case ("superc")
           do i=1,min(flen,Nbath)
              read(unit,*)(&
                   dmft_bath_%e(ispin,1,i),&
                   dmft_bath_%d(ispin,1,i),&
                   (&
                   dmft_bath_%v(ispin,iorb,i),&
                   iorb=1,Norb),&
                   ispin=1,Nspin)
           enddo
        case ("nonsu2")
           do i=1,min(flen,Nbath)
              read(unit,*)(&
                   dmft_bath_%e(ispin,1,i),&
                   (&
                   dmft_bath_%v(ispin,iorb,i),&
                   dmft_bath_%u(ispin,iorb,i),&
                   iorb=1,Norb),&
                   ispin=1,Nspin)
           enddo
        end select
        !
     case ('replica')
        !
        !read number of lambdas
        do ibath=1,Nbath
           read(unit,"(I3)")Nlambdas(ibath)
        enddo
        do ibath=1,Nbath
           !read V
           do ispin=1,Nspin
              read(unit,*)dmft_bath%item(ibath)%v(ispin)
           enddo
           !read lambdas
           read(unit,*)(dmft_bath%item(ibath)%lambda(jo),jo=1,Nlambdas(ibath))
        enddo
        !
        !
     end select
     close(unit)
  endif
end subroutine init_dmft_bath



!+-------------------------------------------------------------------+
!PURPOSE  : write out the bath to a given unit with 
! the following column formatting: 
! [(Ek_iorb,Vk_iorb)_iorb=1,Norb]_ispin=1,Nspin
!+-------------------------------------------------------------------+
subroutine write_dmft_bath(dmft_bath_,unit)
  type(effective_bath) :: dmft_bath_
  integer,optional     :: unit
  integer              :: unit_
  integer              :: i
  integer              :: io,jo,iorb,ispin,isym
  real(8)              :: hybr_aux
  real(8)              :: hrep_aux(Nspin*Norb,Nspin*Norb)
  !
  character(len=64)    :: string_fmt,string_fmt_first
  unit_=LOGfile;if(present(unit))unit_=unit
  if(.not.dmft_bath_%status)stop "write_dmft_bath error: bath not allocated"
  select case(bath_type)
  case default
     !
     select case(ed_mode)
     case default
        write(unit_,"(90(A21,1X))")&
             ((&
             "#Ek_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
             "Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
             iorb=1,Norb),ispin=1,Nspin)
        do i=1,Nbath
           write(unit_,"(90(F21.12,1X))")((&
                dmft_bath_%e(ispin,iorb,i),&
                dmft_bath_%v(ispin,iorb,i),&
                iorb=1,Norb),ispin=1,Nspin)
        enddo
     case ("superc")
        write(unit_,"(90(A21,1X))")&
             ((&
             "#Ek_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
             "Dk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)) ,&
             "Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
             iorb=1,Norb),ispin=1,Nspin)
        do i=1,Nbath
           write(unit_,"(90(F21.12,1X))")((&
                dmft_bath_%e(ispin,iorb,i),&
                dmft_bath_%d(ispin,iorb,i),&
                dmft_bath_%v(ispin,iorb,i),&
                iorb=1,Norb),ispin=1,Nspin)
        enddo
     case ("nonsu2")
        write(unit_,"(90(A21,1X))")&
             ((&
             "#Ek_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
             "Vak_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
             "Vbk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
             iorb=1,Norb), ispin=1,Nspin)
        do i=1,Nbath
           write(unit,"(90(F21.12,1X))")((&
                dmft_bath_%e(ispin,iorb,i),&
                dmft_bath_%v(ispin,iorb,i),&
                dmft_bath_%u(ispin,iorb,i),&
                iorb=1,Norb),ispin=1,Nspin)
        enddo
     end select
     !
  case('hybrid')
     !
     select case(ed_mode)
     case default
        write(unit_,"(90(A21,1X))")(&
             "#Ek_s"//reg(txtfy(ispin)),&
             ("Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),iorb=1,Norb),&
             ispin=1,Nspin)
        do i=1,Nbath
           write(unit_,"(90(F21.12,1X))")(&
                dmft_bath_%e(ispin,1,i),&
                (dmft_bath_%v(ispin,iorb,i),iorb=1,Norb),&
                ispin=1,Nspin)
        enddo
     case ("superc")
        write(unit_,"(90(A21,1X))")(&
             "#Ek_s"//reg(txtfy(ispin)),&
             "Dk_s"//reg(txtfy(ispin)) ,&
             ("Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),iorb=1,Norb),&
             ispin=1,Nspin)
        do i=1,Nbath
           write(unit_,"(90(F21.12,1X))")(&
                dmft_bath_%e(ispin,1,i),&
                dmft_bath_%d(ispin,1,i),&
                (dmft_bath_%v(ispin,iorb,i),iorb=1,Norb),&
                ispin=1,Nspin)
        enddo
     case ("nonsu2")
        write(unit_,"(90(A21,1X))")(&
             "#Ek_s"//reg(txtfy(ispin)),&
             (&
             "Vak_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
             "Vbk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
             iorb=1,Norb),&
             ispin=1,Nspin)
        do i=1,Nbath
           write(unit_,"(90(F21.12,1X))")(&
                dmft_bath_%e(ispin,1,i),    &
                (dmft_bath_%v(ispin,iorb,i),dmft_bath_%u(ispin,iorb,i),iorb=1,Norb),&
                ispin=1,Nspin)
        enddo
     end select
     !
  case ('replica')
     !
     string_fmt      ="(A8,A5,"//str(Nspin*Norb)//"(ES8.4,1X))"
     !
     if(unit_==LOGfile)then
        if(Nspin*Norb.le.8)then
           write(unit_,"(A1)")" "
           write(unit_,"(A8,A3,a5,90(A15,1X))")"V","||"," ","Re(H) | Im(H)"
           do ibath=1,Nbath
              write(unit_,"(A1)")" "
              Hrep_aux   = nn2so_reshape( bath_from_sym(dmft_bath_%item(ibath)%lambda) ,Nspin,Norb)
              do ispin=1,Nspin
                 write(unit_,"(ES8.4)")dmft_bath_%item(ibath)%v(ispin)
              enddo
              do io=1,Nspin*Norb
                 write(unit_,string_fmt) "  "  ,"||  ",(hrep_aux(io,jo),jo=1,Nspin*Norb)
              enddo
           enddo
           write(unit_,"(A1)")" "
        else
           write(LOGfile,"(A)")"Bath matrix too large to print: printing the parameters (including eventual offset)."
           write(unit_,"(A8,A5,90(A8,1X))")"V"," ","lambdas"        
           do ibath=1,Nbath
              do ispin=1,Nspin
                 write(unit_,"(ES8.4)")dmft_bath_%item(ibath)%v(ispin)
              enddo
              write(unit_,"(A8,A5,90(ES8.4,1X))")"","|   ",&
                   (dmft_bath_%item(ibath)%lambda(io),io=1,dmft_bath_%item(ibath)%N_dec)
           enddo
        endif
     else
        do ibath=1,Nbath
           write(unit,"(I3)")dmft_bath_%item(ibath)%N_dec
        enddo
        do ibath=1,Nbath
           do ispin=1,Nspin
              write(unit,*)dmft_bath_%item(ibath)%v(ispin)
           enddo
           write(unit,*)(dmft_bath_%item(ibath)%lambda(jo),jo=1,dmft_bath_%item(ibath)%N_dec)
        enddo
     endif
     !
  end select
end subroutine write_dmft_bath






!+-------------------------------------------------------------------+
!PURPOSE  : save the bath to a given file using the write bath
! procedure and formatting: 
!+-------------------------------------------------------------------+
subroutine save_dmft_bath(dmft_bath_,file,used)
  type(effective_bath)      :: dmft_bath_
  character(len=*),optional :: file
  character(len=256)        :: file_
  logical,optional          :: used
  logical                   :: used_
  character(len=16)         :: extension
  integer                   :: unit_
  if(.not.dmft_bath_%status)stop "save_dmft_bath error: bath is not allocated"
  used_=.false.;if(present(used))used_=used
  extension=".restart";if(used_)extension=".used"
  file_=str(str(Hfile)//str(ed_file_suffix)//str(extension))
  if(present(file))file_=str(file)
  unit_=free_unit()
  open(unit_,file=str(file_))
  call write_dmft_bath(dmft_bath_,unit_)
  close(unit_)
end subroutine save_dmft_bath




!+-------------------------------------------------------------------+
!PURPOSE  : set the bath components from a given user provided 
! bath-array 
!+-------------------------------------------------------------------+
subroutine set_dmft_bath(bath_,dmft_bath_)
  real(8),dimension(:)   :: bath_
  type(effective_bath)   :: dmft_bath_
  integer                :: stride,io,jo,i
  integer                :: iorb,ispin,jorb,jspin,ibath
  logical                :: check
  !
  if(.not.dmft_bath_%status)stop "set_dmft_bath error: bath not allocated"
  check = check_bath_dimension(bath_)
  if(.not.check)stop "set_dmft_bath error: wrong bath dimensions"
  !
  select case(bath_type)
  case default
     !
     select case(ed_mode)
        !
     case default
        stride = 0
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 dmft_bath_%e(ispin,iorb,i) = bath_(io)
              enddo
           enddo
        enddo
        stride = Nspin*Norb*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 dmft_bath_%v(ispin,iorb,i) = bath_(io)
              enddo
           enddo
        enddo
        !
     case ("superc")
        stride = 0
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 dmft_bath_%e(ispin,iorb,i) = bath_(io)
              enddo
           enddo
        enddo
        stride = Nspin*Norb*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 dmft_bath_%d(ispin,iorb,i) = bath_(io)
              enddo
           enddo
        enddo
        stride = 2*Nspin*Norb*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 dmft_bath_%v(ispin,iorb,i) = bath_(io)
              enddo
           enddo
        enddo
        !
     case("nonsu2")
        stride = 0
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 dmft_bath_%e(ispin,iorb,i) = bath_(io)
              enddo
           enddo
        enddo
        stride = Nspin*Norb*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 dmft_bath_%v(ispin,iorb,i) = bath_(io)
              enddo
           enddo
        enddo
        stride = 2*Nspin*Norb*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 dmft_bath_%u(ispin,iorb,i) = bath_(io)
              enddo
           enddo
        enddo
     end select
     !
     !
  case ('hybrid')
     !
     select case(ed_mode)
     case default
        stride = 0
        do ispin=1,Nspin
           do i=1,Nbath
              io = stride + i + (ispin-1)*Nbath
              dmft_bath_%e(ispin,1,i) = bath_(io)
           enddo
        enddo
        stride = Nspin*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                 dmft_bath_%v(ispin,iorb,i) = bath_(io)
              enddo
           enddo
        enddo
        !
     case ("superc")
        stride = 0
        do ispin=1,Nspin
           do i=1,Nbath
              io = stride + i + (ispin-1)*Nbath
              dmft_bath_%e(ispin,1,i) = bath_(io)
           enddo
        enddo
        stride = Nspin*Nbath
        do ispin=1,Nspin
           do i=1,Nbath
              io = stride + i + (ispin-1)*Nbath
              dmft_bath_%d(ispin,1,i) = bath_(io)
           enddo
        enddo
        stride = 2*Nspin*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                 dmft_bath_%v(ispin,iorb,i) = bath_(io)
              enddo
           enddo
        enddo
        !
     case("nonsu2")
        stride = 0
        do ispin=1,Nspin
           do i=1,Nbath
              io = stride + i + (ispin-1)*Nbath
              dmft_bath_%e(ispin,1,i) = bath_(io)
           enddo
        enddo
        stride = Nspin*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                 dmft_bath_%v(ispin,iorb,i) = bath_(io)
              enddo
           enddo
        enddo
        stride = Nspin*Nbath + Nspin*Norb*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                 dmft_bath_%u(ispin,iorb,i) = bath_(io)
              enddo
           enddo
        enddo
     end select
     !
     !
  case ('replica')
     !
     stride = 0
     !Get N_dec
     do ibath=1,Nbath
        stride = stride + 1
        dmft_bath_%item(ibath)%N_dec=NINT(bath_(stride))
     enddo
     !Get N_dec
     !get Lambdas
     do ibath=1,Nbath
        do ispin=1,Nspin
           stride = stride + 1
           dmft_bath_%item(ibath)%v(ispin) = bath_(stride)
        enddo
        dmft_bath_%item(ibath)%lambda=bath_(stride+1 :stride+dmft_bath_%item(ibath)%N_dec)
        stride=stride+dmft_bath_%item(ibath)%N_dec
     enddo
  end select
end subroutine set_dmft_bath



!+-------------------------------------------------------------------+
!PURPOSE  : copy the bath components back to a 1-dim array 
!+-------------------------------------------------------------------+
subroutine get_dmft_bath(dmft_bath_,bath_)
  type(effective_bath)   :: dmft_bath_
  real(8),dimension(:)   :: bath_
  real(8)                :: hrep_aux(Nspin*Norb,Nspin*Norb)
  integer                :: stride,io,jo,i
  integer                :: iorb,ispin,jorb,jspin,ibath,maxspin
  logical                :: check
  if(.not.dmft_bath_%status)stop "get_dmft_bath error: bath not allocated"
  check=check_bath_dimension(bath_)
  if(.not.check)stop "get_dmft_bath error: wrong bath dimensions"
  !
  select case(bath_type)
  case default
     !
     select case(ed_mode)
     case default
        stride = 0
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 bath_(io) = dmft_bath_%e(ispin,iorb,i) 
              enddo
           enddo
        enddo
        stride = Nspin*Norb*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 bath_(io) = dmft_bath_%v(ispin,iorb,i)
              enddo
           enddo
        enddo
        !
     case ("superc")
        stride = 0
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 bath_(io) = dmft_bath_%e(ispin,iorb,i)
              enddo
           enddo
        enddo
        stride = Nspin*Norb*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 bath_(io) =  dmft_bath_%d(ispin,iorb,i)
              enddo
           enddo
        enddo
        stride = 2*Nspin*Norb*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 bath_(io) =  dmft_bath_%v(ispin,iorb,i)
              enddo
           enddo
        enddo
        !
     case("nonsu2")
        stride = 0
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 bath_(io) = dmft_bath_%e(ispin,iorb,i)
              enddo
           enddo
        enddo
        stride = Nspin*Norb*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 bath_(io) = dmft_bath_%v(ispin,iorb,i)
              enddo
           enddo
        enddo
        stride = 2*Nspin*Norb*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 bath_(io) = dmft_bath_%u(ispin,iorb,i)
              enddo
           enddo
        enddo
     end select
     !
  case ('hybrid')
     !
     select case(ed_mode)
     case default
        stride = 0
        do ispin=1,Nspin
           do i=1,Nbath
              io = stride + i + (ispin-1)*Nbath
              bath_(io) =  dmft_bath_%e(ispin,1,i)
           enddo
        enddo
        stride = Nspin*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                 bath_(io) =  dmft_bath_%v(ispin,iorb,i)
              enddo
           enddo
        enddo
        !
     case ("superc")
        stride = 0
        do ispin=1,Nspin
           do i=1,Nbath
              io = stride + i + (ispin-1)*Nbath
              bath_(io) =  dmft_bath_%e(ispin,1,i)
           enddo
        enddo
        stride = Nspin*Nbath
        do ispin=1,Nspin
           do i=1,Nbath
              io = stride + i + (ispin-1)*Nbath
              bath_(io) =  dmft_bath_%d(ispin,1,i)
           enddo
        enddo
        stride = 2*Nspin*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                 bath_(io) =  dmft_bath_%v(ispin,iorb,i)
              enddo
           enddo
        enddo
        !
     case("nonsu2")
        stride = 0
        do ispin=1,Nspin
           do i=1,Nbath
              io = stride + i + (ispin-1)*Nbath
              bath_(io) = dmft_bath_%e(ispin,1,i) 
           enddo
        enddo
        stride = Nspin*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                 bath_(io) = dmft_bath_%v(ispin,iorb,i)
              enddo
           enddo
        enddo
        stride = Nspin*Nbath + Nspin*Norb*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                 bath_(io) = dmft_bath_%u(ispin,iorb,i)
              enddo
           enddo
        enddo
     end select
     !
     !
  case ('replica')
     !
     stride = 0
     do ibath=1,Nbath
        stride = stride + 1
        bath_(stride)=dmft_bath_%item(ibath)%N_dec
     enddo
     do ibath=1,Nbath
        do ispin=1,Nspin
           stride = stride + 1
           bath_(stride)=dmft_bath_%item(ibath)%v(ispin)
        enddo
        bath_(stride+1 : stride+dmft_bath_%item(ibath)%N_dec)=dmft_bath_%item(ibath)%lambda
        stride=stride+dmft_bath_%item(ibath)%N_dec
     enddo
  end select
end subroutine get_dmft_bath





!##################################################################
!
!     W_hyb PROCEDURES
!
!##################################################################
!+-----------------------------------------------------------------------------+!
!PURPOSE: build up the all-spin hybridization matrix W_{ss`}
!+-----------------------------------------------------------------------------+!
function get_Whyb_matrix_1orb(v,u) result(w)
  real(8),dimension(Nspin,Nbath)       :: v,u
  real(8),dimension(Nspin,Nspin,Nbath) :: w
  integer                              :: ispin
  !
  ! if(ed_para)then
  !    do ispin=1,Nspin
  !       w(ispin,ispin,:) = v(1,:)
  !    enddo
  !    w(1,Nspin,:) = u(1,:)
  !    w(Nspin,1,:) = u(1,:)
  ! else
  do ispin=1,Nspin
     w(ispin,ispin,:) = v(ispin,:)
  enddo
  w(1,Nspin,:) = u(1,:)
  w(Nspin,1,:) = u(2,:)
  ! endif
  !
end function get_Whyb_matrix_1orb

function get_Whyb_matrix_Aorb(v,u) result(w)
  real(8),dimension(Nspin,Norb,Nbath)       :: v,u
  real(8),dimension(Nspin,Nspin,Norb,Nbath) :: w
  integer                                   :: ispin
  !
  ! if(ed_para)then
  !    do ispin=1,Nspin
  !       w(ispin,ispin,:,:) = v(1,:,:)
  !    enddo
  !    w(1,Nspin,:,:) = u(1,:,:)
  !    w(Nspin,1,:,:) = u(1,:,:)
  ! else
  do ispin=1,Nspin
     w(ispin,ispin,:,:) = v(ispin,:,:)
  enddo
  w(1,Nspin,:,:) = u(1,:,:)
  w(Nspin,1,:,:) = u(2,:,:)
  ! endif
  !
end function get_Whyb_matrix_Aorb

function get_Whyb_matrix_dmft_bath(dmft_bath_) result(w)
  type(effective_bath)                      :: dmft_bath_
  real(8),dimension(Nspin,Nspin,Norb,Nbath) :: w
  integer                                   :: ispin
  !
  ! if(ed_para)then
  !    do ispin=1,Nspin
  !       w(ispin,ispin,:,:) = dmft_bath_%v(1,:,:)
  !    enddo
  !    w(1,Nspin,:,:) = dmft_bath_%u(1,:,:)
  !    w(Nspin,1,:,:) = dmft_bath_%u(1,:,:)
  ! else
  do ispin=1,Nspin
     w(ispin,ispin,:,:) = dmft_bath_%v(ispin,:,:)
  enddo
  w(1,Nspin,:,:) = dmft_bath_%u(1,:,:)
  w(Nspin,1,:,:) = dmft_bath_%u(2,:,:)
  ! endif
  !
end function get_Whyb_matrix_dmft_bath



