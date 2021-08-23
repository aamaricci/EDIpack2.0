!##################################################################
! THE CALCULATION OF THE \chi^2 FUNCTIONS USE PROCEDURES FURTHER 
! BELOW TO EVALUATE INDEPENDENTLY THE ANDERSON MODEL:
!  - DELTA, 
!  -\GRAD DELTA
!  - G0
! THE LATTER ARE ADAPTED FROM THE PROCEDURES:
! DELTA_BATH_MATS
! GRAD_DELTA_BATH_MATS
! G0 BATH_MATS
! FOR, YOU NEED TO DECOMPOSE THE a INPUT ARRAY INTO ELEMENTS.
!##################################################################


!+-------------------------------------------------------------+
!PURPOSE  : Chi^2 interface for Hybridized bath nonsu2 phase
!+-------------------------------------------------------------+
subroutine chi2_fitgf_hybrid_nonsu2(fg,bath_)
  complex(8),dimension(:,:,:,:,:)    :: fg ![Nspin][Nspin][Norb][Norb][Lmats]
  real(8),dimension(:),intent(inout) :: bath_
  real(8),dimension(:),allocatable   :: array_bath
  integer                            :: i,j,iorb,jorb,ispin,jspin,io,jo
  integer                            :: iter,stride,count,Asize
  real(8)                            :: chi
  logical                            :: check
  type(effective_bath)               :: dmft_bath
  character(len=256)                 :: suffix
  integer                            :: unit
  complex(8),dimension(:,:,:,:,:),allocatable :: fgand ![Nspin][][Norb][][Ldelta]
  !
#ifdef _DEBUG
  if(ed_verbose>2)write(Logfile,"(A)")"DEBUG chi2_fitgf_hybrid_nonsu2: Fit"
#endif
  !
  if(size(fg,1)/=Nspin)stop "chi2_fitgf_hybrid_nonsu2 error: size[fg,1]!=Nspin"
  if(size(fg,2)/=Nspin)stop "chi2_fitgf_hybrid_nonsu2 error: size[fg,2]!=Nspin"
  if(size(fg,3)/=Norb)stop "chi2_fitgf_hybrid_nonsu2 error: size[fg,3]!=Norb"
  if(size(fg,4)/=Norb)stop "chi2_fitgf_hybrid_nonsu2 error: size[fg,4]!=Norb"
  !
  check= check_bath_dimension(bath_)
  if(.not.check)stop "chi2_fitgf_hybrid_nonsu2 error: wrong bath dimensions"
  !
  Ldelta = Lfit ; if(Ldelta>size(fg,5))Ldelta=size(fg,5)
  !
  totNso = Nspin*(Nspin+1)/2 * Norb*(Norb+1)/2
  allocate(getIspin(totNso),getJspin(totNso))
  allocate(getIorb(totNso),getJorb(totNso))
  count=0
  do ispin=1,Nspin
     do jspin=ispin,Nspin
        do iorb=1,Norb
           do jorb=iorb,Norb
              count=count+1
              getIspin(count) = ispin
              getIorb(count)  = iorb
              getJspin(count) = jspin
              getJorb(count)  = jorb
           enddo
        enddo
     enddo
  enddo
  if(totNso/=count)stop "chi2_fitgf_hybrid_nonsu2: error counting the spin-orbitals"
  !
  allocate(Gdelta(totNso,Ldelta))
  allocate(Xdelta(Ldelta))
  allocate(Wdelta(Ldelta))
  !
  Xdelta = pi/beta*(2*arange(1,Ldelta)-1)
  !
  select case(cg_weight)
  case default
     Wdelta=1d0
  case(2)
     Wdelta=1d0*arange(1,Ldelta)
  case(3)
     Wdelta=Xdelta
  end select
  !
  call allocate_dmft_bath(dmft_bath)
  call set_dmft_bath(bath_,dmft_bath)
  !
  do i=1,totNso
     Gdelta(i,1:Ldelta) = fg(getIspin(i),getJspin(i),getIorb(i),getJorb(i),1:Ldelta)
  enddo
  !
  Asize = Nbath + Norb*Nbath + Norb*Nbath
  Asize = Nspin*Asize
  ! if(.not.para_)Asize=Nspin*Asize !fit all spin components
  allocate(array_bath(Asize))
  !
  ! size = Nspin*Nbath + Nspin*Norb*Nbath + Nspin*Norb*Nbath
  stride = 0
  do ispin=1,Nspin
     do i=1,Nbath
        io = stride + i + (ispin-1)*Nbath
        array_bath(io) = dmft_bath%e(ispin,1,i)
     enddo
  enddo
  stride = Nspin*Nbath
  do ispin=1,Nspin
     do iorb=1,Norb
        do i=1,Nbath
           io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
           array_bath(io) = dmft_bath%v(ispin,iorb,i)
        enddo
     enddo
  enddo
  stride =  Nspin*Nbath + Nspin*Norb*Nbath
  do ispin=1,Nspin
     do iorb=1,Norb
        do i=1,Nbath
           io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
           array_bath(io) = dmft_bath%u(ispin,iorb,i)
        enddo
     enddo
  enddo
  !
#ifdef _DEBUG
  if(ed_verbose>3)write(Logfile,"(A)")&
       "DEBUG chi2_fitgf_normal_nonsu2: cg_method:"//str(cg_method)//&
       ", cg_grad:"//str(cg_grad)//&
       ", cg_scheme:"//str(cg_scheme)
#endif
  !
  select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE; 2=CG+
  case default
     select case (cg_scheme)
     case ("weiss")
        call fmin_cg(array_bath,chi2_weiss_hybrid_nonsu2,&
             iter,chi,&
             itmax=cg_niter,&
             ftol=cg_Ftol,  &
             istop=cg_stop, &
             iverbose=(ed_verbose>3))
     case ("delta")
        call fmin_cg(array_bath,chi2_delta_hybrid_nonsu2,&
             iter,chi,&
             itmax=cg_niter,&
             ftol=cg_Ftol,  &
             istop=cg_stop, &
             iverbose=(ed_verbose>3))
     case default
        stop "chi2_fitgf_hybrid_nonsu2 error: cg_scheme != [weiss,delta]"
     end select
     !
  case (1)
     select case (cg_scheme)
     case ("weiss")
        call fmin_cgminimize(array_bath,chi2_weiss_hybrid_nonsu2,&
             iter,chi,itmax=cg_niter,ftol=cg_Ftol,&
             new_version=cg_minimize_ver,&
             hh_par=cg_minimize_hh,&
             iverbose=(ed_verbose>3))
     case ("delta")
        call fmin_cgminimize(array_bath,chi2_delta_hybrid_nonsu2,&
             iter,chi,itmax=cg_niter,ftol=cg_Ftol,&
             new_version=cg_minimize_ver,&
             hh_par=cg_minimize_hh,&
             iverbose=(ed_verbose>3))
     case default
        stop "chi2_fitgf_hybrid_nonsu2 error: cg_scheme != [weiss,delta]"
     end select
     !
  end select
  !
  !
  write(LOGfile,"(A,ES18.9,A,I5,A)")&
       "chi^2|iter"//reg(ed_file_suffix)//'= ',chi," | ",iter,&
       "  <--  All Orbs, All Spins"
  !
  suffix="_ALLorb_ALLspins"//reg(ed_file_suffix)
  unit=free_unit()
  open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
  write(unit,"(ES18.9,1x,I5)") chi,iter
  close(unit)
  !
  ! size = Nspin*Nbath + Nspin*Norb*Nbath + Nspin*Norb*Nbath
  stride = 0
  do ispin=1,Nspin
     do i=1,Nbath
        io = stride + i + (ispin-1)*Nbath
        dmft_bath%e(ispin,1,i) = array_bath(io)
     enddo
  enddo
  stride = Nspin*Nbath
  do ispin=1,Nspin
     do iorb=1,Norb
        do i=1,Nbath
           io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
           dmft_bath%v(ispin,iorb,i) = array_bath(io)
        enddo
     enddo
  enddo
  stride =  Nspin*Nbath + Nspin*Norb*Nbath
  do ispin=1,Nspin
     do iorb=1,Norb
        do i=1,Nbath
           io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
           dmft_bath%u(ispin,iorb,i) = array_bath(io)
        enddo
     enddo
  enddo
  ! dmft_bath%e(Nspin,1,:) = dmft_bath%e(1,1,:)
  ! dmft_bath%v(Nspin,:,:) = dmft_bath%v(1,:,:)
  ! dmft_bath%u(Nspin,:,:) = dmft_bath%u(1,:,:)
  !
  call write_dmft_bath(dmft_bath,LOGfile)
  !
  call save_dmft_bath(dmft_bath)
  !
  allocate(fgand(Nspin,Nspin,Norb,Norb,Ldelta))
  if(cg_scheme=='weiss')then
     fgand = g0and_bath_function(xi*Xdelta(:),dmft_bath)
  else
     fgand = delta_bath_function(xi*Xdelta(:),dmft_bath)
  endif
  call write_fit_result()
  deallocate(fgand)
  !
  call get_dmft_bath(dmft_bath,bath_)
  call deallocate_dmft_bath(dmft_bath)
  deallocate(Gdelta,Xdelta,Wdelta)
  deallocate(getIspin,getJspin)
  deallocate(getIorb,getJorb)
  !
contains
  !
  subroutine write_fit_result()
    integer           :: i,j,s,l,iorb,jorb,ispin,jspin
    !
    do l=1,totNso
       iorb = getIorb(l)
       jorb = getJorb(l)
       ispin = getIspin(l)
       jspin = getJspin(l)
       suffix="_l"//reg(txtfy(iorb))//&
            "_m"//reg(txtfy(jorb))//&
            "_s"//reg(txtfy(ispin))//&
            "_r"//reg(txtfy(jspin))//reg(ed_file_suffix)
       if(cg_scheme=='weiss')then
          open(free_unit(unit),file="fit_weiss"//reg(suffix)//".ed")
       else
          open(free_unit(unit),file="fit_delta"//reg(suffix)//".ed")
       endif
       do i=1,Ldelta
          write(unit,"(5F24.15)")Xdelta(i),&
               dimag(fg(ispin,jspin,iorb,jorb,i)),dimag(fgand(ispin,jspin,iorb,jorb,i)),&
               dreal(fg(ispin,jspin,iorb,jorb,i)),dreal(fgand(ispin,jspin,iorb,jorb,i))
       enddo
       close(unit)
    enddo
  end subroutine write_fit_result
  !
end subroutine chi2_fitgf_hybrid_nonsu2









!##################################################################
! THESE PROCEDURES EVALUATES THE \chi^2 FUNCTIONS TO MINIMIZE. 
!##################################################################
!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function.
!+-------------------------------------------------------------+
function chi2_delta_hybrid_nonsu2(a) result(chi2)
  real(8),dimension(:)                               ::  a
  real(8)                                            ::  chi2
  real(8),dimension(totNso)                          ::  chi2_so
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta) ::  Delta
  real(8),dimension(Ldelta)                          ::  Ctmp
  integer                                            ::  i,l,iorb,jorb,ispin,jspin
  !
#ifdef _DEBUG
  if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_delta_hybrid_nonsu2. a:",a
#endif
  !
  Delta = delta_hybrid_nonsu2(a)
  !
  do l=1,totNso
     iorb = getIorb(l)
     jorb = getJorb(l)
     ispin = getIspin(l)
     jspin = getJspin(l)
     Ctmp=abs(Gdelta(l,:)-Delta(ispin,jspin,iorb,jorb,:))
     chi2_so(l) = sum( Ctmp**cg_pow/Wdelta )
  enddo
  !
  chi2=sum(chi2_so)/Ldelta
#ifdef _DEBUG
  if(ed_verbose>3)write(Logfile,"(A,ES10.2)")"DEBUG chi2_delta_hybrid_nonsu2. Chi**2:",chi2
#endif
  !
end function chi2_delta_hybrid_nonsu2


!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance of G_0_Anderson function 
! The Gradient is not evaluated, so the minimization requires 
! a numerical estimate of the gradient. 
!+-------------------------------------------------------------+
function chi2_weiss_hybrid_nonsu2(a) result(chi2)
  real(8),dimension(:)                               :: a
  real(8)                                            :: chi2
  real(8),dimension(totNso)                          :: chi2_so
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta) :: g0and
  real(8),dimension(Ldelta)                          ::  Ctmp
  integer                                            :: i,l,iorb,jorb,ispin,jspin
  !
#ifdef _DEBUG
  if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG chi2_weiss_hybrid_nonsu2. a:",a
#endif
  !
  g0and = g0and_hybrid_nonsu2(a)
  !
  do l=1,totNso
     iorb = getIorb(l)
     jorb = getJorb(l)
     ispin = getIspin(l)
     jspin = getJspin(l)
     Ctmp=abs(Gdelta(l,:)-g0and(ispin,jspin,iorb,jorb,:))
     chi2_so(l) = sum( Ctmp**cg_pow/Wdelta )
  enddo
  !
  chi2=sum(chi2_so)/Ldelta
#ifdef _DEBUG
  if(ed_verbose>3)write(Logfile,"(A,ES10.2)")"DEBUG chi2_weiss_hybrid_nonsu2. Chi**2:",chi2
#endif
  !
  !
end function chi2_weiss_hybrid_nonsu2






!##################################################################
! THESE PROCEDURES EVALUATES THE 
! - \delta
! - g0
! FUNCTIONS. 
!##################################################################
function delta_hybrid_nonsu2(a) result(Delta)
  real(8),dimension(:)                               :: a
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta) :: Delta
  integer                                            :: ispin,jspin,iorb,jorb
  integer                                            :: i,io,stride,ih
  real(8),dimension(Nspin,Nbath)                     :: hps
  real(8),dimension(Nspin,Norb,Nbath)                :: vops
  real(8),dimension(Nspin,Norb,Nbath)                :: uops
  real(8),dimension(Nhel,Nbath)                      :: ehel
  real(8),dimension(Nhel,Nhel,Norb,Nbath)            :: wohel

  stride = 0
  do ispin=1,Nspin
     do i=1,Nbath
        io = stride + i + (ispin-1)*Nbath
        hps(ispin,i) = a(io)
     enddo
  enddo
  stride = Nspin*Nbath
  do ispin=1,Nspin
     do iorb=1,Norb
        do i=1,Nbath
           io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
           vops(ispin,iorb,i) = a(io)
        enddo
     enddo
  enddo
  stride =  Nspin*Nbath + Nspin*Norb*Nbath
  do ispin=1,Nspin
     do iorb=1,Norb
        do i=1,Nbath
           io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
           uops(ispin,iorb,i) = a(io)
        enddo
     enddo
  enddo
  ! hps(Nspin,:)    = hps(ispin,:)
  ! vops(Nspin,:,:) = vops(ispin,:,:)
  ! uops(Nspin,:,:) = uops(ispin,:,:)
  !
  wohel = get_Whyb_matrix(vops(1:Nspin,1:Norb,1:Nbath),uops(1:Nspin,1:Norb,1:Nbath))
  !
  Delta=zero
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              do i=1,Ldelta
                 do ih=1,Nspin
                    Delta(ispin,jspin,iorb,jorb,i) = Delta(ispin,jspin,iorb,jorb,i) + &
                         sum( wohel(ispin,ih,iorb,:)*wohel(jspin,ih,jorb,:)/(xi*Xdelta(i) - hps(ih,:)) )
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  !
end function delta_hybrid_nonsu2

function g0and_hybrid_nonsu2(a) result(G0and)
  real(8),dimension(:)                               :: a
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta) :: G0and,Delta
  complex(8),dimension(Norb*Nspin,Norb*Nspin)        :: zeta,fgorb
  integer                                            :: i,Nso
  integer                                            :: iorb,jorb,ispin,jspin,io,jo
  !
  Nso = Norb*Nspin
  !
  Delta(:,:,:,:,:) = delta_hybrid_nonsu2(a)
  !
  do i=1,Ldelta
     fgorb=zero
     zeta = (xi*Xdelta(i)+xmu)*zeye(Nso)
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 io = iorb + (ispin-1)*Norb
                 jo = jorb + (jspin-1)*Norb
                 fgorb(io,jo)   = zeta(io,jo) - impHloc(ispin,jspin,iorb,jorb) - Delta(ispin,jspin,iorb,jorb,i)
              enddo
           enddo
        enddo
     enddo
     call inv(fgorb)
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 io = iorb + (ispin-1)*Norb
                 jo = jorb + (jspin-1)*Norb
                 G0and(ispin,jspin,iorb,jorb,i) = fgorb(io,jo)
              enddo
           enddo
        enddo
     enddo
  enddo
  !
end function g0and_hybrid_nonsu2
