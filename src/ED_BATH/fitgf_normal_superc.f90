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
!PURPOSE  : Chi^2 interface for Irreducible bath Superconducting phase
!+-------------------------------------------------------------+
subroutine chi2_fitgf_normal_superc(fg,bath_,ispin,iorb)
  complex(8),dimension(:,:,:,:)      :: fg ![2][Norb][Norb][Lmats]
  real(8),dimension(:),intent(inout) :: bath_
  integer                            :: ispin
  integer,optional                   :: iorb
  real(8),dimension(:),allocatable   :: array_bath
  integer                            :: iter,stride,i,io,j,jorb,Asize
  real(8)                            :: chi
  logical                            :: check
  type(effective_bath)               :: dmft_bath
  character(len=256)                 :: suffix
  integer                            :: unit
  complex(8),dimension(:,:,:,:,:),allocatable :: fgand,ffand ![Nspin][][Norb][][Ldelta]
  !
  if(size(fg,1)/=2)stop "chi2_fitgf_normal_superc error: size[fg,1]!=2"
  if(size(fg,2)/=Norb)stop "chi2_fitgf_normal_superc error: size[fg,2]!=Norb"
  if(size(fg,3)/=Norb)stop "chi2_fitgf_normal_superc error: size[fg,3]!=Norb"
  !
  check= check_bath_dimension(bath_)
  if(.not.check)stop "chi2_fitgf_normal_superc: wrong bath dimensions"
  !
  Ldelta = Lfit ; if(Ldelta>size(fg,4))Ldelta=size(fg,4)
  !
  allocate(Gdelta(1,Ldelta))
  allocate(Fdelta(1,Ldelta))
  allocate(Xdelta(Ldelta))
  allocate(Wdelta(Ldelta))
  !
  Xdelta = pi/beta*dble(2*arange(1,Ldelta)-1)
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
  !
  call allocate_dmft_bath(dmft_bath)
  call set_dmft_bath(bath_,dmft_bath)
  !
  !E_{\s,\a}(:)  [ 1 ][ 1 ][Nbath]
  !D_{\s,\a}(:)  [ 1 ][ 1 ][Nbath]
  !V_{\s,\a}(:)  [ 1 ][ 1 ][Nbath]
  Asize = Nbath + Nbath + Nbath
  allocate(array_bath(Asize))
  !
  do jorb=1,Norb
     if(present(iorb))then
        if(jorb/=iorb)cycle
     endif
     Orb_indx=iorb
     Spin_indx=ispin
     !
     Gdelta(1,1:Ldelta) = fg(1,jorb,jorb,1:Ldelta)
     Fdelta(1,1:Ldelta) = fg(2,jorb,jorb,1:Ldelta)
     !
     !3*Nbath == Nbath + Nbath + Nbath
     stride = 0
     do i=1,Nbath
        io = stride + i
        array_bath(io) = dmft_bath%e(ispin,jorb,i)
     enddo
     stride = Nbath
     do i=1,Nbath
        io = stride + i
        array_bath(io) = dmft_bath%d(ispin,jorb,i)
     enddo
     stride = 2*Nbath
     do i=1,Nbath
        io = stride + i
        array_bath(io) = dmft_bath%v(ispin,jorb,i)
     enddo
     !
     !
     select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE
     case default
        if(cg_grad==0)then
           select case (cg_scheme)
           case ("weiss")
              call fmin_cg(array_bath,chi2_weiss_normal_superc,grad_chi2_weiss_normal_superc,&
                   iter,chi,&
                   itmax=cg_niter,&
                   ftol=cg_Ftol,  &
                   istop=cg_stop, &
                   iverbose=(ed_verbose>3))
           case ("delta")
              call fmin_cg(array_bath,chi2_delta_normal_superc,grad_chi2_delta_normal_superc,&
                   iter,chi,&
                   itmax=cg_niter,&
                   ftol=cg_Ftol,  &
                   istop=cg_stop, &
                   iverbose=(ed_verbose>3))
           case default
              stop "chi2_fitgf_normal_superc error: cg_scheme != [weiss,delta]"
           end select
        else
           select case (cg_scheme)
           case ("weiss")
              call fmin_cg(array_bath,chi2_weiss_normal_superc,&
                   iter,chi,&
                   itmax=cg_niter,&
                   ftol=cg_Ftol,  &
                   istop=cg_stop, &
                   iverbose=(ed_verbose>3))
           case ("delta")
              call fmin_cg(array_bath,chi2_delta_normal_superc,&
                   iter,chi,&
                   itmax=cg_niter,&
                   ftol=cg_Ftol,  &
                   istop=cg_stop, &
                   iverbose=(ed_verbose>3))
           case default
              stop "chi2_fitgf_normal_superc error: cg_scheme != [weiss,delta]"
           end select
        endif
        !
     case (1)
        select case (cg_scheme)
        case ("weiss")
           call fmin_cgminimize(array_bath,chi2_weiss_normal_superc,&
                iter,chi,itmax=cg_niter,ftol=cg_Ftol,&
                new_version=cg_minimize_ver,&
                hh_par=cg_minimize_hh,&
                iverbose=(ed_verbose>3))
        case ("delta")
           call fmin_cgminimize(array_bath,chi2_delta_normal_superc,&                
                iter,chi,itmax=cg_niter,ftol=cg_Ftol,&
                new_version=cg_minimize_ver,&
                hh_par=cg_minimize_hh,&
                iverbose=(ed_verbose>3))
        case default
           stop "chi2_fitgf_normal_superc error: cg_scheme != [weiss,delta]"
        end select
        !
        !
     end select
     !
     write(LOGfile,"(A,ES18.9,A,I5,A)")&
          'chi^2|iter'//reg(ed_file_suffix)//'=',chi," | ",iter,&
          "  <--  Orb"//reg(txtfy(jorb))//" Spin"//reg(txtfy(ispin))
     !
     suffix="_orb"//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
     unit=free_unit()
     open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
     write(unit,"(ES18.9,1x,I5)") chi,iter
     close(unit)
     !
     stride = 0
     do i=1,Nbath
        io = stride + i
        dmft_bath%e(ispin,jorb,i) = array_bath(io)
     enddo
     stride = Nbath
     do i=1,Nbath
        io = stride + i
        dmft_bath%d(ispin,jorb,i) = array_bath(io) 
     enddo
     stride = 2*Nbath
     do i=1,Nbath
        io = stride + i
        dmft_bath%v(ispin,jorb,i) = array_bath(io)
     enddo
     !
  enddo
  call write_dmft_bath(dmft_bath,LOGfile)
  !
  call save_dmft_bath(dmft_bath)

  allocate(fgand(Nspin,Nspin,Norb,Norb,Ldelta))
  allocate(ffand(Nspin,Nspin,Norb,Norb,Ldelta))
  if(cg_scheme=='weiss')then
     fgand = g0and_bath_function(xi*Xdelta(:),dmft_bath)
     ffand = f0and_bath_function(xi*Xdelta(:),dmft_bath)
  else
     fgand = delta_bath_function(xi*Xdelta(:),dmft_bath)
     ffand =fdelta_bath_function(xi*Xdelta(:),dmft_bath)
  endif
  call write_fit_result(ispin)
  deallocate(fgand,ffand)
  call get_dmft_bath(dmft_bath,bath_)
  call deallocate_dmft_bath(dmft_bath)
  deallocate(Gdelta,Fdelta,Xdelta,Wdelta)
  !
contains
  !
  subroutine write_fit_result(ispin)
    integer           :: jorb,ispin,gunit,funit
    do jorb=1,Norb
       if(present(iorb))then
          if(jorb/=iorb)cycle
       endif
       suffix="_l"//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
       if(cg_scheme=='weiss')then
          open(free_unit(gunit),file="fit_weiss"//reg(suffix)//".ed")
          open(free_unit(funit),file="fit_fweiss"//reg(suffix)//".ed")
       else
          open(free_unit(gunit),file="fit_delta"//reg(suffix)//".ed")
          open(free_unit(funit),file="fit_fdelta"//reg(suffix)//".ed")
       endif
       do i=1,Ldelta
          write(gunit,"(5F24.15)")Xdelta(i),&
               dimag(fg(1,jorb,jorb,i)),dimag(fgand(ispin,ispin,jorb,jorb,i)),&
               dreal(fg(1,jorb,jorb,i)),dreal(fgand(ispin,ispin,jorb,jorb,i))
          write(gunit,"(5F24.15)")Xdelta(i),&
               dimag(fg(2,jorb,jorb,i)),dimag(ffand(ispin,ispin,jorb,jorb,i)),&
               dreal(fg(2,jorb,jorb,i)),dreal(ffand(ispin,ispin,jorb,jorb,i))
       enddo
       close(gunit)
       close(funit)
    enddo
  end subroutine write_fit_result
end subroutine chi2_fitgf_normal_superc








! !##################################################################
! ! THESE PROCEDURES EVALUATES THE \chi^2 FUNCTIONS TO MINIMIZE. 
! !##################################################################
!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function 
!         in the SUPERCONDUCTING case.
!+-------------------------------------------------------------+
function chi2_delta_normal_superc(a) result(chi2)
  real(8),dimension(:)           ::  a
  real(8)                        ::  chi2
  complex(8),dimension(2,Ldelta) ::  Delta
  real(8),dimension(Ldelta)      ::  Ctmp,Ftmp
  !
  Delta(:,:) = delta_normal_superc(a)
  !
  Ctmp = abs(Gdelta(1,:)-Delta(1,:))
  Ftmp = abs(Fdelta(1,:)-Delta(2,:))
  chi2 = sum( Ctmp**cg_pow/Wdelta ) + sum( Ftmp**cg_pow/Wdelta )
  chi2 = chi2/Ldelta
  !
end function chi2_delta_normal_superc

!+-------------------------------------------------------------+
!PURPOSE: Evaluate the gradient \Grad\chi^2 of \Delta_Anderson 
! function in the SUPERCONDUCTING case.
!+-------------------------------------------------------------+
function grad_chi2_delta_normal_superc(a) result(dchi2)
  real(8),dimension(:)                   ::  a
  real(8),dimension(size(a))             ::  dchi2
  real(8),dimension(size(a))             ::  df
  complex(8),dimension(2,Ldelta)         ::  Delta
  complex(8),dimension(Ldelta)           ::  Gtmp,Ftmp
  real(8),dimension(Ldelta)              ::  Ctmp,Btmp
  complex(8),dimension(2,Ldelta,size(a)) ::  dDelta
  integer                                ::  j
  !
  Delta(:,:)    = delta_normal_superc(a)
  dDelta(:,:,:) = grad_delta_normal_superc(a)
  !
  Gtmp = Gdelta(1,:)-Delta(1,:)
  Ftmp = Fdelta(1,:)-Delta(2,:)
  !
  Ctmp = abs(Gtmp)**(cg_pow-2)
  Btmp = abs(Ftmp)**(cg_pow-2)
  !
  do j=1,size(a)
     df(j) = &
          sum( dreal(Gtmp)*dreal(dDelta(1,:,j))*Ctmp/Wdelta ) + &
          sum( dimag(Gtmp)*dimag(dDelta(1,:,j))*Ctmp/Wdelta ) + &
          sum( dreal(Ftmp)*dreal(dDelta(2,:,j))*Btmp/Wdelta ) + &
          sum( dimag(Ftmp)*dimag(dDelta(2,:,j))*Btmp/Wdelta )
  enddo
  !
  dchi2 = -cg_pow*df/Ldelta
  !
end function grad_chi2_delta_normal_superc

!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance for G_0 function 
!         in the SUPERCONDUCTING case.
!+-------------------------------------------------------------+
function chi2_weiss_normal_superc(a) result(chi2)
  real(8),dimension(:)           ::  a
  complex(8),dimension(2,Ldelta) ::  g0and
  real(8),dimension(Ldelta)    ::  Gtmp,Ftmp
  real(8)                        ::  chi2
  chi2=0d0
  !
  g0and(:,:)  = g0and_normal_superc(a)
  !
  Gtmp = abs(Gdelta(1,:)-g0and(1,:))
  Ftmp = abs(Fdelta(1,:)-g0and(2,:))
  chi2 =        sum( Gtmp**cg_pow/Wdelta )
  chi2 = chi2 + sum( Ftmp**cg_pow/Wdelta )
  chi2 = chi2/Ldelta
  !
end function chi2_weiss_normal_superc

!+-------------------------------------------------------------+
!PURPOSE: Evaluate the gradient \Grad\chi^2 of 
! \Delta_Anderson function.
!+-------------------------------------------------------------+
function grad_chi2_weiss_normal_superc(a) result(dchi2)
  real(8),dimension(:)                   :: a
  real(8),dimension(size(a))             :: dchi2
  real(8),dimension(size(a))             :: df
  complex(8),dimension(2,Ldelta)         :: g0and
  complex(8),dimension(2,Ldelta,size(a)) :: dg0and
  complex(8),dimension(Ldelta)           :: Gtmp,Ftmp
  real(8),dimension(Ldelta)              :: Ctmp,Btmp
  integer                                :: j
  !
  g0and  = g0and_normal_superc(a)
  dg0and = grad_g0and_normal_superc(a)
  !
  Gtmp = abs(Gdelta(1,:)-g0and(1,:))
  Ftmp = abs(Fdelta(1,:)-g0and(2,:))
  !
  Ctmp = abs(Gtmp)**(cg_pow-2)
  Btmp = abs(Ftmp)**(cg_pow-2)
  !
  do j=1,size(a)
     df(j) = &
          sum( dreal(Gtmp)*dreal(dg0and(1,:,j))*Ctmp/Wdelta ) + &
          sum( dimag(Gtmp)*dimag(dg0and(1,:,j))*Ctmp/Wdelta ) + &
          sum( dreal(Ftmp)*dreal(dg0and(2,:,j))*Btmp/Wdelta ) + &
          sum( dimag(Ftmp)*dimag(dg0and(2,:,j))*Btmp/Wdelta )
  enddo
  !
  dchi2 = -cg_pow*df/Ldelta
  !
end function grad_chi2_weiss_normal_superc







!##################################################################
! THESE PROCEDURES EVALUATES THE 
! - \delta
! - \grad \delta
! - g0
! FUNCTIONS. 
!##################################################################
function delta_normal_superc(a) result(Delta)
  real(8),dimension(:)            :: a
  complex(8),dimension(2,Ldelta)  :: Delta
  integer                         :: i,k,io,stride
  real(8),dimension(Nbath)        :: eps,vps,dps
  real(8),dimension(Nbath)        :: Den
  !
  !\Delta_{aa} = - \sum_k [ V_{a}(k) * V_{a}(k) * (iw_n + E_{a}(k)) / Den(k) ]
  !
  stride = 0
  do i=1,Nbath
     io = stride + i
     eps(i) = a(io)
  enddo
  stride = Nbath
  do i=1,Nbath
     io = stride + i
     dps(i) = a(io) 
  enddo
  stride = 2*Nbath
  do i=1,Nbath
     io = stride + i
     vps(i) = a(io)
  enddo
  !
  do i=1,Ldelta
     Delta(1,i) = -sum( vps(:)*vps(:)*( xi*Xdelta(i) + eps(:) )/(Xdelta(i)**2 + eps(:)**2 + dps(:)**2) )
     Delta(2,i) =  sum( dps(:)*vps(:)*vps(:)/(Xdelta(i)**2 + eps(:)**2 + dps(:)**2) )
  enddo
  !
end function delta_normal_superc

function grad_delta_normal_superc(a) result(dDelta)
  real(8),dimension(:)                   :: a
  complex(8),dimension(2,Ldelta,size(a)) :: dDelta
  integer                                :: i,k,ik,io,stride
  real(8),dimension(Nbath)               :: eps,vps,dps
  real(8),dimension(Ldelta,Nbath)        :: Den
  !
  !\grad_{E_{a}(k)} \Delta_{bb} = -V_{a}(k)*V_{a}(k)*[ 1/den(k) - 2*E_{a}(k)*(iw_n + E_{a}(k))/den(k)**2 ]
  !
  !\grad_{\D_{a}(k)} \Delta_{bb} = V_{a}(k)*V_{a}(k)*\D_{a}(k)*(iw_n + E_{a}(k)) /den(k)**2
  !
  !\grad_{ V_{a}(k)} \Delta_{bb} =  2*V_{a}(k)*(iw_n + E_{a}(k))/den(k)
  !
  !
  !
  !\grad_{E_{a}(k)} \FDelta_{aa} = -2 * V_{a}(k) * V_{a}(k) * E_{a}(k) * \Delta_{a}(k) / Den**2
  !
  !\grad_{\Delta_{a}(k)} \FDelta_{aa} = V_{a}(k) * V_{a}(k) * [ 1/den - 2* \Delta_{a}(k)*\Delta_{a}(k)/den**2 ]
  !
  !\grad_{ V_{a}(k)} \FDelta_{aa} =  2 * V_{a}(k) * \Delta_{a}(k) / den
  !
  stride = 0
  do i=1,Nbath
     io = stride + i
     eps(i) = a(io)
  enddo
  stride = Nbath
  do i=1,Nbath
     io = stride + i
     dps(i) = a(io) 
  enddo
  stride = 2*Nbath
  do i=1,Nbath
     io = stride + i
     vps(i) = a(io)
  enddo
  !
  !Den(k) = (w_n**2 + E_{a}(k)**2 + \D_{a}(k)**2
  forall(i=1:Ldelta,k=1:Nbath)Den(i,k) = Xdelta(i)**2 + eps(k)**2 + dps(k)**2 
  !
  stride = 0
  do k=1,Nbath
     ik = stride + k
     dDelta(1,:,ik) = -vps(k)*vps(k)*(1d0/Den(:,k) - 2d0*eps(k)*( xi*Xdelta(:) + eps(k))/Den(:,k)**2)
  enddo
  stride = Nbath
  do k=1,Nbath
     ik = stride + k
     dDelta(1,:,ik) = 2d0*vps(k)*vps(k)*dps(k)*( xi*Xdelta(:) + eps(k))/Den(:,k)**2
  enddo
  stride = 2*Nbath
  do k=1,Nbath
     ik = stride + k
     dDelta(1,:,ik) = -2d0*vps(k)*( xi*Xdelta(:) + eps(k))/Den(:,k)
  enddo
  !
  !
  stride = 0
  do k=1,Nbath
     ik = stride + k
     dDelta(2,:,ik) = -2d0*vps(k)*vps(k)*eps(k)*dps(k)/Den(:,k)**2
  enddo
  stride = Nbath
  do k=1,Nbath
     ik = stride + k
     dDelta(2,:,ik) = vps(k)*vps(k)*(1d0/Den(:,k) - 2d0*dps(k)*dps(k)/Den(:,k)**2)
  enddo
  stride = 2*Nbath
  do k=1,Nbath
     ik = stride + k
     dDelta(2,:,ik) = 2d0*vps(k)*dps(k)/Den(:,k)
  enddo
  !
end function grad_delta_normal_superc

function g0and_normal_superc(a) result(G0and)
  real(8),dimension(:)            :: a
  complex(8),dimension(2,Ldelta)  :: G0and,Delta
  real(8),dimension(Ldelta)       :: det
  complex(8),dimension(Ldelta)    :: fg,ff
  integer                         :: iorb,ispin
  !
  iorb   = Orb_indx
  ispin  = Spin_indx
  !
  Delta    = delta_normal_superc(a)
  !
  fg(:)    = xi*Xdelta(:) + xmu - impHloc(ispin,ispin,iorb,iorb) -  Delta(1,:)
  ff(:)    =                                                     -  Delta(2,:)
  det(:)   = abs(fg(:))**2 + ff(:)**2
  G0and(1,:) = conjg(fg(:))/det(:)
  G0and(2,:) = ff(:)/det(:)
  !
end function g0and_normal_superc

function grad_g0and_normal_superc(a) result(dG0and)
  real(8),dimension(:)                   :: a
  complex(8),dimension(2,Ldelta)         :: G0and,Delta
  complex(8),dimension(2,Ldelta,size(a)) :: dG0and,dDelta
  integer                                :: i,k,ik,io,stride
  real(8),dimension(Nbath)               :: eps,vps,dps
  integer                                :: iorb,ispin
  real(8),dimension(Ldelta,Nbath)        :: Den
  complex(8),dimension(Ldelta)           :: g0,f0,dD,dC,dDet,zeta
  real(8),dimension(Ldelta)              :: det
  !
  iorb   = Orb_indx
  ispin  = Spin_indx
  !
  Delta  = delta_normal_superc(a)
  dDelta = grad_delta_normal_superc(a)
  !
  zeta= xi*Xdelta(:) + xmu - impHloc(ispin,ispin,iorb,iorb) 
  g0  = zeta - Delta(1,:)
  f0  =      - Delta(2,:)
  Det = abs(g0)**2 + f0**2
  !
  do k=1,size(a)
     dD = conjg(dDelta(1,:,k))
     dC = dDelta(2,:,k)
     dDet = 2*dreal(g0*dD)+2*f0*dC
     dG0and(1,:,k) = -Det*dD + conjg(g0)*dDet
     dG0and(2,:,k) = -Det*dC + f0*dDet
     dG0and(1,:,k) = dG0and(1,:,k)/Det**2
     dG0and(2,:,k) = dG0and(2,:,k)/Det**2
  enddo
end function grad_g0and_normal_superc







