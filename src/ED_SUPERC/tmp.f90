
!################################################################
!################################################################
!################################################################
!################################################################




function get_Gimp_superc(zeta,axis) result(Gf)
  !
  ! Reconstructs the system impurity electrons Green's functions using :f:var:`impgmatrix` to retrieve weights and poles.
  !
  complex(8),dimension(:),intent(in)                     :: zeta
  character(len=*),optional                              :: axis
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(zeta)) :: Gf
  integer                                                :: iorb,jorb,ispin,jspin,i
  logical                                                :: MaskBool
  logical(8),dimension(Nspin,Nspin,Norb,Norb)            :: Hmask
  character(len=1)                                       :: axis_
  complex(8)                                             :: barG(Norb,size(zeta)),auxG(4,size(zeta))
  !
#ifdef _DEBUG
  write(Logfile,"(A)")"DEBUG get_Gimp_superc: Get GFs on a input array zeta"
#endif
  !
  axis_ = 'm' ; if(present(axis))axis_ = axis(1:1) !only for self-consistency, not used here
  !
  if(.not.allocated(impGmatrix))stop "get_Gimp_superc ERROR: impGmatrix not allocated!"
  !
  auxG = zero
  barG = zero
  Gf   = zero
  !
  ispin=1                       !in this channel Nspin=2 is forbidden. check in ED_AUX_FUNX.
  !

  do iorb=1,Norb
     auxG(1:2,:) = get_Gimp_superc_Gdiag(iorb)
     Gf(ispin,ispin,iorb,iorb,:) = auxG(1,:) !this is G_{up,up;iorb,iorb}
     barG(               iorb,:) = auxG(2,:) !this is G_{dw,dw;iorb,iorb}
  enddo



  !
  if(offdiag_gf_flag)then
     Hmask= .true.
     if(.not.ed_all_g)then
        if(bath_type=="replica")Hmask=Hreplica_mask(wdiag=.true.,uplo=.false.)
        if(bath_type=="general")Hmask=Hgeneral_mask(wdiag=.true.,uplo=.false.)
     endif
     do ispin=1,Nspin
        do iorb=1,Norb
           write(LOGfile,*)((Hmask(ispin,jspin,iorb,jorb),jorb=1,Norb),jspin=1,Nspin)
        enddo
     enddo
     do ispin=1,Nspin
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              MaskBool=.true.   
              if(bath_type=="replica".OR.bath_type=="general")MaskBool=Hmask(ispin,ispin,iorb,jorb)
              if(.not.MaskBool)cycle
              Gf(ispin,ispin,iorb,jorb,:) = get_Gimp_normal_component(iorb,jorb,ispin)
           enddo
        enddo
     enddo
     !
     do ispin=1,Nspin
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              MaskBool=.true.   
              if(bath_type=="replica".or.bath_type=="general")MaskBool=Hmask(ispin,ispin,iorb,jorb)
              if(.not.MaskBool)cycle
              Gf(ispin,ispin,iorb,jorb,:) = 0.5d0*(Gf(ispin,ispin,iorb,jorb,:) &
                   - Gf(ispin,ispin,iorb,iorb,:) - Gf(ispin,ispin,jorb,jorb,:))
              Gf(ispin,ispin,jorb,iorb,:) = Gf(ispin,ispin,iorb,jorb,:)
           enddo
        enddo
     enddo
  end if
  !
contains
  !
  function get_superc_Gdiag(iorb) result(Gf)
    complex(8),dimension(2,size(zeta)) :: Gf
    integer,intent(in)                 :: iorb
    integer                            :: Nstates,istate
    integer                            :: Nchannels,ichan
    integer                            :: Nexcs,iexc
    real(8)                            :: peso,de
    !
    Gf=zero
    !
    write(LOGfile,"(A)")"Get G_l"//str(iorb)//"_m"//str(iorb)
    if(.not.allocated(impGmatrix(ispin,ispin,iorb,jorb)%state)) then
       write(LOGfile,*)"get_superc_Gdiag: impGmatrix%state not allocated. Nothing to do"
       return
    endif
    !
    Nstates = size(impGmatrix(1,1,iorb,iorb)%state)
    do istate=1,Nstates
       Nchannels = size(impGmatrix(1,1,iorb,iorb)%state(istate)%channel)     
       do ic=1,Nchannels        
          Nexcs  = size(impGmatrix(1,1,iorb,iorb)%state(istate)%channel(ic)%poles)
          if(Nexcs==0)cycle
          select case(ic)
          case(1,2);ichan=1
          case(3,4);ichan=2
          end select
          do iexc=1,Nexcs
             peso  = impGmatrix(1,1,iorb,iorb)%state(istate)%channel(ic)%weight(iexc)
             de    = impGmatrix(1,1,iorb,iorb)%state(istate)%channel(ic)%poles(iexc)
             Gf(ichan,:)=Gf(ichan,:) + peso/(zeta-de)
          enddo
       enddo
    enddo
    return
  end function get_superc_Gdiag
  !
  function get_superc_Gmix(iorb,jorb) result(Gf)
    complex(8),dimension(size(zeta)) :: Gf
    integer,intent(in)               :: iorb,jorb
    integer                          :: Nstates,istate
    integer                          :: Nchannels,ichan
    integer                          :: Nexcs,iexc
    real(8)                          :: peso,de
    !
    Gf=zero
    !
    write(LOGfile,"(A)")"Get G_l"//str(iorb)//"_m"//str(jorb)
    if(.not.allocated(impGmatrix(ispin,ispin,iorb,jorb)%state)) then
       write(LOGfile,*)"get_superc_Gmix: impGmatrix%state not allocated. Nothing to do"
       return
    endif
    !
    Nstates = size(impGmatrix(1,1,iorb,jorb)%state)
    do istate=1,Nstates
       Nchannels = size(impGmatrix(1,1,iorb,jorb)%state(istate)%channel)     
       do ic=1,Nchannels        
          Nexcs  = size(impGmatrix(1,1,iorb,jorb)%state(istate)%channel(ic)%poles)
          if(Nexcs==0)cycle
          do iexc=1,Nexcs
             peso  = impGmatrix(1,1,iorb,jorb)%state(istate)%channel(ic)%weight(iexc)
             de    = impGmatrix(1,1,iorb,jorb)%state(istate)%channel(ic)%poles(iexc)
             Gf    = Gf + peso/(zeta-de)
          enddo
       enddo
    enddo
    return
  end function get_superc_Gmix
  !

  subroutine  get_superc_Fmix(iorb,jorb)
    integer,intent(in) :: iorb,jorb
    integer            :: Nstates,istate
    integer            :: Nchannels,ic,ichan
    integer            :: Nexcs,iexc
    complex(8)         :: peso
    real(8)            :: de
    !
    write(LOGfile,"(A)")"Get F_l"//str(iorb)//"_m"//str(jorb)
    if(.not.allocated(impGmatrix(ispin,ispin,iorb,jorb)%state)) then
       write(LOGfile,*)"get_superc_Fmix: impGmatrix%state not allocated. Nothing to do"
       return
    endif
    !
    Nstates = size(impGmatrix(Nnambu,Nnambu,iorb,jorb)%state)
    do istate=1,Nstates
       Nchannels = size(impGmatrix(Nnambu,Nnambu,iorb,jorb)%state(istate)%channel)     
       do ic=1,Nchannels        
          Nexcs  = size(impGmatrix(Nnambu,Nnambu,iorb,jorb)%state(istate)%channel(ic)%poles)
          if(Nexcs==0)cycle
          do iexc=1,Nexcs
             peso  = impGmatrix(Nnambu,Nnambu,iorb,jorb)%state(istate)%channel(ic)%weight(iexc)
             de    = impGmatrix(Nnambu,Nnambu,iorb,jorb)%state(istate)%channel(ic)%poles(iexc)
             auxGmats(ichan,:)=auxGmats(ichan,:) + peso/(dcmplx(0d0,wm(i))-de)
             auxGreal(ichan,:)=auxGreal(ichan,:) + peso/(dcmplx(wr(i),eps)-de)
          enddo
       enddo
    enddo
    return
  end subroutine get_superc_Fmix

end function get_Gimp_superc




function get_Sigma_normal(zeta,axis) result(Sigma)
  complex(8),dimension(:),intent(in)                     :: zeta
  character(len=*),optional                              :: axis       !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(zeta)) :: Sigma,invG0,invG
  complex(8),dimension(Norb,Norb)                        :: iGzeta
  character(len=4)                                       :: axis_
  !
  axis_="mats";if(present(axis))axis_=str(axis)
  !
  !Get G0^-1
  invG0 = invg0_bath_function(zeta,dmft_bath,axis_)
  !
  !Get G^-1
  invG  = get_Gimp_normal(zeta)
  !
  !Get Sigma= G0^-1 - G^-1
  do ispin=1,Nspin
     do i=1,size(zeta)
        select case(bath_type)
        case ("normal")
           do iorb=1,Norb
              iGzeta(iorb,iorb) = one/invG(ispin,ispin,iorb,iorb,i)
              invG(ispin,ispin,iorb,iorb,i)=iGzeta(iorb,iorb)
           enddo
        case default !Diagonal in spin
           iGzeta(:,:) = invG(ispin,ispin,:,:,i)
           call inv(iGzeta)
           invG(ispin,ispin,:,:,i)=iGzeta
        end select
     enddo
     Sigma(ispin,ispin,:,:,:) = invG0(ispin,ispin,:,:,:) - invG(ispin,ispin,:,:,:)
  enddo
  !
end function get_Sigma_normal
