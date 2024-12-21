subroutine ed_get_rdm_single(dm,doprint)
  complex(8),dimension(:,:),allocatable,intent(inout) :: dm
  logical               ,intent(in) ,optional         :: doprint
  logical                                             :: doprint_
  doprint_=.false.; if(present(doprint)) doprint_=doprint
  if(.not.allocated(impurity_density_matrix))stop "ERROR: impurity_density_matrix is not allocated"
  !
  if(allocated(dm))deallocate(dm)
  allocate(dm, source=impurity_density_matrix)
  !
  if(doprint_)call ed_print_dm(dm,4**Norb)
  !
end subroutine ed_get_rdm_single




subroutine ed_get_reduced_rdm_global(rdm,orbital_mask,doprint)
  complex(8),dimension(:,:),allocatable,intent(inout) :: rdm
  logical,dimension(Norb),intent(in)                  :: orbital_mask
  logical,intent(in),optional                         :: doprint
  logical                                             :: doprint_
  doprint_=.false.; if(present(doprint)) doprint_=doprint
  !
  select case(ed_mode)
  case default  ;call ed_get_reduced_rdm_normal(rdm,orbital_mask,doprint)
  case("superc");stop "it is not implemented"
  case("nonsu2");stop "it is not implemented"
  end select
end subroutine ed_get_reduced_rdm_global


subroutine ed_get_reduced_rdm_normal(rdm,orbital_mask,doprint)
  complex(8),dimension(:,:),allocatable,intent(inout) :: rdm
  logical,dimension(Norb),intent(in)                  :: orbital_mask
  logical,intent(in),optional                         :: doprint
  logical                                             :: doprint_
  logical                                             :: dotrace_
  integer,dimension(:),allocatable                    :: red_indices
  integer,dimension(:),allocatable                    :: trace_indices
  integer,dimension(:),allocatable                    :: IbUP,IbDW
  integer,dimension(:),allocatable                    :: JbUP,JbDW
  integer                                             :: iUP,iDW
  integer                                             :: jUP,jDW
  integer                                             :: iIMPup,iIMPdw
  integer                                             :: jIMPup,jIMPdw
  integer                                             :: iREDup,iREDdw
  integer                                             :: jREDup,jREDdw
  integer                                             :: iTrUP,iTrDW
  integer                                             :: jTrUP,jTrDW
  integer                                             :: Nred,red_count,trace_count,sign
  real(8)                                             :: IsignUP,IsignDW
  real(8)                                             :: JsignUP,JsignDW
  integer                                             :: i,j,io,jo
  !
  Nred = count(orbital_mask)
  if(Nred<1)stop "ERROR: invalid orbital mask, the reduced system must consist of at least one orbital"  
  if(Nred==Norb)then
     dotrace_ = .FALSE.
  else
     dotrace_ = .TRUE.
     allocate(red_indices(Nred),trace_indices(Norb-Nred))
  endif
  doprint_=.false.; if(present(doprint)) doprint_=doprint
  !
  if(allocated(rdm))deallocate(rdm)
  !
  if(.not.allocated(impurity_density_matrix))stop "ERROR: impurity_density_matrix is not allocated"
  !  
  associate(cdm => impurity_density_matrix)
    if(.not.dotrace_)then
       allocate(rdm, source=cdm) 
    else
       ! Retrieve the requested bit-indices for the reduced/traced system
       red_count   = 0
       trace_count = 0
       do io=1,Norb
          if(orbital_mask(io))then
             red_count              = red_count + 1
             red_indices(red_count) = io
          else
             trace_count                = trace_count + 1
             trace_indices(trace_count) = io
          endif
       enddo
       !
       allocate(rdm(4**Nred,4**Nred))
       rdm = zero
       ! Trace the cdm to the requested subsystem, and store into rdm
       do iUP = 1,2**Norb
          IbUP  = bdecomp(iUP-1,Norb)
          call get_sign(IsignUP,IbUP,red_indices)
          call split_state(IbUp,red_indices,trace_indices,iREDup,iTrUP)
          do iDW = 1,2**Norb
             i = iUP + (iDW-1)*2**Norb
             IbDW = bdecomp(iDW-1,Norb)
             call get_sign(IsignDW,IbDW,red_indices)
             call split_state(IbDw,red_indices,trace_indices,iREDdw,iTrDW)
             do JUP = 1,2**Norb
                JbUP  = bdecomp(Jup-1,Norb)
                call get_sign(JsignUP,JbUP,red_indices)
                call split_state(JbUp,red_indices,trace_indices,jREDup,jTrUP)
                do jDW = 1,2**Norb
                   j = jUP + (jDW-1)*2**Norb
                   JbDW = bdecomp(jDW-1,Norb)
                   call get_sign(JsignDW,JbDW,red_indices)
                   call split_state(JbDw,red_indices,trace_indices,jREDdw,jTrDW)
                   if(jTrUP/=iTrUP.or.jTrDW/=iTrDW)cycle
                   io = (iREDup+1) + iREDdw*2**Nred
                   jo = (jREDup+1) + jREDdw*2**Nred
                   sign = IsignUP * IsignDW * JsignUP * JsignDW
                   !
                   rdm(io,jo) = rdm(io,jo) + cdm(i,j) * sign
                   !
                enddo
             enddo
             ! NB: the spin-factorization of the loops is paramount,
             !     since we build in such way the original cdm and the
             !     ordering of the basis states has to be preserved
             !     in order to trace along the proper matrix elements
          enddo
       enddo
    end if
  end associate
  !
  if(doprint_)call ed_print_dm(rdm,orbital_mask)
  !
contains
  !
  ! Compute the Fermionic sign associated to the required swipes
  subroutine get_sign(sign,state,indices)
    real(8), intent(out)                              :: sign
    integer, intent(in)                               :: state(Norb)
    integer, intent(in)                               :: indices(:)
    integer                                           :: filtered(Norb)
    integer                                           :: N
    integer                                           :: r
    ! FILTER THE STATE TO CONSTRAIN THE SUM
    filtered = state; filtered(indices)=0
    ! PERFORM THE SUM (count permutations)
    N = 0
    do r=1,size(indices)
       N = N + sum(filtered(1:indices(r)))
    enddo
    ! ASSIGN THE SIGN: (-1)^N
    if(mod(N,2)==0)then
       sign = 1
    else
       sign = -1
    endif
  end subroutine get_sign
  !
  ! Extract the reduced and tracing subsystems, given the appropriate bit-indices
  subroutine split_state(state,reduced_indices,tracing_indices,reduced_state,tracing_state)
    integer,intent(in),allocatable                    :: state(:)
    integer,intent(in),allocatable                    :: reduced_indices(:)
    integer,intent(in),allocatable                    :: tracing_indices(:)
    integer,intent(out)                               :: reduced_state
    integer,intent(out)                               :: tracing_state
    integer,allocatable                               :: reduced_ibits(:)
    integer,allocatable                               :: tracing_ibits(:)
    !
    reduced_ibits = state(reduced_indices)
    tracing_ibits = state(tracing_indices)
    !
    reduced_state = bjoin(reduced_ibits,Nred)
    tracing_state = bjoin(tracing_ibits,Norb-Nred)
    !
  end subroutine split_state

end subroutine ed_get_reduced_rdm_normal









