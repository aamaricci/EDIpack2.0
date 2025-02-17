MODULE ED_GFMATRIX
  !
  !Contains definition :f:var:`gfmatrix` data structures which
  !constains all the weights and poles of the impurity Green's
  !functions. This is essentially the result of the DMFT calculation,
  !used to generate the Green's functions and the related Self-energy
  !on any given array of points in the complex plane :f:math:`z\in\mathbb{C}`. 
  !
  USE SF_CONSTANTS
  USE SF_IOTOOLS, only:free_unit,reg,str
  USE ED_SPARSE_MATRIX
  USE ED_INPUT_VARS, only:ed_verbose,logfile
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none
  private

  !The contributions to the GF Kallen-Lehmann sum are stored as
  !GF_{ab,sr}%state%channel%{w,e}.
  !A couple of weight,poles {w,e} is stored for each *channel, corresponding to c,cdg or any
  !their combination thereof as well as for any state |n> of the spectrum such that
  !GF(z) = sum w/z-e
  type GFspectrum
     complex(8),dimension(:),allocatable       :: weight
     real(8),dimension(:),allocatable          :: poles
  end type GFspectrum

  !N_channel = c,cdag,c \pm cdag, c \pm i*cdag, ...
  type GFchannel
     type(GFspectrum),dimension(:),allocatable :: channel 
  end type GFchannel

  !state_list%size = # of state in the spectrum 
  type GFmatrix
     type(GFchannel),dimension(:),allocatable  :: state
     logical                                   :: status=.false.
  end type GFmatrix


  interface allocate_GFmatrix
     module procedure :: allocate_GFmatrix_Nstate
     module procedure :: allocate_GFmatrix_Nchan
     module procedure :: allocate_GFmatrix_Nexc
  end interface allocate_GFmatrix


  interface deallocate_GFmatrix
     module procedure :: deallocate_GFmatrix_single
     module procedure :: deallocate_GFmatrix_all1
     module procedure :: deallocate_GFmatrix_all2
     module procedure :: deallocate_GFmatrix_all3
     module procedure :: deallocate_GFmatrix_all4
     module procedure :: deallocate_GFmatrix_all5
  end interface deallocate_GFmatrix

  interface write_GFmatrix
     module procedure :: write_GFmatrix_single
     module procedure :: write_GFmatrix_all1
     module procedure :: write_GFmatrix_all2
     module procedure :: write_GFmatrix_all3
     module procedure :: write_GFmatrix_all4
     module procedure :: write_GFmatrix_all5
  end interface write_GFmatrix

  interface read_GFmatrix
     module procedure :: read_GFmatrix_single
     module procedure :: read_GFmatrix_all1
     module procedure :: read_GFmatrix_all2
     module procedure :: read_GFmatrix_all3
     module procedure :: read_GFmatrix_all4
     module procedure :: read_GFmatrix_all5
  end interface read_GFmatrix

  !EQUALITY with scalar and function (A=B, A=cmplx)
  interface assignment(=)
     module procedure :: GFmatrix_equal_scalar
     module procedure :: GFmatrix_equal_GFmatrix
  end interface assignment(=)




  public :: GFmatrix
  public :: allocate_GFmatrix
  public :: deallocate_GFmatrix
  public :: write_GFmatrix
  public :: read_GFmatrix
  public :: assignment(=)


contains



  !Allocate the channels in GFmatrix structure
  subroutine allocate_gfmatrix_Nstate(self,Nstate)
    type(GFmatrix) :: self
    integer        :: Nstate
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")"DEBUG allocate_gfmatrix_Nstate: allocate self"
#endif
    if(allocated(self%state))deallocate(self%state)
    allocate(self%state(Nstate))
    self%status=.true.
  end subroutine allocate_gfmatrix_Nstate

  subroutine allocate_gfmatrix_Nchan(self,istate,Nchan)
    type(GFmatrix) :: self
    integer        :: istate,Nchan
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG allocate_gfmatrix_Nchan: allocate self, istate:"//str(istate)
#endif
    if(allocated(self%state(istate)%channel))deallocate(self%state(istate)%channel)
    allocate(self%state(istate)%channel(Nchan))
  end subroutine allocate_gfmatrix_Nchan

  !Allocate the Excitations spectrum at a given channel
  subroutine allocate_gfmatrix_Nexc(self,istate,ichan,Nexc)
    type(GFmatrix) :: self
    integer        :: istate,ichan
    integer        :: Nexc
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A,2I8)")"DEBUG allocate_gfmatrix_Nchan: allocate self, istate:"//str(istate)//", ichan:"//str(ichan)
#endif
    if(allocated(self%state(istate)%channel(ichan)%weight))&
         deallocate(self%state(istate)%channel(ichan)%weight)
    if(allocated(self%state(istate)%channel(ichan)%poles))&
         deallocate(self%state(istate)%channel(ichan)%poles)
    !
    allocate(self%state(istate)%channel(ichan)%weight(Nexc))
    allocate(self%state(istate)%channel(ichan)%poles(Nexc))
  end subroutine allocate_gfmatrix_Nexc





  subroutine deallocate_gfmatrix_single(self)
    type(GFmatrix) :: self
    integer        :: istate,ichan
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG deallocate_gfmatrix_single: deallocate self"
#endif
    if(self%status)then
       do istate=1,size(self%state)
          if(allocated(self%state(istate)%channel))then
             do ichan=1,size(self%state(istate)%channel)
                if(allocated(self%state(istate)%channel(ichan)%weight))&
                     deallocate(self%state(istate)%channel(ichan)%weight)
                !
                if(allocated(self%state(istate)%channel(ichan)%poles))&
                     deallocate(self%state(istate)%channel(ichan)%poles)
             enddo
             deallocate(self%state(istate)%channel)
          endif
       enddo
       deallocate(self%state)       
    endif
    self%status=.false.
  end subroutine deallocate_gfmatrix_single

  subroutine deallocate_gfmatrix_all1(self)
    type(GFmatrix),dimension(:) :: self
    integer                     :: i1
    do i1=1,size(self)
       call deallocate_gfmatrix_single(self(i1))
    enddo
  end subroutine deallocate_gfmatrix_all1

  subroutine deallocate_gfmatrix_all2(self)
    type(GFmatrix),dimension(:,:) :: self
    integer                       :: i1,i2
    do i1=1,size(self,1)
       do i2=1,size(self,2)
          call deallocate_gfmatrix_single(self(i1,i2))
       enddo
    enddo
  end subroutine deallocate_gfmatrix_all2

  subroutine deallocate_gfmatrix_all3(self)
    type(GFmatrix),dimension(:,:,:) :: self
    integer                         :: i1,i2,i3
    do i1=1,size(self,1)
       do i2=1,size(self,2)
          do i3=1,size(self,3)
             call deallocate_gfmatrix_single(self(i1,i2,i3))
          enddo
       enddo
    enddo
  end subroutine deallocate_gfmatrix_all3

  subroutine deallocate_gfmatrix_all4(self)
    type(GFmatrix),dimension(:,:,:,:) :: self
    integer                           :: i1,i2,i3,i4
    do i1=1,size(self,1)
       do i2=1,size(self,2)
          do i3=1,size(self,3)
             do i4=1,size(self,4)
                call deallocate_gfmatrix_single(self(i1,i2,i3,i4))
             enddo
          enddo
       enddo
    enddo
  end subroutine deallocate_gfmatrix_all4

  subroutine deallocate_gfmatrix_all5(self)
    type(GFmatrix),dimension(:,:,:,:,:) :: self
    integer                           :: i1,i2,i3,i4,i5
    do i1=1,size(self,1)
       do i2=1,size(self,2)
          do i3=1,size(self,3)
             do i4=1,size(self,4)
                do i5=1,size(self,5)
                   call deallocate_gfmatrix_single(self(i1,i2,i3,i4,i5))
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine deallocate_gfmatrix_all5










  subroutine GFmatrix_equal_GFmatrix(A,B)
    type(GFmatrix),intent(inout) :: A
    type(GFmatrix),intent(in)    :: B
    integer :: Nstate,Nchan,Nexc,istate,ichan,iexc
    if(.not.B%status)stop "GFmatrix_equal_GFmatrix: B.status=F"
    call deallocate_GFmatrix(A)
    Nstate = size(B%state)
    call allocate_GFmatrix(A,Nstate)
    do istate=1,Nstate
       if(.not.allocated(B%state(istate)%channel))cycle
       Nchan = size(B%state(istate)%channel)
       call allocate_GFmatrix(A,istate,Nchan)
       do ichan=1,Nchan
          if(.not.allocated(B%state(istate)%channel(ichan)%poles))cycle
          Nexc = size(B%state(istate)%channel(ichan)%poles)
          call allocate_GFmatrix(A,istate,ichan,Nexc)
          A%state(istate)%channel(ichan)%weight = B%state(istate)%channel(ichan)%weight
          A%state(istate)%channel(ichan)%poles  = B%state(istate)%channel(ichan)%poles
       enddo
    enddo
  end subroutine GFmatrix_equal_GFmatrix


  subroutine GFmatrix_equal_Scalar(A,B)
    type(GFmatrix),intent(inout) :: A
    complex(8),intent(in)        :: B
    integer                      :: Nstate,Nchan,Nexc,istate,ichan,iexc
    if(.not.A%status)stop "GFmatrix_equal_GFmatrix: A.status=F"
    Nstate = size(A%state)
    do istate=1,Nstate
       if(.not.allocated(A%state(istate)%channel))cycle
       Nchan = size(A%state(istate)%channel)
       do ichan=1,Nchan
          if(.not.allocated(A%state(istate)%channel(ichan)%poles))cycle
          Nexc = size(A%state(istate)%channel(ichan)%poles)
          A%state(istate)%channel(ichan)%weight = B
          A%state(istate)%channel(ichan)%poles  = B
       enddo
    enddo
  end subroutine GFmatrix_equal_Scalar





  !+-------------------------------------------------------------------+
  !PURPOSE  : WRITE GFmatrix to file
  !+-------------------------------------------------------------------+
  subroutine write_gfmatrix_single(self,file)
    class(GFmatrix)    :: self
    character(len=*)   :: file
    integer            :: unit_
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG write_gfmatrix_single: write self"
#endif
    unit_=free_unit()
    open(unit_,file=str(file))
    call write_formatted_gfmatrix(self,unit_)
    close(unit_)
  end subroutine write_gfmatrix_single

  subroutine write_gfmatrix_all1(self,file)
    class(GFmatrix)  :: self(:)
    character(len=*) :: file
    integer          :: unit_
    integer          :: i1
    unit_=free_unit()
    open(unit_,file=str(file))
    do i1=1,size(self)
       call write_formatted_gfmatrix(self(i1),unit_)
    enddo
    close(unit_)
  end subroutine write_gfmatrix_all1

  subroutine write_gfmatrix_all2(self,file)
    class(GFmatrix)  :: self(:,:)
    character(len=*) :: file
    integer          :: unit_
    integer          :: i1,i2
    unit_=free_unit()
    open(unit_,file=str(file))
    do i1=1,size(self,1)
       do i2=1,size(self,2)
          call write_formatted_gfmatrix(self(i1,i2),unit_)
       enddo
    enddo
    close(unit_)
  end subroutine write_gfmatrix_all2

  subroutine write_gfmatrix_all3(self,file)
    class(GFmatrix)  :: self(:,:,:)
    character(len=*) :: file
    integer          :: unit_
    integer          :: i1,i2,i3
    unit_=free_unit()
    open(unit_,file=str(file))
    do i1=1,size(self,1)
       do i2=1,size(self,2)
          do i3=1,size(self,3)
             call write_formatted_gfmatrix(self(i1,i2,i3),unit_)
          enddo
       enddo
    enddo
    close(unit_)
  end subroutine write_gfmatrix_all3

  subroutine write_gfmatrix_all4(self,file)
    class(GFmatrix)  :: self(:,:,:,:)
    character(len=*) :: file
    integer          :: unit_
    integer          :: i1,i2,i3,i4
    unit_=free_unit()
    open(unit_,file=str(file))
    do i1=1,size(self,1)
       do i2=1,size(self,2)
          do i3=1,size(self,3)
             do i4=1,size(self,4)
                call write_formatted_gfmatrix(self(i1,i2,i3,i4),unit_)
             enddo
          enddo
       enddo
    enddo
    close(unit_)
  end subroutine write_gfmatrix_all4

  subroutine write_gfmatrix_all5(self,file)
    class(GFmatrix)  :: self(:,:,:,:,:)
    character(len=*) :: file
    integer          :: unit_
    integer          :: i1,i2,i3,i4,i5
    unit_=free_unit()
    open(unit_,file=str(file))
    do i1=1,size(self,1)
       do i2=1,size(self,2)
          do i3=1,size(self,3)
             do i4=1,size(self,4)
                do i5=1,size(self,5)
                   call write_formatted_gfmatrix(self(i1,i2,i3,i4,i5),unit_)
                enddo
             enddo
          enddo
       enddo
    enddo
    close(unit_)
  end subroutine write_gfmatrix_all5






  !+-------------------------------------------------------------------+
  !PURPOSE  : Read cluster GF from file
  !+-------------------------------------------------------------------+
  subroutine read_gfmatrix_single(self,file)
    class(GFmatrix)  :: self
    character(len=*) :: file
    integer          :: unit_
    integer          :: i1,i2,i3,i4
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG reading_gfmatrix_single: "//str(file)
#endif
    call deallocate_GFmatrix(self)
    unit_=free_unit()
    open(unit_,file=str(file))
    rewind(unit_)
    call read_formatted_gfmatrix(self,unit_)
    close(unit_)
  end subroutine read_gfmatrix_single

  subroutine read_gfmatrix_all1(self,file)
    class(GFmatrix)  :: self(:)
    character(len=*) :: file
    integer          :: unit_
    integer          :: i1
    call deallocate_GFmatrix(self)
    unit_=free_unit()
    open(unit_,file=str(file))
    rewind(unit_)
    do i1=1,size(self)
       call read_formatted_gfmatrix(self(i1),unit_)
    enddo
    close(unit_)
  end subroutine read_gfmatrix_all1

  subroutine read_gfmatrix_all2(self,file)
    class(GFmatrix)  :: self(:,:)
    character(len=*) :: file
    integer          :: unit_
    integer          :: i1,i2
    call deallocate_GFmatrix(self)
    unit_=free_unit()
    open(unit_,file=str(file))
    rewind(unit_)
    do i1=1,size(self,1)
       do i2=1,size(self,2)
          call read_formatted_gfmatrix(self(i1,i2),unit_)
       enddo
    enddo
    close(unit_)
  end subroutine read_gfmatrix_all2

  subroutine read_gfmatrix_all3(self,file)
    class(GFmatrix)  :: self(:,:,:)
    character(len=*) :: file
    integer          :: unit_
    integer          :: i1,i2,i3
    call deallocate_GFmatrix(self)
    unit_=free_unit()
    open(unit_,file=str(file))
    rewind(unit_)
    do i1=1,size(self,1)
       do i2=1,size(self,2)
          do i3=1,size(self,3)
             call read_formatted_gfmatrix(self(i1,i2,i3),unit_)
          enddo
       enddo
    enddo
    close(unit_)
  end subroutine read_gfmatrix_all3

  subroutine read_gfmatrix_all4(self,file)
    class(GFmatrix)  :: self(:,:,:,:)
    character(len=*) :: file
    integer          :: unit_
    integer          :: i1,i2,i3,i4
    call deallocate_GFmatrix(self)
    unit_=free_unit()
    open(unit_,file=str(file))
    rewind(unit_)
    do i1=1,size(self,1)
       do i2=1,size(self,2)
          do i3=1,size(self,3)
             do i4=1,size(self,4)
                call read_formatted_gfmatrix(self(i1,i2,i3,i4),unit_)
             enddo
          enddo
       enddo
    enddo
    close(unit_)
  end subroutine read_gfmatrix_all4

  subroutine read_gfmatrix_all5(self,file)
    class(GFmatrix)  :: self(:,:,:,:,:)
    character(len=*) :: file
    integer          :: unit_
    integer          :: i1,i2,i3,i4,i5
    call deallocate_GFmatrix(self)
    unit_=free_unit()
    open(unit_,file=str(file))
    rewind(unit_)
    do i1=1,size(self,1)
       do i2=1,size(self,2)
          do i3=1,size(self,3)
             do i4=1,size(self,4)
                do i5=1,size(self,5)
                   call read_formatted_gfmatrix(self(i1,i2,i3,i4,i5),unit_)
                enddo
             enddo
          enddo
       enddo
    enddo
    close(unit_)
  end subroutine read_gfmatrix_all5






  !+-------------------------------------------------------------------+
  !PURPOSE  : write overload for GFmatrix type (formatted)
  !+-------------------------------------------------------------------+
  subroutine write_formatted_gfmatrix(dtv, unit)
    class(GFmatrix), intent(in)         :: dtv
    integer, intent(in)                 :: unit
    integer                             :: iexc,Ichan,istate
    integer                             :: Nexc,Nchan,Nstates
    write(unit,*) dtv%status
    if(.not.dtv%status)return
    Nstates = size(dtv%state)
    write(unit,*) Nstates
    do istate=1,Nstates
       Nchan = size(dtv%state(istate)%channel)
       write(unit,*)Nchan
       do ichan=1,Nchan
          write(unit,*) size(dtv%state(istate)%channel(ichan)%poles)
          write(unit,*) dtv%state(istate)%channel(ichan)%poles
          write(unit,*) dtv%state(istate)%channel(ichan)%weight
       enddo
    enddo
    write(unit,*)""
  end subroutine write_formatted_gfmatrix

  !+-------------------------------------------------------------------+
  !PURPOSE  : read overload for GFmatrix type (formatted)
  !+-------------------------------------------------------------------+
  subroutine read_formatted_gfmatrix(dtv, unit)
    class(GFmatrix), intent(inout)                :: dtv
    integer, intent(in)                           :: unit
    logical                                       :: alloc
    integer                                       :: ichan,Nchan,Nlanc,istate,Nstates
    !
    read(unit,*) alloc
    if(.not.alloc)return
    read(unit,*)Nstates
    call allocate_GFmatrix(dtv,Nstate=Nstates)
    do istate=1,Nstates
       read(unit,*)Nchan
       call allocate_GFmatrix(dtv,istate=istate,Nchan=Nchan)
       do ichan=1,Nchan
          read(unit,*)Nlanc
          call allocate_GFmatrix(dtv,istate=istate,ichan=ichan,Nexc=Nlanc)
          read(unit,*) dtv%state(istate)%channel(ichan)%poles
          read(unit,*) dtv%state(istate)%channel(ichan)%weight
       enddo
    enddo
    !
  end subroutine read_formatted_gfmatrix


END MODULE ED_GFMATRIX
