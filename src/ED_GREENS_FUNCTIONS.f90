MODULE ED_GREENS_FUNCTIONS
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: str,free_unit,reg,free_units,txtfy,to_lower
  USE SF_LINALG,  only: inv,eigh,eye
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  ! USE ED_IO                     !< this contains the routine to print GF,Sigma and G0
  USE ED_EIGENSPACE
  USE ED_BATH
  USE ED_SETUP
  USE ED_HAMILTONIAN
  USE ED_AUX_FUNX
  !
  USE ED_GF_NORMAL
  USE ED_GF_SUPERC
  USE ED_GF_NONSU2
  !
  implicit none
  private 

  public :: buildGf_impurity

  public :: get_G_impurity
  public :: get_F_impurity
  public :: get_Sigma_impurity
  public :: get_Self_impurity

contains


  !+------------------------------------------------------------------+
  ! GF CALCULATIONS
  !+------------------------------------------------------------------+
  subroutine buildGF_impurity()
    ! 
    ! Build the quantum impurity electrons Green's functions :math:`\hat{G}` , the self-energy :math:`\hat{\Sigma}` and the phonons Green's function :math:`\hat{D}` , calling the correct procedure according to the value of :f:var:`ed_mode` .
    ! Write the results on file according to the values of input variables :f:var:`ed_print_g` , :f:var:`ed_print_sigma` and :f:var:`ed_print_g0` .
    !
    ! * :code:`normal` : :f:func:`build_gf_normal` and :f:func:`build_sigma_normal`
    ! * :code:`superc` : :f:func:`build_gf_superc` and :f:func:`build_sigma_superc`
    ! * :code:`nonsu2` : :f:func:`build_gf_nonsu2` and :f:func:`build_sigma_nonsu2`
    !
    !
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG build_GF: build GFs"
#endif
    call allocate_grids
    !
    call deallocate_GFmatrix(impGmatrix)
    !
    impDmats_ph=zero
    impDreal_ph=zero
    !
    write(LOGfile,"(A)")"Get impurity Greens functions:"
    select case(ed_mode)
    case default  ;call build_Gimp_normal()
    case("superc");call build_Gimp_superc()
    case("nonsu2");call build_Gimp_nonsu2()
    end select
    !
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG build_GF: writing results"
    write(Logfile,"(A)")""
#endif
    if(MPIMASTER)then
       call print_impGmatrix()
       if(ed_print_Sigma)        call print_Sigma()
       if(ed_print_G)            call print_impG()
       if(ed_print_G0)           call print_impG0()
       if(ed_print_G.AND.DimPh>1)call print_impD()
    endif
    !
    call deallocate_grids
    !
  end subroutine buildGF_impurity





  function get_G_impurity(zeta,axis) result(self)
    complex(8),dimension(:),intent(in)                     :: zeta
    character(len=*),optional                              :: axis
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(zeta)) :: self
    character(len=1)                                       :: axis_
    axis_ = 'm' ; if(present(axis))axis_ = axis(1:1)
    select case(ed_mode)
    case default  ;stop "get_Gimp error: not a valid ed_mode"
    case("normal");self = get_Gimp_normal(zeta,axis_)
    case("superc");self = get_Gimp_superc(zeta,axis_)
    case("nonsu2");self = get_Gimp_nonsu2(zeta,axis_)
    end select
  end function get_G_impurity
  !
  function get_F_impurity(zeta,axis) result(self)
    complex(8),dimension(:),intent(in)                     :: zeta
    character(len=*),optional                              :: axis
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(zeta)) :: self
    character(len=1)                                       :: axis_
    axis_ = 'm' ; if(present(axis))axis_ = axis(1:1)
    select case(ed_mode)
    case default  ;stop "get_Fimp error: not a valid ed_mode"
    case("superc");self = get_Fimp_superc(zeta,axis_)
    end select
  end function get_F_impurity





  function get_Sigma_impurity(zeta,axis) result(self)
    complex(8),dimension(:),intent(in)                     :: zeta
    character(len=*),optional                              :: axis
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(zeta)) :: self
    character(len=1)                                       :: axis_
    axis_ = 'm' ; if(present(axis))axis_ = axis(1:1)
    select case(ed_mode)
    case default  ;stop "get_Sigma error: not a valid ed_mode"
    case("normal");self = get_Sigma_normal(zeta,axis_)
    case("superc");self = get_Sigma_superc(zeta,axis_)
    case("nonsu2");self = get_Sigma_nonsu2(zeta,axis_)
    end select
  end function get_Sigma_impurity
  !
  function get_Self_impurity(zeta,axis) result(self)
    complex(8),dimension(:),intent(in)                     :: zeta
    character(len=*),optional                              :: axis
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(zeta)) :: self
    character(len=1)                                       :: axis_
    axis_ = 'm' ; if(present(axis))axis_ = axis(1:1)
    select case(ed_mode)
    case default  ;stop "get_Self error: not a valid ed_mode"
    case("superc");self = get_Self_superc(zeta,axis_)
    end select
  end function get_Self_impurity











  !##################################################################
  !##################################################################
  !                    PRINT FUNCTIONS
  !##################################################################
  !##################################################################
  subroutine print_impGmatrix(file)
    !This subroutine prints weights and poles of the impurity Green's function by calling :f:func:`write_GFmatrix`. These are stored
    !one a file named :code:`"file"//str(ed_file_suffix)//.restart"` taking into account the value of the global variable :f:var:`ed_file_suffix` ,
    !which is :code:`"_ineq_Nineq"` padded with 4 zeros in the case of inequivalent sites, as per documentation
    !
    character(len=*),optional :: file !filename prefix (default :code:`gfmatrix`)
    character(len=256)        :: file_
    if(.not.allocated(impGmatrix))stop "ED_PRINT_IMPGFMATRIX ERROR: impGmatrix not allocated!"
    file_="gfmatrix";if(present(file))file_=str(file)
    call write_GFmatrix(impGmatrix,str(file_)//str(ed_file_suffix)//".restart")
  end subroutine print_impGmatrix





  subroutine print_Sigma
    !This subroutine print the impurity self-energy on plain text files in the execution folder.
    !The files are formatted like :math:`[\omega,\mathrm{Im}\Sigma,\mathrm{Re}\Sigma]` .
    !One file per self-energy component, with the name
    !
    !  * :code:`"impSigma_l"//str(iorb)[str(jorb)]//_s"//str(ispin)[str(jspin)]"_iw"//reg(ed_file_suffix)//".ed"` normal self-energy, Matsubara axis
    !  * :code:`"impSigma_l"//str(iorb)[str(jorb)]//_s"//str(ispin)[str(jspin)]"_realw"//reg(ed_file_suffix)//".ed"` normal self-energy, real frequency axis
    !  * :code:`"impSelf_l"//str(iorb)[str(jorb)]//_s"//str(ispin)"_iw"//reg(ed_file_suffix)//".ed"` anomalous self-energy, Matsubara axis
    !  * :code:`"impSelf_l"//str(iorb)[str(jorb)]//_s"//str(ispin)"_realw"//reg(ed_file_suffix)//".ed"` anomalous self-energy, real frequency axis
    !
    !The variable :f:var:`ed_file_suffix` is :code:`"_ineq_Nineq"` padded with 4 zeros in the case of inequivalent sites, as per documentation.
    !
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: impSmats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: impSAmats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal) :: impSreal
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal) :: impSAreal
    call allocate_grids
    select case(ed_mode)
    case ("normal");
       impSmats  = get_Sigma_impurity(dcmplx(0d0,wm(:)),axis='m')
       impSreal  = get_Sigma_impurity(dcmplx(wr(:),eps),axis='r')
       call Gprint_normal("_impSigma",impSmats,'m')
       call Gprint_normal("_impSigma",impSreal,'r')
    case ("superc");
       impSmats  = get_Sigma_impurity(dcmplx(0d0,wm(:)),axis='m')
       impSreal  = get_Sigma_impurity(dcmplx(wr(:),eps),axis='r')
       impSAmats = get_Self_impurity(dcmplx(0d0,wm(:)),axis='m')
       impSAreal = get_Self_impurity(dcmplx(wr(:),eps),axis='r')
       call Gprint_superc("_impSigma",impSmats,'m')
       call Gprint_superc("_impSelf",impSAmats,'m')
       call Gprint_superc("_impSigma",impSreal,'r')
       call Gprint_superc("_impSelf",impSAreal,'r')
    case ("nonsu2");
       impSmats  = get_Sigma_impurity(dcmplx(0d0,wm(:)),axis='m')
       impSreal  = get_Sigma_impurity(dcmplx(wr(:),eps),axis='r')
       call Gprint_nonsu2("_impSigma",impSmats,'m')
       call Gprint_nonsu2("_impSigma",impSreal,'r')
    case default;stop "ed_print_Sigma error: ed_mode not valid"
    end select
    call deallocate_grids
  end subroutine print_Sigma






  subroutine print_impG
    !This subroutine print the impurity Green's function on plain text files in the execution folder.
    !The files are formatted like :math:`[\omega,\mathrm{Im}G,\mathrm{Re}G]` .
    !One file per Green'sfunction component, with the name
    !
    !  * :code:`"impG_l"//str(iorb)[str(jorb)]//_s"//str(ispin)[str(jspin)]"_iw"//reg(ed_file_suffix)//".ed"` normal G, Matsubara axis
    !  * :code:`"impG_l"//str(iorb)[str(jorb)]//_s"//str(ispin)[str(jspin)]"_realw"//reg(ed_file_suffix)//".ed"` normal G, real frequency axis
    !  * :code:`"impF_l"//str(iorb)[str(jorb)]//_s"//str(ispin)"_iw"//reg(ed_file_suffix)//".ed"` anomalous G, Matsubara axis
    !  * :code:`"impF_l"//str(iorb)[str(jorb)]//_s"//str(ispin)"_realw"//reg(ed_file_suffix)//".ed"` anomalous G, real frequency axis
    !
    !The variable :f:var:`ed_file_suffix` is :code:`"_ineq_Nineq"` padded with 4 zeros in the case of inequivalent sites, as per documentation.
    !
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: impGmats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: impFmats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal) :: impGreal
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal) :: impFreal
    call allocate_grids
    select case(ed_mode)
    case ("normal");
       impGmats  = get_G_impurity(dcmplx(0d0,wm(:)),axis='m')
       impGreal  = get_G_impurity(dcmplx(wr(:),eps),axis='r')
       call Gprint_normal("_impG",impGmats,'m')
       call Gprint_normal("_impG",impGreal,'r')
    case ("superc");
       impGmats  = get_G_impurity(dcmplx(0d0,wm(:)),axis='m')
       impGreal  = get_G_impurity(dcmplx(wr(:),eps),axis='r')
       impFmats = get_F_impurity(dcmplx(0d0,wm(:)),axis='m')
       impFreal = get_F_impurity(dcmplx(wr(:),eps),axis='r')
       call Gprint_superc("_impG",impGmats,'m')
       call Gprint_superc("_impF",impFmats,'m')
       call Gprint_superc("_impG",impGreal,'r')
       call Gprint_superc("_impF",impFreal,'r')
    case ("nonsu2");
       impGmats  = get_G_impurity(dcmplx(0d0,wm(:)),axis='m')
       impGreal  = get_G_impurity(dcmplx(wr(:),eps),axis='r')
       call Gprint_nonsu2("_impG",impGmats,'m')
       call Gprint_nonsu2("_impG",impGreal,'r')
    case default;stop "ed_print_impG error: ed_mode not valid"
    end select
    call deallocate_grids
  end subroutine Print_ImpG





  subroutine ed_print_impG0
    !This subroutine print the non-interacting impurity Green's function on plain text files in the execution folder.
    !The files are formatted like :math:`[\omega,\mathrm{Im}G_{0},\mathrm{Re}G_{0}]` .
    !One file per Green's function component, with the name
    !
    !  * :code:`"impG0_l"//str(iorb)[str(jorb)]//_s"//str(ispin)[str(jspin)]"_iw"//reg(ed_file_suffix)//".ed"` normal G, Matsubara axis
    !  * :code:`"impG0_l"//str(iorb)[str(jorb)]//_s"//str(ispin)[str(jspin)]"_realw"//reg(ed_file_suffix)//".ed"` normal G, real frequency axis
    !  * :code:`"impF0_l"//str(iorb)[str(jorb)]//_s"//str(ispin)"_iw"//reg(ed_file_suffix)//".ed"` anomalous G, Matsubara axis
    !  * :code:`"impF0_l"//str(iorb)[str(jorb)]//_s"//str(ispin)"_realw"//reg(ed_file_suffix)//".ed"` anomalous G, real frequency axis
    !
    !The variable :f:var:`ed_file_suffix` is :code:`"_ineq_Nineq"` padded with 4 zeros in the case of inequivalent sites, as per documentation.
    !
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: impG0mats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: impF0mats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal) :: impG0real
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal) :: impF0real
    call allocate_grids
    select case(ed_mode)
    case ("normal");
       impG0mats  = g0and_bath_function(dcmplx(0d0,wm(:)),dmft_bath,axis='mats')
       impG0real  = g0and_bath_function(dcmplx(wr(:),eps),dmft_bath,axis='real')
       call Gprint_normal("_impG0",impG0mats,'m')
       call Gprint_normal("_impG0",impG0real,'r')
    case ("superc");
       impG0mats  = g0and_bath_function(dcmplx(0d0,wm(:)),dmft_bath,axis='mats')
       impF0mats = f0and_bath_function(dcmplx(0d0,wm(:)),dmft_bath,axis='mats')
       impG0real  = g0and_bath_function(dcmplx(wr(:),eps),dmft_bath,axis='real')
       impF0real = f0and_bath_function(dcmplx(wr(:),eps),dmft_bath,axis='real')
       call Gprint_superc("_impG0",impG0mats,'m')
       call Gprint_superc("_impF0",impF0mats,'m')
       call Gprint_superc("_impG0",impG0real,'r')
       call Gprint_superc("_impF0",impF0real,'r')
    case ("nonsu2");
       impG0mats  = g0and_bath_function(dcmplx(0d0,wm(:)),dmft_bath,axis='mats')
       impG0real  = g0and_bath_function(dcmplx(wr(:),eps),dmft_bath,axis='real')
       call Gprint_nonsu2("_impG0",impG0mats,'m')
       call Gprint_nonsu2("_impG0",impG0real,'r')
    case default;stop "ed_print_impG0 error: ed_mode not valid"
    end select
    call deallocate_grids
  end subroutine Ed_Print_ImpG0






  subroutine ed_print_impD
    !This subroutine print the impurity phonon self-energy on the files
    !  * :code:`"impDph_iw.ed"`  matsubara axis
    !  * :code:`impDph_realw.ed"` real frequency axis
    !
    call allocate_grids()
    !Print the impurity functions:
    call splot("impDph_iw.ed"   ,vm,impDmats_ph(:))
    call splot("impDph_realw.ed",vr,impDreal_ph(:))
    call deallocate_grids()
  end subroutine Ed_Print_ImpD





  !################################################################
  !################################################################
  !################################################################
  !################################################################







  subroutine Gprint_Normal(file,Self,axis)
    character(len=*)                 :: file
    complex(8),dimension(:,:,:,:,:)  :: Self
    character(len=1)                 :: axis
    integer                          :: L,i,ispin,isign,unit(2),iorb,jorb
    character(len=20)                :: suffix
    integer,dimension(:),allocatable :: getIorb,getJorb
    integer                          :: totNorb
    !
    L = size(self,5)
    call assert_shape(self,[Nspin,Nspin,Norb,Norb,L],"Gprint_Normal","Self")
    !
    select case(bath_type)
    case default                !Diagonal in both spin and orbital
       totNorb=Norb
       allocate(getIorb(totNorb),getJorb(totNorb))
       l=0
       do iorb=1,Norb
          L=l+1
          getIorb(l)=iorb
          getJorb(l)=iorb
       enddo
       totNorb=l
    case ("hybrid","replica","general")             !Diagonal in spin only. Full Orbital structure
       totNorb=Norb*(Norb+1)/2
       allocate(getIorb(totNorb),getJorb(totNorb))
       l=0
       do iorb=1,Norb
          do jorb=iorb,Norb
             l=l+1
             getIorb(l)=iorb
             getJorb(l)=jorb
          enddo
       enddo
    end select
    if(l/=totNorb)stop "print_gf_normal error counting the orbitals"
    !!
    !Print the impurity functions:
    do ispin=1,Nspin
       do l=1,totNorb
          iorb=getIorb(l)
          jorb=getJorb(l)
          suffix="_l"//str(iorb)//str(jorb)//"_s"//str(ispin)
          select case(to_lower(axis))
          case default;stop "Gprint_Normal error: axis not supported"
          case("m");call splot(reg(file)//reg(suffix)//"_iw"//reg(ed_file_suffix)//".ed"   ,wm,Self(ispin,ispin,iorb,jorb,:))
          case("r");call splot(reg(file)//reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed",wr,Self(ispin,ispin,iorb,jorb,:))
          end select
       enddo
    enddo
    !
  end subroutine Gprint_Normal





  subroutine Gprint_superc(file,Self,axis)
    character(len=*)                 :: file
    complex(8),dimension(:,:,:,:,:)  :: Self
    character(len=1)                 :: axis
    integer                          :: i,ispin,unit(4),iorb,jorb,isign
    character(len=20)                :: suffix
    integer,dimension(:),allocatable :: getIorb,getJorb
    integer                          :: totNorb,L
    !
    L = size(self,5)
    call assert_shape(self,[Nspin,Nspin,Norb,Norb,L],"Gprint_Superc","Self")
    !
    select case(bath_type)
    case default
       totNorb=Norb
       allocate(getIorb(Norb),getJorb(Norb))
       l=0
       do iorb=1,Norb
          l=l+1
          getIorb(l)=iorb
          getJorb(l)=iorb
       enddo
    case ("hybrid","replica","general")
       totNorb=Norb*(Norb+1)/2
       allocate(getIorb(totNorb),getJorb(totNorb))
       l=0
       do iorb=1,Norb
          do jorb=iorb,Norb
             l=l+1
             getIorb(l)=iorb
             getJorb(l)=jorb
          enddo
       enddo
    end select
    if(l/=totNorb)stop "print_gf_superc error counting the orbitals"
    !!
    !!PRINT OUT GF:
    do ispin=1,Nspin
       do l=1,totNorb
          iorb=getIorb(l)
          jorb=getJorb(l)
          suffix="_l"//str(iorb)//str(jorb)//"_s"//str(ispin)
          select case(to_lower(axis))
          case default;stop "Gprint_Normal error: axis not supported"
          case("m");call splot(reg(file)//reg(suffix)//"_iw"//reg(ed_file_suffix)//".ed"   ,wm,Self(ispin,ispin,iorb,jorb,:))
          case("r");call splot(reg(file)//reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed",wr,Self(ispin,ispin,iorb,jorb,:))
          end select
       enddo
    enddo
    !
  end subroutine Gprint_superc


  subroutine Gprint_nonsu2(file,Self,axis)
    character(len=*)                 :: file
    complex(8),dimension(:,:,:,:,:)  :: Self
    character(len=1)                 :: axis
    integer                          :: i,isign,unit(2),iorb,jorb,ispin,jspin
    integer,dimension(:),allocatable :: getIorb,getJorb,getIspin,getJspin
    integer                          :: totNso,totNorb,totNspin,l,io,jo
    character(len=20)                :: suffix
    !
    L = size(self,5)
    call assert_shape(self,[Nspin,Nspin,Norb,Norb,L],"Gprint_nonSU2","Self")
    !
    select case(bath_type)
    case default
       totNorb =Norb
       totNspin=Nspin*(Nspin+1)/2
       totNso  =totNorb*totNspin
       allocate(getIorb(totNso),getJorb(totNso),getIspin(totNso),getJspin(totNso))
       l=0
       do iorb=1,Norb
          do ispin=1,Nspin
             do jspin=ispin,Nspin
                l=l+1
                getIorb(l)=iorb
                getIspin(l)=ispin
                getJorb(l)=iorb
                getJspin(l)=jspin
             enddo
          enddo
       enddo
    case ("hybrid","replica","general")
       totNso  = (Norb*Nspin)**2
       allocate(getIorb(totNso),getJorb(totNso),getIspin(totNso),getJspin(totNso))
       l=0
       do iorb=1,Norb
          do jorb=1,Norb
             do ispin=1,Nspin
                do jspin=1,Nspin
                   l=l+1
                   getIorb(l)=iorb
                   getIspin(l)=ispin
                   getJorb(l)=jorb
                   getJspin(l)=jspin
                enddo
             enddo
          enddo
       enddo
    end select
    if(l/=totNso)stop "print_gf_nonsu2 error counting the spin-orbitals"
    !!
    !!PRINT OUT GF:
    do l=1,totNso
       iorb=getIorb(l)
       jorb=getJorb(l)
       ispin=getIspin(l)
       jspin=getJspin(l)
       !
       suffix="_l"//str(iorb)//str(jorb)//"_s"//str(ispin)//str(jspin)
       select case(to_lower(axis))
       case default;stop "Gprint_Normal error: axis not supported"
       case("m");call splot(reg(file)//reg(suffix)//"_iw"//reg(ed_file_suffix)//".ed"   ,wm,Self(ispin,ispin,iorb,jorb,:))
       case("r");call splot(reg(file)//reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed",wr,Self(ispin,ispin,iorb,jorb,:))
       end select
    enddo
    !
  end subroutine Gprint_nonsu2


end MODULE ED_GREENS_FUNCTIONS










