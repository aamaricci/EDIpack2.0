MODULE ED_BATH_FUNCTIONS
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,txtfy,str
  USE SF_LINALG, only: eye,inv,zeye,inv_her
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_BATH
  USE ED_AUX_FUNX
  implicit none

  private


  !
  !\DELTA, THE HYBRIDIZATION FUNCTION
  interface delta_bath_function
     module procedure delta_bath_array
  end interface delta_bath_function
  interface Fdelta_bath_function
     module procedure Fdelta_bath_array
  end interface Fdelta_bath_function
  !
  !NON-INTERACTING GREEN'S FUNCTION 
  interface g0and_bath_function
     module procedure g0and_bath_array
  end interface g0and_bath_function
  interface f0and_bath_function
     module procedure f0and_bath_array
  end interface f0and_bath_function
  !
  !INVERSE NON-INTERACTING GREEN'S FUNCTION 
  interface invg0_bath_function
     module procedure invg0_bath_array
  end interface invg0_bath_function
  interface invf0_bath_function
     module procedure invf0_bath_array
  end interface invf0_bath_function

  public :: delta_bath_function,fdelta_bath_function
  public :: g0and_bath_function,f0and_bath_function
  public :: invg0_bath_function,invf0_bath_function






contains



  function delta_bath_array(x,dmft_bath_,axis) result(Delta)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    character(len=*),optional                           :: axis    
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Delta
    integer                                             :: i,ih,L
    integer                                             :: iorb,jorb,ispin,jspin,ibath
    integer                                             :: io,jo
    real(8),dimension(Nbath)                            :: eps,dps,vps
    real(8),dimension(Norb,Nbath)                       :: vops
    !
    real(8),dimension(Nspin,Nbath)                      :: ehel
    real(8),dimension(Nspin,Nspin,Nbath)                :: whel
    real(8),dimension(Nspin,Nspin,Norb,Nbath)           :: wohel
    !
    complex(8),dimension(Nspin*Norb,Nspin*Norb)         :: invH_k
    complex(8),dimension(Nspin,Nspin,Norb,Norb)         :: invH_knn
    character(len=4)                                    :: axis_
    !
    axis_="mats";if(present(axis))axis_=str(axis)
    !
    Delta=zero
    !
    L = size(x)
    !
    select case(bath_type)
       !NORMAL:
       !\Delta_{aa} = \sum_k [ V_{a}(k) * V_{a}(k)/(x - E_{a}(k)) ]
       !SUPERC:
       !IF(mats):
       ! \Delta_{aa} = - \sum_k [ V_{a}(k) * V_{a}(k) * (iw_n + E_{a}(k)) / Den(k) ]
       !ELSE:
       ! \Delta_{aa}^{ss} = - \sum_k [ V_{a}(k) * V_{a}(k) * (w+i\h + E_{a}(k)) / ((w+i\h)*(-w-i\h) + E(k)**2 + \D(k)**2 ]
       !NONSU2:
       !\Delta_{aa}^{ss`} = \sum_h \sum_k [ W_{a}^{sh}(k) * W_{a}^{s`h}(k)/(x - H_{a}^{h}(k))]
    case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
       select case(ed_mode)
       case default
          do ispin=1,Nspin
             do iorb=1,Norb
                eps = dmft_bath_%e(ispin,iorb,1:Nbath)
                vps = dmft_bath_%v(ispin,iorb,1:Nbath)
                do i=1,L
                   Delta(ispin,ispin,iorb,iorb,i) = sum( vps(:)*vps(:)/(x(i) - eps(:)) )
                enddo
             enddo
          enddo
          !
       case ("superc")
          do ispin=1,Nspin
             do iorb=1,Norb
                eps = dmft_bath_%e(ispin,iorb,1:Nbath)
                dps = dmft_bath_%d(ispin,iorb,1:Nbath)
                vps = dmft_bath_%v(ispin,iorb,1:Nbath)
                select case(axis_)
                case default
                   do i=1,L
                      Delta(ispin,ispin,iorb,iorb,i) = -sum( vps(:)*vps(:)*(x(i) + eps(:))/(dimag(x(i))**2 + eps(:)**2 + dps(:)**2) )
                   enddo
                case ("real")
                   do i=1,L
                      Delta(ispin,ispin,iorb,iorb,i) = -sum( vps(:)*vps(:)*(x(i) + eps(:))/(x(i)*(-x(i)) + eps(:)**2 + dps(:)**2) )
                   enddo
                end select
             enddo
          enddo
          !
       case ("nonsu2")
          do iorb=1,Norb
             ehel = dmft_bath_%e(1:Nspin,iorb,1:Nbath)
             whel = get_Whyb_matrix(dmft_bath_%v(1:Nspin,iorb,1:Nbath),dmft_bath_%u(1:Nspin,iorb,1:Nbath))
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do i=1,L
                      do ih=1,Nspin
                         Delta(ispin,jspin,iorb,iorb,i) = Delta(ispin,jspin,iorb,iorb,i) + &
                              sum( whel(ispin,ih,:)*whel(jspin,ih,:)/(x(i) - ehel(ih,:)) )
                      enddo
                   enddo
                enddo
             enddo
          enddo
          !
       end select
       !
       !NORMAL:
       !\Delta_{ab} = \sum_k [ V_{a}(k) * V_{b}(k)/(iw_n - E(k)) ]
       !SUPERC:
       !IF(MATS):
       ! \Delta_{ab} = - \sum_k [ V_{a}(k) * V_{b}(k) * (iw_n + E(k)) / Den(k) ]
       !ELSE:
       ! \Delta_{ab} = - \sum_k [ V_{a}(k) * V_{b}(k) * (w+i\h + E(k)) / ((w+i\h)*(-w-i\h) + E(k)**2 + \D(k)**2 ]
       !NONSU2:
       !\Delta_{ab}^{ss`} = \sum_h \sum_k [ W_{a}^{sh}(k) * W_{b}^{s`h}(k)/(x - H^{h}(k))]
    case ("hybrid")
       select case(ed_mode)
       case default
          do ispin=1,Nspin
             eps  = dmft_bath_%e(ispin,1     ,1:Nbath)
             vops = dmft_bath_%v(ispin,1:Norb,1:Nbath)
             do iorb=1,Norb
                do jorb=1,Norb
                   do i=1,L
                      Delta(ispin,ispin,iorb,jorb,i) = sum( vops(iorb,:)*vops(jorb,:)/(x(i) - eps(:)) )
                   enddo
                enddo
             enddo
          enddo
          !
       case ("superc")
          do ispin=1,Nspin
             eps  = dmft_bath_%e(ispin,1     ,1:Nbath)
             dps  = dmft_bath_%d(ispin,1     ,1:Nbath)
             vops = dmft_bath_%v(ispin,1:Norb,1:Nbath)
             do iorb=1,Norb
                do jorb=1,Norb
                   select case(axis_)
                   case default
                      do i=1,L
                         Delta(ispin,ispin,iorb,jorb,i) = &
                              -sum( vops(iorb,:)*vops(jorb,:)*(x(i) + eps(:))/(dimag(x(i))**2 + eps(:)**2 + dps(:)**2) )
                      enddo
                   case ("real")
                      do i=1,L
                         Delta(ispin,ispin,iorb,jorb,i) = &
                              -sum( vops(iorb,:)*vops(jorb,:)*(x(i) + eps(:))/((x(i)*(-x(i)) + eps(:)**2 + dps(:)**2)) )
                      enddo
                   end select
                enddo
             enddo
          enddo
          !
       case ("nonsu2")
          ehel  = dmft_bath_%e(1:Nspin,1,1:Nbath)
          wohel = get_Whyb_matrix(dmft_bath_%v(1:Nspin,1:Norb,1:Nbath),dmft_bath_%u(1:Nspin,1:Norb,1:Nbath))
          do iorb=1,Norb
             do jorb=1,Norb
                do ispin=1,Nspin
                   do jspin=1,Nspin
                      do i=1,L
                         do ih=1,Nspin
                            Delta(ispin,jspin,iorb,jorb,i) = Delta(ispin,jspin,iorb,jorb,i) + &
                                 sum( wohel(ispin,ih,iorb,:)*wohel(jspin,ih,jorb,:)/(x(i) - ehel(ih,:)) )
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
          !
       end select
       !
    case ("replica")
       invH_k=zero
       do i=1,L
          do ibath=1,Nbath
             invH_knn = Hreplica_build(dmft_bath_%item(ibath)%lambda)
             invH_k   = nn2so_reshape(invH_knn,Nspin,Norb)
             invH_k   = zeye(Nspin*Norb)*x(i) - invH_k
             call inv(invH_k)
             invH_knn = so2nn_reshape(invH_k,Nspin,Norb)
             Delta(:,:,:,:,i)=Delta(:,:,:,:,i) + &
                  dmft_bath_%item(ibath)%v*invH_knn*dmft_bath_%item(ibath)%v
          enddo
          !
       enddo
       !
    end select
  end function delta_bath_array


  !ANOMALous:
  function fdelta_bath_array(x,dmft_bath_,axis) result(Fdelta)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Fdelta
    integer                                             :: iorb,ispin,jorb
    real(8),dimension(Nbath)                            :: eps,dps,vps
    real(8),dimension(Norb,Nbath)                       :: vops
    integer                                             :: i,L
    character(len=*),optional                           :: axis    
    character(len=4)                                    :: axis_
    !
    axis_="mats";if(present(axis))axis_=str(axis)
    !
    Fdelta=zero
    !
    L = size(x)
    !
    select case(bath_type)
    case ("replica")
       stop "Fdelta_bath_mats error: called with bath_type=replica"
    case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
       select case(ed_mode)
       case default
          stop "Fdelta_bath_mats error: called with ed_mode=normal/nonsu2, bath_type=normal"
          !
          !IF(MATS):          
          !\FDelta_{aa} = \sum_k [ \Delta_{a}(k) * V_{a}(k) * V_{a}(k)  / Den(k) ]
          !ELSE:
          !\FDelta_{aa} = \sum_k [ \Delta_{a}(k) * V_{a}(k) * V_{a}(k)  / ((w+i\h)*(-w-i\h) + E(k)**2 + \D(k)**2) ]
       case ("superc")
          do ispin=1,Nspin
             do iorb=1,Norb
                eps = dmft_bath_%e(ispin,iorb,1:Nbath)
                dps = dmft_bath_%d(ispin,iorb,1:Nbath)
                vps = dmft_bath_%v(ispin,iorb,1:Nbath)
                select case(axis_)
                case default
                   do i=1,L
                      Fdelta(ispin,ispin,iorb,iorb,i) = sum( dps(:)*vps(:)*vps(:)/(dimag(x(i))**2 + eps(:)**2 + dps(:)**2) )
                   enddo
                case ("real")
                   do i=1,L
                      Fdelta(ispin,ispin,iorb,iorb,i) = sum( dps(:)*vps(:)*vps(:)/( x(i)*(-x(i)) + eps(:)**2 + dps(:)**2) )
                   enddo
                end select
             enddo
          enddo
          !
       end select
       !
    case ("hybrid")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
       select case(ed_mode)
       case default
          stop "Fdelta_bath_mats error: called with ed_mode=normal/nonsu2, bath_type=hybrid"
          !
       case ("superc")
          !
          !IF(MATS):          
          !\FDelta_{aa} = \sum_k [ \Delta_{a}(k) * V_{a}(k) * V_{a}(k)  / Den(k) ]
          !ELSE:
          !\FDelta_{aa} = \sum_k [ \Delta_{a}(k) * V_{a}(k) * V_{a}(k)  / ((w+i\h)*(-w-i\h) + E(k)**2 + \D(k)**2) ]
          do ispin=1,Nspin
             eps  = dmft_bath_%e(ispin,1     ,1:Nbath)
             dps  = dmft_bath_%d(ispin,1     ,1:Nbath)
             vops = dmft_bath_%v(ispin,1:Norb,1:Nbath)
             do iorb=1,Norb
                do jorb=1,Norb
                   select case(axis_)
                   case default
                      do i=1,L
                         Fdelta(ispin,ispin,iorb,jorb,i) = &
                              -sum( dps(:)*vops(iorb,:)*vops(jorb,:)/(dimag(x(i))**2+eps(:)**2+dps(:)**2))
                      enddo
                   case ("real")
                      do i=1,L
                         Fdelta(ispin,ispin,iorb,jorb,i) = &
                              -sum( dps(:)*vops(iorb,:)*vops(jorb,:)/(x(i)*(-x(i)) + eps(:)**2 + dps(:)**2) )
                      enddo
                   end select
                enddo
             enddo
          enddo
          !
       end select
    end select
  end function fdelta_bath_array







  function g0and_bath_array(x,dmft_bath_,axis) result(G0and)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and,Delta,Fdelta
    integer                                             :: iorb,jorb,ispin,jspin,io,jo,Nso,i,L
    real(8),dimension(size(x))                          :: det
    complex(8),dimension(size(x))                       :: fg,ff
    complex(8),dimension(:,:),allocatable               :: fgorb,zeta
    character(len=*),optional                           :: axis    
    character(len=4)                                    :: axis_
    !
    axis_="mats";if(present(axis))axis_=str(axis)
    !
    G0and = zero
    !
    L=size(x)
    !
    select case(bath_type)
    case default
       select case(ed_mode)
       case default
          Delta = delta_bath_array(x,dmft_bath_)
          do ispin=1,Nspin
             do iorb=1,Norb
                fg(:) = x(:) + xmu - impHloc(ispin,ispin,iorb,iorb) - Delta(ispin,ispin,iorb,iorb,:)
                G0and(ispin,ispin,iorb,iorb,:) = one/fg(:)
             enddo
          enddo
          !
       case ("superc")
          Delta =  delta_bath_array(x,dmft_bath_,axis_)
          Fdelta= fdelta_bath_array(x,dmft_bath_,axis_)
          select case(axis_)
          case default
             do ispin=1,Nspin
                do iorb=1,Norb
                   fg(:)  = x(:) + xmu - impHloc(ispin,ispin,iorb,iorb) -  Delta(ispin,ispin,iorb,iorb,:)
                   ff(:)  =                                             - Fdelta(ispin,ispin,iorb,iorb,:)
                   det(:) = abs(fg(:))**2 + ff(:)*ff(:)
                   G0and(ispin,ispin,iorb,iorb,:) = conjg(fg(:))/det(:)
                enddo
             enddo
          case("real")
             do ispin=1,Nspin
                do iorb=1,Norb
                   fg(:)  =  dreal(x(:)) + xmu - impHloc(ispin,ispin,iorb,iorb) -  Delta(ispin,ispin,iorb,iorb,:)
                   ff(:)  =                                              - Fdelta(ispin,ispin,iorb,iorb,:)
                   det(:) = -fg(:)*conjg(fg(L:1:-1)) - ff(:)*ff(:)
                   G0and(ispin,ispin,iorb,iorb,:) = conjg(fg(L:1:-1))/det(:)
                enddo
             enddo
          end select
          !
       case ("nonsu2")
          allocate(fgorb(Nspin,Nspin),zeta(Nspin,Nspin))
          Delta = delta_bath_array(x,dmft_bath_)
          do i=1,L
             zeta  = (x(i) + xmu)*zeye(Nspin)
             fgorb = zero
             do iorb=1,Norb
                do ispin=1,Nspin
                   do jspin=1,Nspin
                      fgorb(ispin,jspin) = zeta(ispin,jspin) - impHloc(ispin,jspin,iorb,iorb) - Delta(ispin,jspin,iorb,iorb,i)
                   enddo
                enddo
                call inv(fgorb)
                do ispin=1,Nspin
                   do jspin=1,Nspin
                      G0and(ispin,jspin,iorb,iorb,i) = fgorb(ispin,jspin)
                   enddo
                enddo
             enddo
          enddo
          deallocate(fgorb,zeta)
          !
       end select
       !
       !
    case ("hybrid","replica")
       select case(ed_mode)
       case default
          allocate(fgorb(Norb,Norb),zeta(Norb,Norb))
          Delta = delta_bath_array(x,dmft_bath_)
          do ispin=1,Nspin
             do i=1,L
                fgorb = (x(i)+xmu)*zeye(Norb) - impHloc(ispin,ispin,:,:) - Delta(ispin,ispin,:,:,i)
                call inv(fgorb)
                G0and(ispin,ispin,:,:,i)=fgorb
             enddo
          enddo
          deallocate(fgorb,zeta)
          !
       case ("superc")
          allocate(fgorb(2*Norb,2*Norb),zeta(2*Norb,2*Norb))
          Delta =  delta_bath_array(x,dmft_bath_,axis_)
          Fdelta= fdelta_bath_array(x,dmft_bath_,axis_)
          select case(axis_)
          case default
             do ispin=1,Nspin
                do i=1,L
                   zeta = zero
                   fgorb= zero
                   do iorb=1,Norb
                      zeta(iorb,iorb)           = x(i) + xmu
                      zeta(iorb+Norb,iorb+Norb) = x(i) - xmu
                   enddo
                   do iorb=1,Norb
                      do jorb=1,Norb
                         fgorb(iorb,jorb)           = zeta(iorb,jorb)           - impHloc(ispin,ispin,iorb,jorb)  - Delta(ispin,ispin,iorb,jorb,i)
                         fgorb(iorb,jorb+Norb)      = zeta(iorb,jorb+Norb)                                        - Fdelta(ispin,ispin,iorb,jorb,i)
                         fgorb(iorb+Norb,jorb)      = zeta(iorb+Norb,jorb)                                        - Fdelta(ispin,ispin,iorb,jorb,i)
                         fgorb(iorb+Norb,jorb+Norb) = zeta(iorb+Norb,jorb+Norb) + impHloc(ispin,ispin,iorb,jorb)  + conjg( Delta(ispin,ispin,iorb,jorb,i) )
                      enddo
                   enddo
                   call inv(fgorb)
                   G0and(ispin,ispin,:,:,i) = fgorb(1:Norb,1:Norb)
                enddo
             enddo
             !
          case("real")
             do ispin=1,Nspin
                do i=1,L
                   zeta = zero
                   fgorb= zero
                   do iorb=1,Norb
                      zeta(iorb,iorb)           =        x(i)     + xmu
                      zeta(iorb+Norb,iorb+Norb) = -conjg(x(L-i+1) + xmu)
                   enddo
                   do iorb=1,Norb
                      do jorb=1,Norb
                         fgorb(iorb,jorb)           = zeta(iorb,jorb)           - impHloc(ispin,ispin,iorb,jorb)  - Delta(ispin,ispin,iorb,jorb,i)
                         fgorb(iorb,jorb+Norb)      = zeta(iorb,jorb+Norb)                                        - Fdelta(ispin,ispin,iorb,jorb,i)
                         fgorb(iorb+Norb,jorb)      = zeta(iorb+Norb,jorb)                                        - Fdelta(ispin,ispin,iorb,jorb,i)
                         fgorb(iorb+Norb,jorb+Norb) = zeta(iorb+Norb,jorb+Norb) + impHloc(ispin,ispin,iorb,jorb)  + conjg( Delta(ispin,ispin,iorb,jorb,L-i+1) )
                      enddo
                   enddo
                   call inv(fgorb)
                   G0and(ispin,ispin,:,:,i) = fgorb(1:Norb,1:Norb)
                enddo
             enddo
          end select
          deallocate(fgorb,zeta)
          !
       case ("nonsu2")
          !
          Nso=Nspin*Norb
          allocate(fgorb(Nso,Nso),zeta(Nso,Nso));fgorb=zero
          Delta = delta_bath_array(x,dmft_bath_)
          do i=1,L
             zeta  = (x(i) + xmu)*zeye(Nso)
             fgorb = zero
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         io = iorb + (ispin-1)*Norb
                         jo = jorb + (jspin-1)*Norb
                         fgorb(io,jo) = zeta(io,jo) - impHloc(ispin,jspin,iorb,jorb) - Delta(ispin,jspin,iorb,jorb,i)
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
          deallocate(fgorb,zeta)
          !
       end select
       !
    end select
  end function g0and_bath_array


  !ANOMALous:
  function f0and_bath_array(x,dmft_bath_,axis) result(F0and)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: F0and,Delta,Fdelta
    integer                                             :: iorb,jorb,ispin,i,L
    real(8),dimension(size(x))                          :: det
    complex(8),dimension(size(x))                       :: fg,ff
    complex(8),dimension(:,:),allocatable               :: fgorb,zeta
    character(len=*),optional                           :: axis    
    character(len=4)                                    :: axis_
    !
    axis_="mats";if(present(axis))axis_=str(axis)
    !
    F0and=zero
    !
    L = size(x)
    !
    select case(bath_type)
    case default
       select case(ed_mode)
       case default
          stop "F0and_bath_mats error: called with ed_mode=normal/nonsu2, bath_type=normal"
          !
       case ("superc")
          Delta =  delta_bath_array(x,dmft_bath_,axis_)
          Fdelta= fdelta_bath_array(x,dmft_bath_,axis_)
          select case(axis_)
          case default          
             do ispin=1,Nspin
                do iorb=1,Norb
                   fg(:) = x(:) + xmu - impHloc(ispin,ispin,iorb,iorb) -  Delta(ispin,ispin,iorb,iorb,:)
                   ff(:) =                                              - Fdelta(ispin,ispin,iorb,iorb,:)
                   det(:)= abs(fg(:))**2 + ff(:)*ff(:)
                   F0and(ispin,ispin,iorb,iorb,:) = ff(:)/det(:)
                enddo
             enddo
          case("real")
             do ispin=1,Nspin
                do iorb=1,Norb
                   fg(:)  = dreal(x(:)) + xmu - impHloc(ispin,ispin,iorb,iorb) -  Delta(ispin,ispin,iorb,iorb,:)
                   ff(:)  =                                                    - Fdelta(ispin,ispin,iorb,iorb,:)
                   det(:) = fg(:)*conjg(fg(L:1:-1)) + ff(:)*ff(:)
                   F0and(ispin,ispin,iorb,iorb,:) = ff(:)/det(:)
                enddo
             enddo
          end select
       end select
       !
       !
    case ("hybrid")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
       select case(ed_mode)
       case default
          stop "F0and_bath_mats error: called with ed_mode=normal/nonsu2, bath_type=hybrid"
          !
       case ("superc")
          allocate(fgorb(2*Norb,2*Norb),zeta(2*Norb,2*Norb))
          Delta =  delta_bath_array(x,dmft_bath_,axis_)
          Fdelta= fdelta_bath_array(x,dmft_bath_,axis_)
          select case(axis_)
          case default 
             do ispin=1,Nspin
                do i=1,L
                   zeta = zero
                   fgorb= zero
                   do iorb=1,Norb
                      zeta(iorb,iorb)           = x(i) + xmu
                      zeta(iorb+Norb,iorb+Norb) = x(i) - xmu
                   enddo
                   do iorb=1,Norb
                      do jorb=1,Norb
                         fgorb(iorb,jorb)           = zeta(iorb,jorb)           - impHloc(ispin,ispin,iorb,jorb)  - Delta(ispin,ispin,iorb,jorb,i)
                         fgorb(iorb,jorb+Norb)      = zeta(iorb,jorb+Norb)                                        - Fdelta(ispin,ispin,iorb,jorb,i)
                         fgorb(iorb+Norb,jorb)      = zeta(iorb+Norb,jorb)                                        - Fdelta(ispin,ispin,iorb,jorb,i)
                         fgorb(iorb+Norb,jorb+Norb) = zeta(iorb+Norb,jorb+Norb) + impHloc(ispin,ispin,iorb,jorb)  + conjg( Delta(ispin,ispin,iorb,jorb,i) )
                      enddo
                   enddo
                   call inv(fgorb)
                   F0and(ispin,ispin,:,:,i) = fgorb(1:Norb,1+Norb:Norb+Norb)
                enddo
             enddo
          case("real")
             do ispin=1,Nspin
                do i=1,L
                   zeta = zero
                   fgorb= zero
                   do iorb=1,Norb
                      zeta(iorb,iorb)           =        x(i)     + xmu
                      zeta(iorb+Norb,iorb+Norb) = -conjg(x(L-i+1) + xmu)
                   enddo
                   do iorb=1,Norb
                      do jorb=1,Norb
                         fgorb(iorb,jorb)           = zeta(iorb,jorb)           - impHloc(ispin,ispin,iorb,jorb)  - Delta(ispin,ispin,iorb,jorb,i)
                         fgorb(iorb,jorb+Norb)      = zeta(iorb,jorb+Norb)                                        - Fdelta(ispin,ispin,iorb,jorb,i)
                         fgorb(iorb+Norb,jorb)      = zeta(iorb+Norb,jorb)                                        - Fdelta(ispin,ispin,iorb,jorb,i)
                         fgorb(iorb+Norb,jorb+Norb) = zeta(iorb+Norb,jorb+Norb) + impHloc(ispin,ispin,iorb,jorb)  + conjg( Delta(ispin,ispin,iorb,jorb,L-i+1) )
                      enddo
                   enddo
                   call inv(fgorb)
                   F0and(ispin,ispin,:,:,i) = fgorb(1:Norb,1+Norb:Norb+Norb)
                enddo
             enddo
          end select
          deallocate(fgorb,zeta)
       end select
    end select
  end function f0and_bath_array











  function invg0_bath_array(x,dmft_bath_,axis) result(G0and)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and,Delta
    complex(8),dimension(:,:),allocatable               :: fgorb,zeta
    integer                                             :: i,iorb,jorb,ispin,jspin,io,jo,Nso,L
    character(len=*),optional                           :: axis    
    character(len=4)                                    :: axis_
    !
    axis_="mats";if(present(axis))axis_=str(axis)
    !
    G0and = zero
    !
    L=size(x)
    !
    select case(bath_type)
    case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
       !
       select case(ed_mode)
       case default
          Delta = delta_bath_array(x,dmft_bath_)
          do ispin=1,Nspin
             do iorb=1,Norb
                G0and(ispin,ispin,iorb,iorb,:) = x(:) + xmu - impHloc(ispin,ispin,iorb,iorb) - Delta(ispin,ispin,iorb,iorb,:)
             enddo
          enddo
          !
       case ("superc")
          Delta =  delta_bath_array(x,dmft_bath_,axis_)
          select case(axis_)
          case default
             do ispin=1,Nspin
                do iorb=1,Norb
                   G0and(ispin,ispin,iorb,iorb,:)  =  x(:) + xmu - impHloc(ispin,ispin,iorb,iorb) -  Delta(ispin,ispin,iorb,iorb,:)
                enddo
             enddo
          case("real")
             do ispin=1,Nspin
                do iorb=1,Norb
                   G0and(ispin,ispin,iorb,iorb,:) = dreal(x(:)) + xmu - impHloc(ispin,ispin,iorb,iorb) - Delta(ispin,ispin,iorb,iorb,:)
                enddo
             enddo
          end select
          !
       case ("nonsu2")
          Delta = delta_bath_array(x,dmft_bath_)
          allocate(zeta(Nspin,Nspin))
          do i=1,L
             zeta  = (x(i) + xmu)*zeye(Nspin)
             do iorb=1,Norb
                do ispin=1,Nspin
                   do jspin=1,Nspin
                      G0and(ispin,jspin,iorb,iorb,i) = zeta(ispin,jspin) - impHloc(ispin,jspin,iorb,iorb) - Delta(ispin,jspin,iorb,iorb,i)
                   enddo
                enddo
             enddo
          enddo
          deallocate(zeta)
       end select
       !
       !
    case ("hybrid","replica")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
       !
       select case(ed_mode)
       case default
          Delta = delta_bath_array(x,dmft_bath_)
          do ispin=1,Nspin
             do i=1,L
                G0and(ispin,ispin,:,:,i) = (x(i)+xmu)*zeye(Norb)-impHloc(ispin,ispin,:,:)-Delta(ispin,ispin,:,:,i)
             enddo
          enddo
          !
       case ("superc")
          allocate(zeta(2*Norb,2*Norb))
          Delta =  delta_bath_array(x,dmft_bath_,axis_)
          select case(axis_)
          case default
             do ispin=1,Nspin
                do i=1,L
                   zeta = zero
                   do iorb=1,Norb
                      zeta(iorb,iorb)           = x(i) + xmu
                      zeta(iorb+Norb,iorb+Norb) = x(i) - xmu
                   enddo
                   do iorb=1,Norb
                      do jorb=1,Norb
                         G0and(ispin,ispin,iorb,jorb,i) = zeta(iorb,jorb) - impHloc(ispin,ispin,iorb,jorb) - Delta(ispin,ispin,iorb,jorb,i)
                      enddo
                   enddo
                enddo
             enddo
          case("real")
             do ispin=1,Nspin
                do i=1,L
                   zeta = ((x(i))  + xmu)*eye(Norb)
                   do iorb=1,Norb
                      do jorb=1,Norb
                         G0and(ispin,ispin,iorb,jorb,i) = zeta(iorb,jorb) - impHloc(ispin,ispin,iorb,jorb)  - Delta(ispin,ispin,iorb,jorb,i)
                      enddo
                   enddo
                enddo
             enddo
          end select
          deallocate(zeta)
          !
       case ("nonsu2")
          Nso=Nspin*Norb
          allocate(zeta(Nso,Nso))
          Delta = delta_bath_array(x,dmft_bath_)
          do i=1,L
             zeta  = (x(i) + xmu)*zeye(Nso)
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         io = iorb + (ispin-1)*Norb
                         jo = jorb + (jspin-1)*Norb
                         G0and(ispin,jspin,iorb,jorb,i) = zeta(io,jo) - impHloc(ispin,jspin,iorb,jorb) - Delta(ispin,jspin,iorb,jorb,i)
                      enddo
                   enddo
                enddo
             enddo
          enddo
          deallocate(zeta)
       end select
       !
    end select
    !
  end function invg0_bath_array




  !ANOMALous:
  function invf0_bath_array(x,dmft_bath_,axis) result(F0and)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: F0and,Fdelta
    integer                                             :: iorb,jorb,ispin,L
    character(len=*),optional                           :: axis    
    character(len=4)                                    :: axis_
    !
    axis_="mats";if(present(axis))axis_=str(axis)
    !
    F0and=zero
    !
    L = size(x)
    !
    select case(bath_type)
    case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
       !
       select case(ed_mode)
       case default
          stop "Invf0_bath_mats error: called with ed_mode=normal/nonsu2, bath_type=normal"
          !
       case ("superc")
          Fdelta= fdelta_bath_array(x,dmft_bath_,axis_)
          do ispin=1,Nspin
             do iorb=1,Norb
                F0and(ispin,ispin,iorb,iorb,:) = -Fdelta(ispin,ispin,iorb,iorb,:)
             enddo
          enddo
       end select
       !
       !
    case ("hybrid")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
       select case(ed_mode)
       case default
          stop "Invf0_bath_mats error: called with ed_mode=normal/nonsu2, bath_type=hybrid"
          !
       case ("superc")
          Fdelta= fdelta_bath_array(x,dmft_bath_,axis_)
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   F0and(ispin,ispin,iorb,jorb,:) = -Fdelta(ispin,ispin,iorb,jorb,:)
                enddo
             enddo
          enddo
          !
       end select
       !
    end select
  end function invf0_bath_array





END MODULE ED_BATH_FUNCTIONS
