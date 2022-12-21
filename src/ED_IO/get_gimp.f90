subroutine ed_get_gimp_site(self,axis,type)
  complex(8),dimension(..),intent(inout) :: self
  character(len=*),optional              :: axis
  character(len=*),optional              :: type
  character(len=1)                       :: axis_
  character(len=1)                       :: type_
  axis_='m';if(present(axis))axis_=trim(axis)
  type_='n';if(present(type))type_=trim(type)
  select case(type_)
  case default; stop "ed_get_gimp ERROR: type is neither Normal, nor Anomalous"
  case ('n','N')
     select case(axis_)
     case default;stop "ed_get_gimp ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        select rank(self)
        rank default; stop "ed_get_gimp ERROR: self has a wrong rank"
        rank (3)
        call assert_shape(self,[Nspin*Norb,Nspin*Norb,Lmats],'ed_get_gimp','self')
        self = nn2so_reshape(impGmats,Nspin,Norb,Lmats)
        rank (5)
        call assert_shape(self,[Nspin,Nspin,Norb,Norb,Lmats],'ed_get_gimp','self')
        self = impGmats
        end select
     case('r','R')
        select rank(self)
        rank default; stop "ed_get_gimp ERROR: self has a wrong rank"
        rank (3)
        call assert_shape(self,[Nspin*Norb,Nspin*Norb,Lreal],'ed_get_gimp','self')
        self = nn2so_reshape(impGreal,Nspin,Norb,Lmats)
        rank (5)
        call assert_shape(self,[Nspin,Nspin,Norb,Norb,Lreal],'ed_get_gimp','self')
        self = impGreal
        end select
     end select
  case('a','A')
     select case(axis_)
     case default;stop "ed_get_gimp ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        select rank(self)
        rank default; stop "ed_get_gimp ERROR: self has a wrong rank"
        rank (3)
        call assert_shape(self,[Nspin*Norb,Nspin*Norb,Lmats],'ed_get_gimp','self')
        self = nn2so_reshape(impFmats,Nspin,Norb,Lmats)
        rank (5)
        call assert_shape(self,[Nspin,Nspin,Norb,Norb,Lmats],'ed_get_gimp','self')
        self = impFmats
        end select
     case('r','R')
        select rank(self)
        rank default; stop "ed_get_gimp ERROR: self has a wrong rank"
        rank (3)
        call assert_shape(self,[Nspin*Norb,Nspin*Norb,Lreal],'ed_get_gimp','self')
        self = nn2so_reshape(impFreal,Nspin,Norb,Lmats)
        rank (5)
        call assert_shape(self,[Nspin,Nspin,Norb,Norb,Lreal],'ed_get_gimp','self')
        self = impFreal
        end select
     end select
  end select
end subroutine ed_get_gimp_site



subroutine ed_get_gimp_lattice(self,nlat,axis,type)
  complex(8),dimension(..),intent(inout) :: self
  integer,intent(in)                     :: nlat
  character(len=*),optional              :: axis
  character(len=*),optional              :: type
  character(len=1)                       :: axis_
  character(len=1)                       :: type_
  integer                                :: ilat
  axis_='m';if(present(axis))axis_=trim(axis)
  type_='n';if(present(type))type_=trim(type)
  if(Nlat/=size(Gmats_ineq,1))stop "ERROR ed_get_gimp: wrong Nlat"
  select case(type_)
  case default; stop "ed_get_gimp ERROR: type is neither Normal, nor Anomalous"
  case ('n','N')
     select case(axis_)
     case default;stop "ed_get_gimp ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        select rank(self)
        rank default; stop "ed_get_gimp ERROR: self has a wrong rank"
        rank (3)
        call assert_shape(self,[Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lmats],'ed_get_gimp','self')
        self = nnn2lso_reshape(Gmats_ineq,Nlat,Nspin,Norb,Lmats)
        rank (4)
        call assert_shape(self,[Nlat,Nspin*Norb,Nspin*Norb,Lmats],'ed_get_gimp','self')
        do ilat=1,Nlat
           self(ilat,:,:,:) = nn2so_reshape(Gmats_ineq(ilat,:,:,:,:,:),Nspin,Norb,Lmats)
        enddo
        rank (6)
        call assert_shape(self,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],'ed_get_gimp','self')
        self = Gmats_ineq
        end select
     case('r','R')
        select rank(self)
        rank default; stop "ed_get_gimp ERROR: self has a wrong rank"
        rank (3)
        call assert_shape(self,[Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lreal],'ed_get_gimp','self')
        self = nnn2lso_reshape(Freal_ineq,Nlat,Nspin,Norb,Lreal)
        rank (4)
        call assert_shape(self,[Nlat,Nspin*Norb,Nspin*Norb,Lreal],'ed_get_gimp','self')
        do ilat=1,Nlat
           self(ilat,:,:,:) = nn2so_reshape(Freal_ineq(ilat,:,:,:,:,:),Nspin,Norb,Lreal)
        enddo
        rank (6)
        call assert_shape(self,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],'ed_get_gimp','self')
        self = Freal_ineq
        end select
     end select
  case('a','A')
     select case(axis_)
     case default;stop "ed_get_gimp ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        select rank(self)
        rank default; stop "ed_get_gimp ERROR: self has a wrong rank"
        rank (3)
        call assert_shape(self,[Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lmats],'ed_get_gimp','self')
        self = nnn2lso_reshape(Fmats_ineq,Nlat,Nspin,Norb,Lmats)
        rank (4)
        call assert_shape(self,[Nlat,Nspin*Norb,Nspin*Norb,Lmats],'ed_get_gimp','self')
        do ilat=1,Nlat
           self(ilat,:,:,:) = nn2so_reshape(Fmats_ineq(ilat,:,:,:,:,:),Nspin,Norb,Lmats)
        enddo
        rank (6)
        call assert_shape(self,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],'ed_get_gimp','self')
        self = Fmats_ineq
        end select
     case('r','R')
        select rank(self)
        rank default; stop "ed_get_gimp ERROR: self has a wrong rank"
        rank (3)
        call assert_shape(self,[Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lreal],'ed_get_gimp','self')
        self = nnn2lso_reshape(Freal_ineq,Nlat,Nspin,Norb,Lreal)
        rank (4)
        call assert_shape(self,[Nlat,Nspin*Norb,Nspin*Norb,Lreal],'ed_get_gimp','self')
        do ilat=1,Nlat
           self(ilat,:,:,:) = nn2so_reshape(Freal_ineq(ilat,:,:,:,:,:),Nspin,Norb,Lreal)
        enddo
        rank (6)
        call assert_shape(self,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],'ed_get_gimp','self')
        self = Freal_ineq
        end select
     end select
  end select
end subroutine ed_get_gimp_lattice



