subroutine ed_get_sigma_site(self,axis,type)
  complex(8),dimension(..),intent(inout) :: self
  character(len=*),optional              :: axis
  character(len=*),optional              :: type
  character(len=1)                       :: axis_
  character(len=1)                       :: type_
  axis_='m';if(present(axis))axis_=trim(axis)
  type_='n';if(present(type))type_=trim(type)
  select case(type_)
  case default; stop "ed_get_sigma ERROR: type is neither Normal, nor Anomalous"
  case ('n','N')
     select case(axis_)
     case default;stop "ed_get_sigma ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        select rank(self)
        rank default; stop "ed_get_sigma ERROR: self has a wrong rank"
        rank (3)
        call assert_shape(self,[Nspin*Norb,Nspin*Norb,Lmats],'ed_get_sigma','self')
        self = nn2so_reshape(impSmats,Nspin,Norb,Lmats)
        rank (5)
        call assert_shape(self,[Nspin,Nspin,Norb,Norb,Lmats],'ed_get_sigma','self')
        self = impSmats
        end select
     case('r','R')
        select rank(self)
        rank default; stop "ed_get_sigma ERROR: self has a wrong rank"
        rank (3)
        call assert_shape(self,[Nspin*Norb,Nspin*Norb,Lreal],'ed_get_sigma','self')
        self = nn2so_reshape(impSreal,Nspin,Norb,Lmats)
        rank (5)
        call assert_shape(self,[Nspin,Nspin,Norb,Norb,Lreal],'ed_get_sigma','self')
        self = impSreal
        end select
     end select
  case('a','A')
     select case(axis_)
     case default;stop "ed_get_sigma ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        select rank(self)
        rank default; stop "ed_get_sigma ERROR: self has a wrong rank"
        rank (3)
        call assert_shape(self,[Nspin*Norb,Nspin*Norb,Lmats],'ed_get_sigma','self')
        self = nn2so_reshape(impSAmats,Nspin,Norb,Lmats)
        rank (5)
        call assert_shape(self,[Nspin,Nspin,Norb,Norb,Lmats],'ed_get_sigma','self')
        self = impSmats
        end select
     case('r','R')
        select rank(self)
        rank default; stop "ed_get_sigma ERROR: self has a wrong rank"
        rank (3)
        call assert_shape(self,[Nspin*Norb,Nspin*Norb,Lreal],'ed_get_sigma','self')
        self = nn2so_reshape(impSAreal,Nspin,Norb,Lmats)
        rank (5)
        call assert_shape(self,[Nspin,Nspin,Norb,Norb,Lreal],'ed_get_sigma','self')
        self = impSAreal
        end select
     end select
  end select
end subroutine ed_get_sigma_site



subroutine ed_get_sigma_lattice(self,nlat,axis,type)
  complex(8),dimension(..),intent(inout) :: self
  integer,intent(in)                     :: nlat
  character(len=*),optional              :: axis
  character(len=*),optional              :: type
  character(len=1)                       :: axis_
  character(len=1)                       :: type_
  integer                                :: ilat
  axis_='m';if(present(axis))axis_=trim(axis)
  type_='n';if(present(type))type_=trim(type)
  if(Nlat/=size(Smats_ineq,1))stop "ERROR ed_get_sigma: wrong Nlat"
  select case(type_)
  case default; stop "ed_get_sigma ERROR: type is neither Normal, nor Anomalous"
  case ('n','N')
     select case(axis_)
     case default;stop "ed_get_sigma ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        select rank(self)
        rank default; stop "ed_get_sigma ERROR: self has a wrong rank"
        rank (3)
        call assert_shape(self,[Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lmats],'ed_get_sigma','self')
        self = nnn2lso_reshape(Smats_ineq,Nlat,Nspin,Norb,Lmats)
        rank (4)
        call assert_shape(self,[Nlat,Nspin*Norb,Nspin*Norb,Lmats],'ed_get_sigma','self')
        do ilat=1,Nlat
           self(ilat,:,:,:) = nn2so_reshape(Smats_ineq(ilat,:,:,:,:,:),Nspin,Norb,Lmats)
        enddo
        rank (6)
        call assert_shape(self,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],'ed_get_sigma','self')
        self = Smats_ineq
        end select
     case('r','R')
        select rank(self)
        rank default; stop "ed_get_sigma ERROR: self has a wrong rank"
        rank (3)
        call assert_shape(self,[Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lreal],'ed_get_sigma','self')
        self = nnn2lso_reshape(Sreal_ineq,Nlat,Nspin,Norb,Lreal)
        rank (4)
        call assert_shape(self,[Nlat,Nspin*Norb,Nspin*Norb,Lreal],'ed_get_sigma','self')
        do ilat=1,Nlat
           self(ilat,:,:,:) = nn2so_reshape(Sreal_ineq(ilat,:,:,:,:,:),Nspin,Norb,Lreal)
        enddo
        rank (6)
        call assert_shape(self,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],'ed_get_sigma','self')
        self = Sreal_ineq
        end select
     end select
  case('a','A')
     select case(axis_)
     case default;stop "ed_get_sigma ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        select rank(self)
        rank default; stop "ed_get_sigma ERROR: self has a wrong rank"
        rank (3)
        call assert_shape(self,[Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lmats],'ed_get_sigma','self')
        self = nnn2lso_reshape(SAmats_ineq,Nlat,Nspin,Norb,Lmats)
        rank (4)
        call assert_shape(self,[Nlat,Nspin*Norb,Nspin*Norb,Lmats],'ed_get_sigma','self')
        do ilat=1,Nlat
           self(ilat,:,:,:) = nn2so_reshape(SAmats_ineq(ilat,:,:,:,:,:),Nspin,Norb,Lmats)
        enddo
        rank (6)
        call assert_shape(self,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],'ed_get_sigma','self')
        self = SAmats_ineq
        end select
     case('r','R')
        select rank(self)
        rank default; stop "ed_get_sigma ERROR: self has a wrong rank"
        rank (3)
        call assert_shape(self,[Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lreal],'ed_get_sigma','self')
        self = nnn2lso_reshape(SAreal_ineq,Nlat,Nspin,Norb,Lreal)
        rank (4)
        call assert_shape(self,[Nlat,Nspin*Norb,Nspin*Norb,Lreal],'ed_get_sigma','self')
        do ilat=1,Nlat
           self(ilat,:,:,:) = nn2so_reshape(SAreal_ineq(ilat,:,:,:,:,:),Nspin,Norb,Lreal)
        enddo
        rank (6)
        call assert_shape(self,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],'ed_get_sigma','self')
        self = SAreal_ineq
        end select
     end select
  end select
end subroutine ed_get_sigma_lattice



