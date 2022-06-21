[1mdiff --git a/src/ED_BATH/fitgf_replica.f90 b/src/ED_BATH/fitgf_replica.f90[m
[1mindex c0c5319..00b1a14 100644[m
[1m--- a/src/ED_BATH/fitgf_replica.f90[m
[1m+++ b/src/ED_BATH/fitgf_replica.f90[m
[36m@@ -98,6 +98,7 @@[m [msubroutine chi2_fitgf_replica(fg,bath_)[m
                 ftol=cg_Ftol,  &[m
                 istop=cg_stop, &[m
                 iverbose=(ed_verbose>3))[m
[32m+[m[32m           print*,"done using"[m
         case ("delta")[m
            call fmin_cg(array_bath,chi2_delta_replica,grad_chi2_delta_replica,&[m
                 iter,chi,&[m
[36m@@ -117,6 +118,7 @@[m [msubroutine chi2_fitgf_replica(fg,bath_)[m
                 ftol=cg_Ftol,  &[m
                 istop=cg_stop, &[m
                 iverbose=(ed_verbose>3))[m
[32m+[m[32m           print*, "done using"[m
         case ("delta")[m
            call fmin_cg(array_bath,chi2_delta_replica,&[m
                 iter,chi,&[m
[1mdiff --git a/src/ED_NORMAL/ED_CHI_EXCT.f90 b/src/ED_NORMAL/ED_CHI_EXCT.f90[m
[1mindex f30d124..b276406 100644[m
[1m--- a/src/ED_NORMAL/ED_CHI_EXCT.f90[m
[1m+++ b/src/ED_NORMAL/ED_CHI_EXCT.f90[m
[36m@@ -235,7 +235,7 @@[m [mcontains[m
              do k=1,sectorK%Dim[m
                 call apply_op_CDG(k,i,sgn,ipos,ialfa,2,sectorK,sectorI)[m
                 if(sgn==0.OR.k==0)cycle[m
[31m-                vvinit(i) = sgn*vvinit_tmp(k)[m
[32m+[m[32m                vvinit(i) = -sgn*vvinit_tmp(k)[m[41m [m
              enddo[m
              deallocate(vvinit_tmp)[m
              call delete_sector(sectorK)[m
[36m@@ -256,7 +256,7 @@[m [mcontains[m
              do k=1,sectorK%Dim[m
                 call apply_op_CDG(k,i,sgn,ipos,ialfa,1,sectorK,sectorI)[m
                 if(sgn==0.OR.k==0)cycle[m
[31m-                vvinit(i) =  sgn*vvinit_tmp(k) - vvinit(i)[m
[32m+[m[32m                vvinit(i) =  sgn*vvinit_tmp(k) + vvinit(i)[m
              enddo[m
              deallocate(vvinit_tmp)[m
              call delete_sector(sectorK)[m
[1mdiff --git a/src/ED_VERSION.f90 b/src/ED_VERSION.f90[m
[1mindex 2a91a69..ca31b61 100644[m
[1m--- a/src/ED_VERSION.f90[m
[1m+++ b/src/ED_VERSION.f90[m
[36m@@ -1,5 +1,5 @@[m
 MODULE ED_VERSION[m
   implicit none[m
   !GIT VERSION[m
[31m-  character(len=41),parameter,public :: version = "7bea97ad297ff4dd5eea0bc9714f12891f6771f6"[m
[32m+[m[32m  character(len=41),parameter,public :: version = "b8a019ce1dd9fcb17b622b9e8e49480d2f633ada"[m
 END MODULE ED_VERSION[m
