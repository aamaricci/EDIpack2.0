from ctypes import *
import numpy as np
import os,sys
import types

#sigma
def get_sigma(self,ilat=None,ishape=None,axis="m",typ="n"):
    ed_get_sigma_site_n3 = self.library.ed_get_sigma_site_n3
    ed_get_sigma_site_n3.argtypes =[np.ctypeslib.ndpointer(dtype=complex,ndim=3, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                           c_char_p,
                                           c_char_p]       
    ed_get_sigma_site_n3.restype = None

    ed_get_sigma_site_n5 = self.library.ed_get_sigma_site_n5
    ed_get_sigma_site_n5.argtypes =[np.ctypeslib.ndpointer(dtype=complex,ndim=5, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                           c_char_p,
                                           c_char_p]       
    ed_get_sigma_site_n5.restype = None
    
    ed_get_sigma_lattice_n3 = self.library.ed_get_sigma_lattice_n3
    ed_get_sigma_lattice_n3.argtypes =[np.ctypeslib.ndpointer(dtype=complex,ndim=3, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                           c_int,
                                           c_char_p,
                                           c_char_p]       
    ed_get_sigma_lattice_n3.restype = None

    ed_get_sigma_lattice_n4 = self.library.ed_get_sigma_lattice_n4
    ed_get_sigma_lattice_n4.argtypes =[np.ctypeslib.ndpointer(dtype=complex,ndim=4, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                           c_int,
                                           c_char_p,
                                           c_char_p]       
    ed_get_sigma_lattice_n4.restype = None
    
    ed_get_sigma_lattice_n6 = self.library.ed_get_sigma_lattice_n6
    ed_get_sigma_lattice_n4.argtypes =[np.ctypeslib.ndpointer(dtype=complex,ndim=6, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                           c_int,
                                           c_char_p,
                                           c_char_p]       
    ed_get_sigma_lattice_n6.restype = None
    
    norb_aux = c_int.in_dll(self.library, "Norb").value
    nspin_aux = c_int.in_dll(self.library, "Nspin").value
    if axis == "m" or axis == "M":
        nfreq = c_int.in_dll(self.library, "Lmats").value
    elif axis == "r" or axis == "R":
        nfreq = c_int.in_dll(self.library, "Lreal").value
    else:
        raise ValueError("axis can only be 'm' or 'r'")
    if ishape is None:
            ishape = self.dim_hloc + 1
    
    if self.Nineq == 0:
        if ilat is not None:
            raise ValueError("ilat is not defined in single-impurity DMFT")
        if ishape==3:
            Sigma = np.zeros([nspin_aux*norb_aux,nspin_aux*norb_aux,nfreq],dtype=complex,order="F")
            DimSigma = np.asarray([nspin_aux*norb_aux,nspin_aux*norb_aux,nfreq],dtype=np.int64,order="F")
            ed_get_sigma_site_n3(Sigma,DimSigma,c_char_p(axis.encode()),c_char_p(typ.encode()))
        elif ishape==5:
            Sigma = np.zeros([nspin_aux,nspin_aux,norb_aux,norb_aux,nfreq],dtype=complex,order="F")
            DimSigma = np.asarray([nspin_aux,nspin_aux,norb_aux,norb_aux,nfreq],dtype=np.int64,order="F")
            ed_get_sigma_site_n5(Sigma,DimSigma,c_char_p(axis.encode()),c_char_p(typ.encode()))
        else:
            raise ValueError('Shape(array) != 3,5 in get_sigma_site')
        return Sigma
    else:
        if ishape==3:
            Sigma = np.zeros([self.Nineq*nspin_aux*norb_aux,self.Nineq*nspin_aux*norb_aux,nfreq],dtype=complex,order="F")
            DimSigma = np.asarray([self.Nineq*nspin_aux*norb_aux,self.Nineq*nspin_aux*norb_aux,nfreq],dtype=np.int64,order="F")
            ed_get_sigma_site_n3(Sigma,DimSigma,self.Nineq,c_char_p(axis.encode()),c_char_p(typ.encode()))
        elif ishape==4:
            Sigma = np.zeros([self.Nineq,nspin_aux*norb_aux,nspin_aux*norb_aux,nfreq],dtype=complex,order="F")
            DimSigma = np.asarray([self.Nineq,nspin_aux*norb_aux,nspin_aux*norb_aux,nfreq],dtype=np.int64,order="F")
            ed_get_sigma_site_n4(Sigma,DimSigma,self.Nineq,c_char_p(axis.encode()),c_char_p(typ.encode()))
        elif ishape==6:
            Sigma = np.zeros([self.Nineq,nspin_aux,nspin_aux,norb_aux,norb_aux,nfreq],dtype=complex,order="F")
            DimSigma = np.asarray([self.Nineq,nspin_aux,nspin_aux,norb_aux,norb_aux,nfreq],dtype=np.int64,order="F")
            ed_get_sigma_site_n6(Sigma,DimSigma,self.Nineq,c_char_p(axis.encode()),c_char_p(typ.encode()))
        else:
            raise ValueError('Shape(array) != 3,4,6 in get_sigma_lattice')
        if ilat is not None and ishape != 3:
            return Sigma[ilat]
        else:
            return Sigma



#gimp
def get_gimp(self,ilat=None,ishape=None,axis="m",typ="n"):
    ed_get_gimp_site_n3 = self.library.ed_get_gimp_site_n3
    ed_get_gimp_site_n3.argtypes =[np.ctypeslib.ndpointer(dtype=complex,ndim=3, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                           c_char_p,
                                           c_char_p]       
    ed_get_gimp_site_n3.restype = None

    ed_get_gimp_site_n5 = self.library.ed_get_gimp_site_n5
    ed_get_gimp_site_n5.argtypes =[np.ctypeslib.ndpointer(dtype=complex,ndim=5, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                           c_char_p,
                                           c_char_p]       
    ed_get_gimp_site_n5.restype = None
    
    ed_get_gimp_lattice_n3 = self.library.ed_get_gimp_lattice_n3
    ed_get_gimp_lattice_n3.argtypes =[np.ctypeslib.ndpointer(dtype=complex,ndim=3, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                           c_int,
                                           c_char_p,
                                           c_char_p]       
    ed_get_gimp_lattice_n3.restype = None

    ed_get_gimp_lattice_n4 = self.library.ed_get_gimp_lattice_n4
    ed_get_gimp_lattice_n4.argtypes =[np.ctypeslib.ndpointer(dtype=complex,ndim=4, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                           c_int,
                                           c_char_p,
                                           c_char_p]       
    ed_get_gimp_lattice_n4.restype = None
    
    ed_get_gimp_lattice_n6 = self.library.ed_get_gimp_lattice_n6
    ed_get_gimp_lattice_n4.argtypes =[np.ctypeslib.ndpointer(dtype=complex,ndim=6, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                           c_int,
                                           c_char_p,
                                           c_char_p]       
    ed_get_gimp_lattice_n6.restype = None

    norb_aux = c_int.in_dll(self.library, "Norb").value
    nspin_aux = c_int.in_dll(self.library, "Nspin").value
    if axis == "m" or axis == "M":
        nfreq = c_int.in_dll(self.library, "Lmats").value
    elif axis == "r" or axis == "R":
        nfreq = c_int.in_dll(self.library, "Lreal").value
    else:
        raise ValueError("axis can only be 'm' or 'r'") 
    if ishape is None:
        ishape = self.dim_hloc + 1
    
    if self.Nineq == 0:
        if ilat is not None:
            raise ValueError("ilat is not defined in single-impurity DMFT")
        if ishape==3:
            gimp = np.zeros([nspin_aux*norb_aux,nspin_aux*norb_aux,nfreq],dtype=complex,order="F")
            Dimgimp = np.asarray([nspin_aux*norb_aux,nspin_aux*norb_aux,nfreq],dtype=np.int64,order="F")
            ed_get_gimp_site_n3(gimp,Dimgimp,c_char_p(axis.encode()),c_char_p(typ.encode()))
        elif ishape==5:
            gimp = np.zeros([nspin_aux,nspin_aux,norb_aux,norb_aux,nfreq],dtype=complex,order="F")
            Dimgimp = np.asarray([nspin_aux,nspin_aux,norb_aux,norb_aux,nfreq],dtype=np.int64,order="F")
            ed_get_gimp_site_n5(gimp,Dimgimp,c_char_p(axis.encode()),c_char_p(typ.encode()))
        else:
            raise ValueError('Shape(array) != 3,5 in get_gimp_site')
        return gimp
    else:
        if ishape==3:
            gimp = np.zeros([self.Nineq*nspin_aux*norb_aux,self.Nineq*nspin_aux*norb_aux,nfreq],dtype=complex,order="F")
            Dimgimp = np.asarray([self.Nineq*nspin_aux*norb_aux,self.Nineq*nspin_aux*norb_aux,nfreq],dtype=np.int64,order="F")
            ed_get_gimp_site_n3(gimp,Dimgimp,self.Nineq,c_char_p(axis.encode()),c_char_p(typ.encode()))
        elif ishape==4:
            gimp = np.zeros([self.Nineq,nspin_aux*norb_aux,nspin_aux*norb_aux,nfreq],dtype=complex,order="F")
            Dimgimp = np.asarray([self.Nineq,nspin_aux*norb_aux,nspin_aux*norb_aux,nfreq],dtype=np.int64,order="F")
            ed_get_gimp_site_n4(gimp,Dimgimp,self.Nineq,c_char_p(axis.encode()),c_char_p(typ.encode()))
        elif ishape==6:
            Dimgimp = np.zeros([self.Nineq,nspin_aux,nspin_aux,norb_aux,norb_aux,nfreq],dtype=complex,order="F")
            Dimgimp = np.asarray([self.Nineq,nspin_aux,nspin_aux,norb_aux,norb_aux,nfreq],dtype=np.int64,order="F")
            ed_get_gimp_site_n6(gimp,Dimgimp,self.Nineq,c_char_p(axis.encode()),c_char_p(typ.encode()))
        else:
            raise ValueError('Shape(array) != 3,4,6 in get_gimp_lattice')
        if ilat is not None and ishape != 3:
            return gimp[ilat]
        else:
            return gimp
    
    
#observables

#density
def get_dens(self,ilat=None,iorb=None):

    aux_norb = c_int.in_dll(self.library, "Norb").value
    
    ed_get_dens_n1_wrap = self.library.ed_get_dens_n1
    ed_get_dens_n1_wrap.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS')]
    ed_get_dens_n1_wrap.restype = None
    
    ed_get_dens_n2_wrap = self.library.ed_get_dens_n2
    ed_get_dens_n2_wrap.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS'),
                                    c_int]
    ed_get_dens_n2_wrap.restype = None

    if self.Nineq == 0:
        densvec = np.zeros(aux_norb,dtype=float,order="F")
        ed_get_dens_n1_wrap(densvec)
        
        if ilat is not None:
            raise ValueError("ilat cannot be none for single-impurity DMFT")
        elif iorb is not None:
            return densvec[iorb]
        else:
            return densvec
    else:
        densvec = np.zeros([self.Nineq,aux_norb],dtype=float,order="F")
        ed_get_dens_n2_wrap(densvec,self.Nineq)
        
        if ilat is not None and iorb is not None:
            return densvec[ilat,iorb]
        elif ilat is None and iorb is not None:
            return densvec[:,iorb]
        elif ilat is not None and iorb is None:
            return densvec[ilat,:]
        else:
            return densvec
            
#magnetization
def get_mag(self,icomp=None,ilat=None,iorb=None):

    if icomp =="x" or icomp == "X":
        icomp = 0
    elif icomp =="y" or icomp == "Y":
        icomp = 1
    elif icomp == "z" or icomp == "Z":
        icomp = 2

    aux_norb = c_int.in_dll(self.library, "Norb").value
    
    ed_get_mag_n2_wrap = self.library.ed_get_mag_n2
    ed_get_mag_n2_wrap.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS')]
    ed_get_mag_n2_wrap.restype = None
    
    ed_get_mag_n3_wrap = self.library.ed_get_mag_n3
    ed_get_mag_n3_wrap.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=3, flags='F_CONTIGUOUS'),
                                    c_int]
    ed_get_mag_n3_wrap.restype = None

    if self.Nineq == 0:
        magvec = np.zeros([3,aux_norb],dtype=float,order="F")
        ed_get_mag_n2_wrap(magvec)
        
        if ilat is not None:
            raise ValueError("ilat cannot be none for single-impurity DMFT")
        elif iorb is not None and icomp is not None:
            return magvec[icomp,iorb]
        elif iorb is not None and icomp is None:
            return magvec[:,iorb]
        elif iorb is None and icomp is not None:
            return magvec[icomp,:]
        elif iorb is None and icomp is None:
            return magvec
    else:
        magvec = np.zeros([self.Nineq,3,aux_norb],dtype=float,order="F")
        ed_get_mag_n3_wrap(magvec,self.Nineq)
        
        if ilat is not None:
            if iorb is not None and icomp is not None:
                return magvec[ilat,icomp,iorb]
            if iorb is None and icomp is not None:
                return magvec[ilat,icomp,:]
            if iorb is not None and icomp is None:
                return magvec[ilat,:,iorb]
            if iorb is None and icomp is None:
                return magvec[ilat,:,:]
        else:
            if iorb is not None and icomp is not None:
                return magvec[:,icomp,iorb]
            if iorb is None and icomp is not None:
                return magvec[:,icomp,:]
            if iorb is not None and icomp is None:
                return magvec[:,:,iorb]
            if iorb is None and icomp is None:
                return magvec
            
#double occupation
def get_docc(self,ilat=None,iorb=None):

    aux_norb = c_int.in_dll(self.library, "Norb").value
    
    ed_get_docc_n1_wrap = self.library.ed_get_docc_n1
    ed_get_docc_n1_wrap.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS')]
    ed_get_docc_n1_wrap.restype = None
    
    ed_get_docc_n2_wrap = self.library.ed_get_docc_n2
    ed_get_docc_n2_wrap.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS'),
                                    c_int]
    ed_get_docc_n2_wrap.restype = None

    if self.Nineq == 0:
        doccvec = np.zeros(aux_norb,dtype=float,order="F")
        ed_get_docc_n1_wrap(doccvec)
        
        if ilat is not None:
            raise ValueError("ilat cannot be none for single-impurity DMFT")
        elif iorb is not None:
            return doccvec[iorb]
        else:
            return doccvec
    else:
        doccvec = np.zeros([self.Nineq,aux_norb],dtype=float,order="F")
        ed_get_docc_n2_wrap(doccvec,self.Nineq)
        
        if ilat is not None and iorb is not None:
            return doccvec[ilat,iorb]
        elif ilat is None and iorb is not None:
            return doccvec[:,iorb]
        elif ilat is not None and iorb is None:
            return doccvec[ilat,:]
        else:
            return doccvec
            
#energy
def get_eimp(self,ilat=None,ikind=None):
  
    ed_get_eimp_n1_wrap = self.library.ed_get_eimp_n1
    ed_get_eimp_n1_wrap.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS')]
    ed_get_eimp_n1_wrap.restype = None
    
    ed_get_eimp_n2_wrap = self.library.ed_get_eimp_n2
    ed_get_eimp_n2_wrap.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS'),
                                    c_int]
    ed_get_eimp_n2_wrap.restype = None

    if self.Nineq == 0:
        eimp_vec = np.zeros(4,dtype=float,order="F")
        ed_get_eimp_n1_wrap(eimp_vec)
        
        if ilat is not None:
            raise ValueError("ilat cannot be none for single-impurity DMFT")
        elif ikind is not None:
            return eimp_vec[ikind]
        else:
            return eimp_vec
    else:
        eimp_vec = np.zeros([self.Nineq,4],dtype=float,order="F")
        ed_get_eimp_n2_wrap(eimp_vec,self.Nineq)
        
        if ilat is not None and ikind is not None:
            return eimp_vec[ilat,ikind]
        elif ilat is None and ikind is not None:
            return eimp_vec[:,ikind]
        elif ilat is not None and ikind is None:
            return eimp_vec[ilat,:]
        else:
            return eimp_vec
            
#build Sigma
def build_sigma(self,zeta,ilat=None,ishape=None,typ="n"):
    ed_build_sigma_site_n3 = self.library.build_sigma_single_n3
    ed_build_sigma_site_n3.argtypes =[np.ctypeslib.ndpointer(dtype=complex,ndim=1, flags='F_CONTIGUOUS'),
                                      c_int,
                                      np.ctypeslib.ndpointer(dtype=complex,ndim=3, flags='F_CONTIGUOUS'),
                                      np.ctypeslib.ndpointer(dtype=complex,ndim=3, flags='F_CONTIGUOUS')
                                      ]       
    ed_build_sigma_site_n3.restype = None

    ed_build_sigma_site_n5 = self.library.build_sigma_single_n5
    ed_build_sigma_site_n5.argtypes =[np.ctypeslib.ndpointer(dtype=complex,ndim=1, flags='F_CONTIGUOUS'),
                                      c_int,
                                      np.ctypeslib.ndpointer(dtype=complex,ndim=5, flags='F_CONTIGUOUS'),
                                      np.ctypeslib.ndpointer(dtype=complex,ndim=5, flags='F_CONTIGUOUS')
                                      ]       
    ed_build_sigma_site_n5.restype = None
    
    ed_build_sigma_lattice_n3 = self.library.build_sigma_ineq_n3
    ed_build_sigma_lattice_n3.argtypes =[np.ctypeslib.ndpointer(dtype=complex,ndim=1, flags='F_CONTIGUOUS'),
                                         c_int,
                                         c_int,
                                         np.ctypeslib.ndpointer(dtype=complex,ndim=3, flags='F_CONTIGUOUS'),
                                         np.ctypeslib.ndpointer(dtype=complex,ndim=3, flags='F_CONTIGUOUS')
                                         ]       
    ed_build_sigma_lattice_n3.restype = None

    ed_build_sigma_lattice_n4 = self.library.build_sigma_ineq_n4
    ed_build_sigma_lattice_n4.argtypes =[np.ctypeslib.ndpointer(dtype=complex,ndim=1, flags='F_CONTIGUOUS'),
                                         c_int,
                                         c_int,
                                         np.ctypeslib.ndpointer(dtype=complex,ndim=4, flags='F_CONTIGUOUS'),
                                         np.ctypeslib.ndpointer(dtype=complex,ndim=4, flags='F_CONTIGUOUS')
                                         ]     
    ed_build_sigma_lattice_n4.restype = None
    
    ed_build_sigma_lattice_n6 = self.library.build_sigma_ineq_n6
    ed_build_sigma_lattice_n4.argtypes =[np.ctypeslib.ndpointer(dtype=complex,ndim=1, flags='F_CONTIGUOUS'),
                                         c_int,
                                         c_int,
                                         np.ctypeslib.ndpointer(dtype=complex,ndim=6, flags='F_CONTIGUOUS'),
                                         np.ctypeslib.ndpointer(dtype=complex,ndim=6, flags='F_CONTIGUOUS')
                                         ]        
    ed_build_sigma_lattice_n6.restype = None
    
    norb_aux = c_int.in_dll(self.library, "Norb").value
    nspin_aux = c_int.in_dll(self.library, "Nspin").value
    
    if np.isscalar(zeta):
        zeta=[zeta]
    zeta = np.asarray(zeta,order='F')
    nfreq = np.shape(zeta)[0]
    
    if ishape is None:
            ishape = self.dim_hloc + 1
    
    if self.Nineq == 0:
        if ilat is not None:
            raise ValueError("ilat is not defined in single-impurity DMFT")
        if ishape==3:
            Sigma_normal = np.zeros([nspin_aux*norb_aux,nspin_aux*norb_aux,nfreq],dtype=complex,order="F")
            Sigma_anomalous = np.zeros_like(Sigma_normal)
            ed_build_sigma_site_n3(zeta,nfreq,Sigma_normal,Sigma_anomalous)
        elif ishape==5:
            Sigma_normal = np.zeros([nspin_aux,nspin_aux,norb_aux,norb_aux,nfreq],dtype=complex,order="F")
            Sigma_anomalous = np.zeros_like(Sigma_normal)
            ed_build_sigma_site_n5(zeta,nfreq,Sigma_normal,Sigma_anomalous)
        else:
            raise ValueError('Shape(array) != 3,5 in build_sigma_site')
        if typ=="n":
            return Sigma_normal
        else:
            return Sigma_anomalous
    else:
        if ishape==3:
            Sigma_normal = np.zeros([self.Nineq*nspin_aux*norb_aux,self.Nineq*nspin_aux*norb_aux,nfreq],dtype=complex,order="F")
            Sigma_anomalous = np.zeros_like(Sigma_normal)
            ed_build_sigma_site_n3(zeta,nfreq,self.Nineq,Sigma_normal,Sigma_anomalous)
        elif ishape==4:
            Sigma_normal = np.zeros([self.Nineq,nspin_aux*norb_aux,nspin_aux*norb_aux,nfreq],dtype=complex,order="F")
            Sigma_anomalous = np.zeros_like(Sigma_normal)
            ed_build_sigma_site_n4(zeta,nfreq,self.Nineq,Sigma_normal,Sigma_anomalous)
        elif ishape==6:
            Sigma_normal = np.zeros([self.Nineq,nspin_aux,nspin_aux,norb_aux,norb_aux,nfreq],dtype=complex,order="F")
            Sigma_anomalous = np.zeros_like(Sigma_normal)
            ed_build_sigma_site_n6(zeta,nfreq,self.Nineq,Sigma_normal,Sigma_anomalous)
        else:
            raise ValueError('Shape(array) != 3,4,6 in build_sigma_lattice')
        if ilat is not None and ishape != 3:
            if typ == "n":
                return Sigma_normal[ilat]
            else:
                return Sigma_anomalous[ilat]
        else:
            if typ == "n":
                return Sigma_normal
            else:
                return Sigma_anomalous

#build gimp
def build_gimp(self,zeta,ilat=None,ishape=None,typ="n"):
    ed_build_gimp_site_n3 = self.library.build_gimp_single_n3
    ed_build_gimp_site_n3.argtypes =[np.ctypeslib.ndpointer(dtype=complex,ndim=1, flags='F_CONTIGUOUS'),
                                      c_int,
                                      np.ctypeslib.ndpointer(dtype=complex,ndim=3, flags='F_CONTIGUOUS'),
                                      np.ctypeslib.ndpointer(dtype=complex,ndim=3, flags='F_CONTIGUOUS')
                                      ]       
    ed_build_gimp_site_n3.restype = None

    ed_build_gimp_site_n5 = self.library.build_gimp_single_n5
    ed_build_gimp_site_n5.argtypes =[np.ctypeslib.ndpointer(dtype=complex,ndim=1, flags='F_CONTIGUOUS'),
                                      c_int,
                                      np.ctypeslib.ndpointer(dtype=complex,ndim=5, flags='F_CONTIGUOUS'),
                                      np.ctypeslib.ndpointer(dtype=complex,ndim=5, flags='F_CONTIGUOUS')
                                      ]       
    ed_build_gimp_site_n5.restype = None
    
    ed_build_gimp_lattice_n3 = self.library.build_gimp_ineq_n3
    ed_build_gimp_lattice_n3.argtypes =[np.ctypeslib.ndpointer(dtype=complex,ndim=1, flags='F_CONTIGUOUS'),
                                         c_int,
                                         c_int,
                                         np.ctypeslib.ndpointer(dtype=complex,ndim=3, flags='F_CONTIGUOUS'),
                                         np.ctypeslib.ndpointer(dtype=complex,ndim=3, flags='F_CONTIGUOUS')
                                         ]       
    ed_build_gimp_lattice_n3.restype = None

    ed_build_gimp_lattice_n4 = self.library.build_gimp_ineq_n4
    ed_build_gimp_lattice_n4.argtypes =[np.ctypeslib.ndpointer(dtype=complex,ndim=1, flags='F_CONTIGUOUS'),
                                         c_int,
                                         c_int,
                                         np.ctypeslib.ndpointer(dtype=complex,ndim=4, flags='F_CONTIGUOUS'),
                                         np.ctypeslib.ndpointer(dtype=complex,ndim=4, flags='F_CONTIGUOUS')
                                         ]     
    ed_build_gimp_lattice_n4.restype = None
    
    ed_build_gimp_lattice_n6 = self.library.build_gimp_ineq_n6
    ed_build_gimp_lattice_n4.argtypes =[np.ctypeslib.ndpointer(dtype=complex,ndim=1, flags='F_CONTIGUOUS'),
                                         c_int,
                                         c_int,
                                         np.ctypeslib.ndpointer(dtype=complex,ndim=6, flags='F_CONTIGUOUS'),
                                         np.ctypeslib.ndpointer(dtype=complex,ndim=6, flags='F_CONTIGUOUS')
                                         ]        
    ed_build_gimp_lattice_n6.restype = None
    
    norb_aux = c_int.in_dll(self.library, "Norb").value
    nspin_aux = c_int.in_dll(self.library, "Nspin").value
    
    if np.isscalar(zeta):
        zeta=[zeta]
    zeta = np.asarray(zeta,order='F')
    nfreq = np.shape(zeta)[0]
    
    if ishape is None:
            ishape = self.dim_hloc + 1
    
    if self.Nineq == 0:
        if ilat is not None:
            raise ValueError("ilat is not defined in single-impurity DMFT")
        if ishape==3:
            gimp_normal = np.zeros([nspin_aux*norb_aux,nspin_aux*norb_aux,nfreq],dtype=complex,order="F")
            gimp_anomalous = np.zeros_like(gimp_normal)
            ed_build_gimp_site_n3(zeta,nfreq,gimp_normal,gimp_anomalous)
        elif ishape==5:
            gimp_normal = np.zeros([nspin_aux,nspin_aux,norb_aux,norb_aux,nfreq],dtype=complex,order="F")
            gimp_anomalous = np.zeros_like(gimp_normal)
            ed_build_gimp_site_n5(zeta,nfreq,gimp_normal,gimp_anomalous)
        else:
            raise ValueError('Shape(array) != 3,5 in build_gimp_site')
        if typ=="n":
            return gimp_normal
        else:
            return gimp_anomalous
    else:
        if ishape==3:
            gimp_normal = np.zeros([self.Nineq*nspin_aux*norb_aux,self.Nineq*nspin_aux*norb_aux,nfreq],dtype=complex,order="F")
            gimp_anomalous = np.zeros_like(gimp_normal)
            ed_build_gimp_site_n3(zeta,nfreq,self.Nineq,gimp_normal,gimp_anomalous)
        elif ishape==4:
            gimp_normal = np.zeros([self.Nineq,nspin_aux*norb_aux,nspin_aux*norb_aux,nfreq],dtype=complex,order="F")
            gimp_anomalous = np.zeros_like(gimp_normal)
            ed_build_gimp_site_n4(zeta,nfreq,self.Nineq,gimp_normal,gimp_anomalous)
        elif ishape==6:
            gimp_normal = np.zeros([self.Nineq,nspin_aux,nspin_aux,norb_aux,norb_aux,nfreq],dtype=complex,order="F")
            gimp_anomalous = np.zeros_like(gimp_normal)
            ed_build_gimp_site_n6(zeta,nfreq,self.Nineq,gimp_normal,gimp_anomalous)
        else:
            raise ValueError('Shape(array) != 3,4,6 in build_gimp_lattice')
        if ilat is not None and ishape != 3:
            if typ == "n":
                return gimp_normal[ilat]
            else:
                return gimp_anomalous[ilat]
        else:
            if typ == "n":
                return gimp_normal
            else:
                return gimp_anomalous