from ctypes import *
import numpy as np
import os,sys
import types

#sigma
def get_sigma(self,Sigma,Nlat=-1,axis="m",typ="n"):
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
    
    DimSigma = np.asarray(np.shape(Sigma),dtype=np.int64,order="F")
    if Nlat < 0:
        if len(DimSigma)==3:
            ed_get_sigma_site_n3(Sigma,DimSigma,c_char_p(axis.encode()),c_char_p(typ.encode()))
        if len(DimSigma)==5:
            ed_get_sigma_site_n5(Sigma,DimSigma,c_char_p(axis.encode()),c_char_p(typ.encode()))
        else:
            raise ValueError('Shape(array) != 3 in get_bath_component')
    else:
        if len(DimSigma)==3:
            ed_get_sigma_site_n3(Sigma,DimSigma,Nlat,c_char_p(axis.encode()),c_char_p(typ.encode()))
        if len(DimSigma)==4:
            ed_get_sigma_site_n4(Sigma,DimSigma,Nlat,c_char_p(axis.encode()),c_char_p(typ.encode()))
        if len(DimSigma)==6:
            ed_get_sigma_site_n6(Sigma,DimSigma,Nlat,c_char_p(axis.encode()),c_char_p(typ.encode()))
        else:
            raise ValueError('Shape(array) != 3 in get_bath_component')
    return Sigma


#gimp
def get_gimp(self,gimp,Nlat=-1,axis="m",typ="n"):
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
    
    Dimgimp = np.asarray(np.shape(gimp),dtype=np.int64,order="F")
    if Nlat < 0:
        if len(Dimgimp)==3:
            ed_get_gimp_site_n3(gimp,Dimgimp,c_char_p(axis.encode()),c_char_p(typ.encode()))
        if len(Dimgimp)==5:
            ed_get_gimp_site_n5(gimp,Dimgimp,c_char_p(axis.encode()),c_char_p(typ.encode()))
        else:
            raise ValueError('Shape(array) != 3 in get_bath_component')
    else:
        if len(Dimgimp)==3:
            ed_get_gimp_site_n3(gimp,Dimgimp,Nlat,c_char_p(axis.encode()),c_char_p(typ.encode()))
        if len(Dimgimp)==4:
            ed_get_gimp_site_n4(gimp,Dimgimp,Nlat,c_char_p(axis.encode()),c_char_p(typ.encode()))
        if len(Dimgimp)==6:
            ed_get_gimp_site_n6(gimp,Dimgimp,Nlat,c_char_p(axis.encode()),c_char_p(typ.encode()))
        else:
            raise ValueError('Shape(array) != 3 in get_bath_component')
    return gimp