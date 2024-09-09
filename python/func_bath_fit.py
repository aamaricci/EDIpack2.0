from ctypes import *
import numpy as np
import os,sys
import types

def chi2_fitgf(self,*args,ispin=0,iorb=-1,fmpi=True):
#single normal
    chi2_fitgf_single_normal_n3 = self.library.chi2_fitgf_single_normal_n3
    chi2_fitgf_single_normal_n3.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=3, flags='F_CONTIGUOUS'),
                                            np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                            np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),
                                            np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                            c_int,
                                            c_int,
                                            c_int] 
    chi2_fitgf_single_normal_n3.restype = None

    chi2_fitgf_single_normal_n5 = self.library.chi2_fitgf_single_normal_n5
    chi2_fitgf_single_normal_n5.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=5, flags='F_CONTIGUOUS'),
                                            np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                            np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),
                                            np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                            c_int,
                                            c_int,
                                            c_int] 
    chi2_fitgf_single_normal_n5.restype = None

    #single superc
    chi2_fitgf_single_superc_n3 = self.library.chi2_fitgf_single_superc_n3
    chi2_fitgf_single_superc_n3.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=3, flags='F_CONTIGUOUS'),
                                            np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                            np.ctypeslib.ndpointer(dtype=complex,ndim=3, flags='F_CONTIGUOUS'),
                                            np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                            np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),
                                            np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                            c_int,
                                            c_int,
                                            c_int] 
    chi2_fitgf_single_superc_n3.restype = None

    chi2_fitgf_single_superc_n5 = self.library.chi2_fitgf_single_superc_n5
    chi2_fitgf_single_superc_n5.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=5, flags='F_CONTIGUOUS'),
                                            np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                            np.ctypeslib.ndpointer(dtype=complex,ndim=5, flags='F_CONTIGUOUS'),
                                            np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                            np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),
                                            np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                            c_int,
                                            c_int,
                                            c_int] 
    chi2_fitgf_single_superc_n5.restype = None

    #lattice normal
    chi2_fitgf_lattice_normal_n3 = self.library.chi2_fitgf_lattice_normal_n3
    chi2_fitgf_lattice_normal_n3.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=3, flags='F_CONTIGUOUS'),
                                             np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                             np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS'),
                                             np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                             c_int] 
    chi2_fitgf_lattice_normal_n3.restype = None

    chi2_fitgf_lattice_normal_n4 = self.library.chi2_fitgf_lattice_normal_n4
    chi2_fitgf_lattice_normal_n4.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=4, flags='F_CONTIGUOUS'),
                                             np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                             np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS'),
                                             np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                             c_int] 
    chi2_fitgf_lattice_normal_n4.restype = None

    chi2_fitgf_lattice_normal_n6 = self.library.chi2_fitgf_lattice_normal_n6
    chi2_fitgf_lattice_normal_n6.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=6, flags='F_CONTIGUOUS'),
                                             np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                             np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS'),
                                             np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                             c_int] 
    chi2_fitgf_lattice_normal_n6.restype = None

    #lattice superc
    chi2_fitgf_lattice_superc_n3 = self.library.chi2_fitgf_lattice_superc_n3
    chi2_fitgf_lattice_superc_n3.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=3, flags='F_CONTIGUOUS'),
                                             np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                             np.ctypeslib.ndpointer(dtype=complex,ndim=3, flags='F_CONTIGUOUS'),
                                             np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                             np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS'),
                                             np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                             c_int] 
    chi2_fitgf_lattice_superc_n3.restype = None

    chi2_fitgf_lattice_superc_n4 = self.library.chi2_fitgf_lattice_superc_n4
    chi2_fitgf_lattice_superc_n4.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=4, flags='F_CONTIGUOUS'),
                                             np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                             np.ctypeslib.ndpointer(dtype=complex,ndim=4, flags='F_CONTIGUOUS'),
                                             np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                             np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS'),
                                             np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                             c_int] 
    chi2_fitgf_lattice_superc_n4.restype = None

    chi2_fitgf_lattice_superc_n6 = self.library.chi2_fitgf_lattice_superc_n6
    chi2_fitgf_lattice_superc_n6.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=6, flags='F_CONTIGUOUS'),
                                             np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                             np.ctypeslib.ndpointer(dtype=complex,ndim=6, flags='F_CONTIGUOUS'),
                                             np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                             np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS'),
                                             np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                             c_int] 
    chi2_fitgf_lattice_superc_n6.restype = None
    
    #main function
    ispin = ispin + 1
    iorb = iorb + 1
    if len(args) == 2: #normal
        g = args[0]
        bath = args[1]
        dim_g = np.asarray(np.shape(g),dtype=np.int64,order="F")
        dim_bath = np.asarray(np.shape(bath),dtype=np.int64,order="F")
        if len(dim_bath) == 1: #single
            if len(dim_g) == 3:
                chi2_fitgf_single_normal_n3(g,dim_g,bath,dim_bath,ispin,iorb,fmpi)
            elif len(dim_g) == 5:
                chi2_fitgf_single_normal_n5(g,dim_g,bath,dim_bath,ispin,iorb,fmpi)
            else:
                raise ValueError("chi_fitgf_normal: takes dim(g) = 3 or 5")
        elif len(dim_bath) == 2: #lattice
            if len(dim_g) == 3:
                chi2_fitgf_lattice_normal_n3(g,dim_g,bath,dim_bath,ispin)
            if len(dim_g) == 4:
                chi2_fitgf_lattice_normal_n4(g,dim_g,bath,dim_bath,ispin)
            elif len(dim_g) == 6:
                chi2_fitgf_lattice_normal_n6(g,dim_g,bath,dim_bath,ispin)
            else:
                raise ValueError("chi_fitgf_normal: takes dim(g) = 3 or 5")
        else:
            raise ValueError("chi_fitgf_normal: takes dim(bath) = 1 or 2")
    elif len(args) == 3: #superc
        g = args[0]
        f = args[1]
        bath = args[2]
        dim_g = np.asarray(np.shape(g),dtype=np.int64,order="F")
        dim_f = np.asarray(np.shape(g),dtype=np.int64,order="F")
        dim_bath = np.asarray(np.shape(bath),dtype=np.int64,order="F")
        if len(dim_bath) == 1: #single
            if len(dim_g) == 3:
                chi2_fitgf_single_superc_n3(g,dim_g,f,dim_f,bath,dim_bath,ispin,iorb,fmpi)
            elif len(dim_g) == 5:
                chi2_fitgf_single_superc_n5(g,dim_g,f,dim_f,bath,dim_bath,ispin,iorb,fmpi)
            else:
                raise ValueError("chi_fitgf_superc: takes dim(g,f) = 3 or 5")
        elif len(dim_bath) == 2: #lattice
            if len(dim_g) == 3:
                chi2_fitgf_lattice_superc_n3(g,dim_g,f,dim_f,bath,dim_bath,ispin)
            if len(dim_g) == 4:
                chi2_fitgf_lattice_superc_n4(g,dim_g,f,dim_f,bath,dim_bath,ispin)
            elif len(dim_g) == 6:
                chi2_fitgf_lattice_superc_n6(g,dim_g,f,dim_f,bath,dim_bath,ispin)
            else:
                raise ValueError("chi_fitgf_superc: takes dim(g,f) = 3 or 5")
        else:
            raise ValueError("chi_fitgf_superc: takes dim(bath) = 1 or 2")
    else:
        raise ValueError("chi_fitgf: takes g,bath or g,f,bath")
    return bath