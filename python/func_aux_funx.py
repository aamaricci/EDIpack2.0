from ctypes import *
import numpy as np
import os,sys
import types

#set_hloc

def set_hloc(self,hloc,Nlat=None):
    ed_set_Hloc_single_N2 = self.library.ed_set_Hloc_single_N2
    ed_set_Hloc_single_N2.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=2, flags='F_CONTIGUOUS'),
                                np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS')] 
    ed_set_Hloc_single_N2.restype = None

    ed_set_Hloc_single_N4 = self.library.ed_set_Hloc_single_N4
    ed_set_Hloc_single_N4.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=4, flags='F_CONTIGUOUS'),
                                np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS')] 
    ed_set_Hloc_single_N4.restype = None

    ed_set_Hloc_lattice_N2 = self.library.ed_set_Hloc_lattice_N2
    ed_set_Hloc_lattice_N2.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=2, flags='F_CONTIGUOUS'),
                                np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                c_int] 
    ed_set_Hloc_lattice_N2.restype = None

    ed_set_Hloc_lattice_N3 = self.library.ed_set_Hloc_lattice_N3
    ed_set_Hloc_lattice_N3.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=3, flags='F_CONTIGUOUS'),
                                np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                c_int] 
    ed_set_Hloc_lattice_N3.restype = None

    ed_set_Hloc_lattice_N5 = self.library.ed_set_Hloc_lattice_N5
    ed_set_Hloc_lattice_N5.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=5, flags='F_CONTIGUOUS'),
                                np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                c_int] 
    ed_set_Hloc_lattice_N5.restype = None
    
    try:
        hloc = np.asarray(hloc,order="F")
        dim_hloc = np.asarray(np.shape(hloc),dtype=np.int64,order="F")
        self.dim_hloc = len(dim_hloc)
    except:
        raise ValueError("In Edipack2.0, set_Hloc needs an Hloc defined")
    
    if(Nlat is not None):
        if len(dim_hloc) == 2:
            ed_set_Hloc_lattice_N2(hloc,dim_hloc,Nlat)
        elif len(dim_hloc) == 3:
            ed_set_Hloc_lattice_N3(hloc,dim_hloc,Nlat)
        elif len(dim_hloc) == 5:
            ed_set_Hloc_lattice_N5(hloc,dim_hloc,Nlat)
        else:
            raise ValueError ("ed_set_Hloc_lattice: dimension must be 2,3 or 5")
    else:
        if len(dim_hloc) == 2:
            ed_set_Hloc_single_N2(hloc,dim_hloc)
        elif len(dim_hloc) == 4:
            ed_set_Hloc_single_N4(hloc,dim_hloc)
        else:
            raise ValueError ("ed_set_Hloc_site: dimension must be 2 or 4")
    return ;


#search_variable
def search_variable(self,var,ntmp,converged):
    search_variable_wrap = self.library.search_variable
    search_variable_wrap.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),
                                     np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),
                                     np.ctypeslib.ndpointer(dtype=int,ndim=1, flags='F_CONTIGUOUS')]
    search_variable_wrap.restype = None
    var = np.asarray([var])
    ntmp = np.asarray([ntmp])
    converged = np.asarray([converged])
    conv_int=int(converged)
    search_variable_wrap(var,ntmp,converged)
    if conv_int[0]==0:
        converged=False
    else:
        converged=True
    return var[0],conv_bool

#check_convergence
def check_convergence(self,func,threshold,N1,N2):
    check_convergence_wrap = self.library.check_convergence
    check_convergence_wrap.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=1, flags='F_CONTIGUOUS'),
                                       c_int,
                                       c_double,
                                       c_int,
                                       c_int,
                                       np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),
                                       np.ctypeslib.ndpointer(dtype=int,ndim=1, flags='F_CONTIGUOUS')]
    check_convergence_wrap.restype = None
    err=np.asarray([1.0])
    converged=np.asarray([0])
    func=np.asarray(func,order="F")
    dim_func=np.shape(func)
    check_convergence_wrap(func,dim_func[0],threshold,N1,N2,err,converged)
    if converged[0]==0:
        conv_bool=False
    else:
        conv_bool=True
    return err[0],conv_bool
