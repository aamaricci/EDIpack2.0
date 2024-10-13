from ctypes import *
import numpy as np
import os,sys
import types

#set_hloc

def set_hloc(self,hloc,Nlat=None):
    """
       This function sets the local Hamiltonian of the impurity problem. 
        
       :type hloc: np.array(dtype=complex)
       :param hloc: Local Hamiltonian matrix. This can have the following shapes:
       
        * :code:`[ed.Nspin*ed.Norb,ed.Nspin*ed.Norb]`: single-impurity case, 2-dimensional array
        * :code:`[ed.Nspin,ed.Nspin,ed.Norb,ed.Norb]`: single-impurity case, 4-dimensional array
        * :code:`[Nlat*ed.Nspin*ed.Norb,Nlat*ed.Nspin*ed.Norb]`: real-space DMFT case, 2-dimensional array.
        * :code:`[Nlat,Nspin*ed.Norb,ed.Nspin*ed.Norb]`: single-impurity case, 3-dimensional array.
        * :code:`[Nlat,ed.Nspin,ed.Nspin,ed.Norb,ed.Norb]`: single-impurity case, 5-dimensional array.
       
        The array is ordered in F convention inside the function.
        
        **Note**: the way the EDIpack2 library passes from 1 comulative to 2 or 3 running indices is, from slower to faster: ``lat``, ``spin``, ``orb``
        
       :type Nlat: int
       :param Nlat: Number of inequivalent sites for real-space DMFT. The function will raise a ValueError if the dimensions of ``hloc`` are inconsistent with the presence or absence of Nlat. 
        The EDIpack2 library will check the correctness of the dimensions of ``hloc`` and terminate execution if inconsistent.
       
       :raise ValueError: If hloc is not provided or has the wrong shape
       
       :return: Nothing
       :rtype: None
    """
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
    """
     
     This function checks the value of the read density :code:`ntmp` against the desired value :code:`ed.nread` (if different from zero) and adjusts :code:`var` accordingly (in a monotonous way).
   
     :type var: float
     :param var: the variable to be adjusted (usually :code:`ed.xmu`)

     :type ntmp: float
     :param ntmp: the density value at the given iteration
   
     :type converged: bool
     :param converged: whether the DMFT loop has achieved a sufficiently small error independently on the density
   
     :return: 
      - the new value of :code:`var`
      - a boolean signifying convergence
     :rtype: float, bool
     
    """
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
    """
    
    This function checks the variation of a given quantity (Weiss field, Delta, ...) against the one for the previous step. It is used to determined whether the DMFT loop has converged. If a maximum number of loops is exceeded, returns True with a warning.

    :type func: np.array(dtype=complex) 
    :param func: the quantity to be checked. It is one-dimensional, with its length being a number of frequencies
   
    :type threshold: float 
    :param threshold: the error threshold
   
    :type N1: int
    :param N1: minimum number of loops

    :type N2: int
    :param N2: maximum number of loops
   
    :return: 
     - the error
     - a boolean signifying convergence
    :rtype: float, bool
    
    """
    
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
