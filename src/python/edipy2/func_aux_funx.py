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
       
        * [ :data:`Nspin` :math:`\\cdot` :data:`Norb` , :data:`Nspin` :math:`\\cdot` :data:`Norb` ]: single-impurity case, 2-dimensional array
        * [ :data:`Nspin` , :data:`Nspin` , :data:`Norb` , :data:`Norb` ]: single-impurity case, 4-dimensional array
        * [ :code:`Nlat` :math:`\\cdot` :data:`Nspin` :math:`\\cdot` :data:`Norb` , :code:`Nlat` :math:`\\cdot` :data:`Nspin` :math:`\\cdot` :data:`Norb` ]: real-space DMFT case, 2-dimensional array.
        * [ :code:`Nlat` , :data:`Nspin` :math:`\\cdot` :data:`Norb` , :data:`Nspin` :math:`\\cdot` :data:`Norb` ]: single-impurity case, 3-dimensional array.
        * [ :code:`Nlat` , :data:`Nspin` , :data:`Nspin` , :data:`Norb` , :data:`Norb` ]: single-impurity case, 5-dimensional array.
       
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
     
     This function checks the value of the read density :code:`ntmp` against the desired value :data:`nread` (if different from zero) and adjusts :code:`var` accordingly (in a monotonous way).
   
     :type var: float
     :param var: the variable to be adjusted (usually :data:`xmu` )

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
def check_convergence(self,func,threshold,N1=None,N2=None):
    """
    
    This function checks the variation of a given quantity (Weiss field, Delta, ...) against the one for the previous step. It is used to determined whether the DMFT loop has converged. If a maximum number of loops is exceeded, returns True with a warning.

    :type func: np.array(dtype=complex) 
    :param func: the quantity to be checked. It is one-dimensional, with its length being a number of frequencies
   
    :type threshold: float 
    :param threshold: the error threshold
   
    :type N1: int
    :param N1: minimum number of loops (default = :data:`edipy2.global_env.Nsuccess`)

    :type N2: int
    :param N2: maximum number of loops (default = :data:`edipy2.global_env.Nloop`)
   
    :return: 
     - the error
     - a boolean signifying convergence
    :rtype: float, bool
    
    """
    try:
        import mpi4py
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
    except:
        rank = 0
        
    func = np.asarray(func)
    err = 1.0
    conv_bool = False
    outfile = "error.err"
    
    #if N1 and/or N2 are None, set them to the input variables
    if N1 is None:
        N1 = c_int.in_dll(self.library, "Nsuccess").value
    if N2 is None:
        N2 = c_int.in_dll(self.library, "Nloop").value
    
    #if first loop, allocate old function as method
    if not hasattr(self,"oldfunc"):
        self.oldfunc = np.zeros_like(func,dtype=complex)
        self.whichiter = 1
        self.gooditer = 0
    
    #only the master does the calculation      
    if rank == 0:
        errvec = np.real(np.sum(abs(func - self.oldfunc),axis = -1) / np.sum(abs(func),axis = -1))
        #first iteration
        if self.whichiter == 1:
            errvec = np.ones_like(errvec)        
        #remove nan compoments, if some component is divided by zero
        if np.prod(np.shape(errvec)) > 1:
            errvec = errvec[~np.isnan(errvec)]       
        errmax = np.max(errvec)
        errmin = np.min(errvec)
        err = np.average(errvec)
        self.oldfunc = np.copy(func)
        conv_bool = ((err < threshold) and (self.gooditer >= N1) and (self.whichiter <= N2)) or (self.whichiter > N2)
        if err < threshold:
            self.gooditer += 1
        self.whichiter += 1
        
        #write out
        with open(outfile, "a") as file:
            file.write(f"{self.whichiter} {err:.6e}\n")
        if np.prod(np.shape(errvec)) > 1:
            with open(outfile+".max", "a") as file:
                file.write(f"{self.whichiter} {errmax:.6e}\n")
            with open(outfile+".min", "a") as file:
                file.write(f"{self.whichiter} {errmin:.6e}\n")
            with open(outfile+".distribution", "a") as file:
                file.write(f"{self.whichiter}" + " ".join([f"{x:.6e}" for x in errvec.flatten()]) + "\n")
        
        #print convergence message:
        if conv_bool:
            colorprefix= self.BOLD + self.GREEN
        elif (err < threshold) and (self.gooditer <= N1):
            colorprefix= self.BOLD + self.YELLOW
        else:
            colorprefix= self.BOLD + self.RED
        
        print("\n")
        if self.whichiter <= N2:
            if np.prod(np.shape(errvec)) > 1:
                print(colorprefix + "max error=" + self.COLOREND + f"{errmax:.6e}")
            print(colorprefix + "    error=" + self.COLOREND + f"{err:.6e}")
            if np.prod(np.shape(errvec)) > 1:
                print(colorprefix + "min error=" + self.COLOREND + f"{errmin:.6e}\n")
        else:
            print("Not converged after "+str(N2)+" iterations.")
            with open("ERROR.ERR", "w") as file:
                file.write("Not converged after "+str(N2)+" iterations.")
        
        #delete methods after convergence or Nmax exceeded
        if conv_bool:
            del self.oldfunc
            del self.gooditer
            del self.whichiter
    
    #pass to other cores:   
    try:
        conv_bool = comm.bcast(conv_bool, root=0)
        err = comm.bcast(err, root=0)
    except:
        pass
    return err,conv_bool
