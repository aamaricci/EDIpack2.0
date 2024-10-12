from ctypes import *
import numpy as np
import os,sys
import types

# init_solver
def init_solver(self,*args):
    if not args:
        Nb = self.library.get_bath_dimension()
        bath = np.zeros(Nb,dtype='float',order='F')
    elif isinstance(args[0], np.ndarray) or isinstance(args[0], list):
        if isinstance(args[0][0], int):
            bath = np.zeros(args[0],dtype='float',order='F')
        elif isinstance(args[0][0], float):
            bath = args[0]
    elif isinstance(args[0], int):
        bath = np.zeros(args[0],dtype='float',order='F')
    else:
        print(type(args[0]))
        raise ValueError("init_solver: wrong input type.")
        
    init_solver_site = self.library.init_solver_site
    init_solver_site.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),
                                 np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS')]  
    init_solver_site.restype = None

    init_solver_ineq = self.library.init_solver_ineq
    init_solver_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS'),
                                 np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS')]  
    init_solver_ineq.restype = None
    
    dim_bath=np.asarray(np.shape(bath),dtype=np.int64,order="F")
    
    if len(dim_bath)<2:
        init_solver_site(bath,dim_bath)
        self.Nineq = 0
    else:
        init_solver_ineq(bath,dim_bath)
        self.Nineq = np.shape(bath)[0] #save number of inequivalent sites: this is useful when we want to know if we are in lattice case or not
    
    return bath

# `solve`.
def solve(self,bath,sflag=True,iflag=True,fmpi=True,mpi_lanc=False):
    solve_site = self.library.solve_site
    solve_site.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),
                           np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                           c_int,
                           c_int]  
    solve_site.restype = None

    # Define the function signature for the Fortran function `solve_ineq`.
    solve_ineq = self.library.solve_ineq
    solve_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS'),
                           np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                           c_int,
                           c_int]                      
    solve_ineq.restype = None  
  
    dim_bath=np.asarray(np.shape(bath),dtype=np.int64,order="F")
    
    if len(dim_bath)<2:
        solve_site(bath,dim_bath,sflag,fmpi)
    else:
        solve_ineq(bath,dim_bath,mpi_lanc,iflag)
        
#finalize solver
def finalize_solver(self):
    finalize_solver_wrapper = self.library.finalize_solver
    finalize_solver_wrapper.argtypes = [c_int]
    finalize_solver_wrapper.restype = None
    if self.Nineq is None:
        print("ED environment is not initialized yet")
        return ;
    else:
        finalize_solver_wrapper(self.Nineq)
        self.Nineq = None
        print("ED environment finalized")
        return ;
