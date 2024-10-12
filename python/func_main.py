from ctypes import *
import numpy as np
import os,sys
import types

# init_solver
def init_solver(self,bath=None,Nb=None,Nlat=None):
    if bath is None:
        if Nb is None and Nlat is None:
            Nb = self.library.get_bath_dimension()
            bath = np.zeros(Nb,dtype='float',order='F')
        elif Nb is None and Nlat is not None:
            Nb = self.library.get_bath_dimension()
            bath = np.zeros((Nlat,Nb),dtype='float',order='F')
        elif Nb is not None and Nlat is None:
            bath = np.zeros(Nb,dtype='float',order='F')
        elif Nb is not None and Nlat is not None:
            bath = np.zeros((Nlat,Nb),dtype='float',order='F')
    else:
        if Nb is not None or Nlat is not None:
            print("INIT_SOLVER WARNING: Bath vector provided, Nb and Nlat are discarded")

        
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
