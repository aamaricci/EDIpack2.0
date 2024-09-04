from ctypes import *
import numpy as np
import os,sys
import types

def init_solver(self,bath):
    init_solver_site = self.library.init_solver_site
    init_solver_site.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),
                                 np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS')]  
    init_solver_site.restype = None

    init_solver_ineq = self.library.init_solver_ineq
    init_solver_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS'),
                                 np.ctypeslib.ndpointer(dtype=np.int64,ndim=2, flags='F_CONTIGUOUS')]  
    init_solver_ineq.restype = None
    
    dim_bath=np.asarray(np.shape(bath),dtype=np.int64,order="F")
    
    if len(dim_bath)<2:
        init_solver_site(bath,dim_bath)
    else:
        init_solver_ineq(bath,dim_bath)

# Define the function signature for the Fortran function `solve_site`.
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