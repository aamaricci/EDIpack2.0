from ctypes import *
import numpy as np
import os, sys
import types


# init_solver
def init_solver(self, bath=None, Nb=None, Nlat=None):
    """

       This function initializes the ED environment for the impurity problem \
       solution, and sets the bath reading it from the ``hamiltonian.restart``\
       file or initializing it in a symmetric way.
       The function can take different argument combinations. 
       
       If no input is provided, a single-impurity bath is allocated, with \
       dimension given by :func:`get_bath_dimension`.

        
       :type bath: np.array(dtype=float) **or** [float]
       :param bath: If a bath array is provided, it has to be a numpy array \
       or a tuple of floats. It has to have one or two dimensions. If it has \
       one dimension, that must be the same as specified by :func:`get_bath_dimension`. \
       If it has two dimensions, the first has to be the number of inequivalent \
       sites for real-space DMFT, the second must be in agreement with \
       :func:`get_bath_dimension`. If ``Nlat`` or ``Nb`` are provided, \
       this overrides them. If the provided vector is not in agreement \
       with the global system parameters, EDIpack2 will exit with an error.
        The array is ordered in F convention inside the function.
            
       :type Nb: int 
       :param Nb: This sets the bath vector length for each single impurity \
       problem. It has to be in agreement with :func:`get_bath_dimension`. \
       When this parameter alone is provided, a numpy array of this length \
       will be initialized.
        
       :type Nlat: int 
       :param Nlat: This sets the number of inequivalent sites for \
       real-space DMFT. If this parameter alone is provided, \
       :func:`get_bath_dimension` is invoked to determine the bath vector \
       length Nb for each impurity. A ``[Nlat,Nb]`` vector is then allocated.
       
         
       :return: An array of floats that contains the bath parameters for the \
       impurity problem. This is a required input of :func:`solve` and \
       :func:`chi2_fitgf`. Its elements are ordered differently depending \
       on the bath geometry. They are (de)compactified for user interaction \
       via :func:`bath_inspect`. Specific symmetrization operations are \
       implemented and listed in the :ref:`bath` section.
       :rtype: np.array(dtype=float) 
    """
    if bath is None:
        if Nb is None and Nlat is None:
            Nb = self.library.get_bath_dimension()
            bath = np.zeros(Nb, dtype="float", order="F")
        elif Nb is None and Nlat is not None:
            Nb = self.library.get_bath_dimension()
            bath = np.zeros((Nlat, Nb), dtype="float", order="F")
        elif Nb is not None and Nlat is None:
            bath = np.zeros(Nb, dtype="float", order="F")
        elif Nb is not None and Nlat is not None:
            bath = np.zeros((Nlat, Nb), dtype="float", order="F")
    else:
        if Nb is not None or Nlat is not None:
            print(
                "INIT_SOLVER WARNING: Bath vector provided, Nb and Nlat are discarded"
            )

    init_solver_site = self.library.init_solver_site
    init_solver_site.argtypes = [
        np.ctypeslib.ndpointer(dtype=float, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
    ]
    init_solver_site.restype = None

    init_solver_ineq = self.library.init_solver_ineq
    init_solver_ineq.argtypes = [
        np.ctypeslib.ndpointer(dtype=float, ndim=2, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
    ]
    init_solver_ineq.restype = None

    dim_bath = np.asarray(np.shape(bath), dtype=np.int64, order="F")

    if len(dim_bath) < 2:
        init_solver_site(bath, dim_bath)
        self.Nineq = 0
    else:
        init_solver_ineq(bath, dim_bath)
        self.Nineq = np.shape(bath)[
            0
        ]  # save number of inequivalent sites: this is useful when we want to know if we are in lattice case or not

    return bath


# `solve`.
def solve(self, bath, flag_gf=True, flag_mpi=True, mpi_lanc=False):
    """
       This function solves the impurity problem and calculates the \
       observables, Green's function and self-energy.

       :type bath: np.array(dtype=float) 
       :param bath: The bath array returned by  :func:`init_solver`. \
       If the bath dimensions are inconsistent with the global properties \
       of the problem, EDIpack2 will exit with an error.
       
       :type flag_gf: bool
       :param flag_gf: for single-impurity DMFT, if :code:`False`, it disables \
       the calculation of the Green's function and susceptibilities
       
       :type flag_mpi: bool
       :param flag_mpi: if :code:`False`, for single-impurity DMFT, it disables \
       MPI for the ED routine, if the communicator is used elsewhere
       
       :type mpi_lanc: bool
       :param mpi_lanc: if :code:`True`, for real-space DMFT sets the MPI \
       parallelization for the ED routine. By default it is :code:`False`, \
       and each inequivalent site is solved serially by a different core.
            
       :return: Nothing
       :rtype: None
    """
    solve_site = self.library.solve_site
    solve_site.argtypes = [
        np.ctypeslib.ndpointer(dtype=float, ndim=1, flags="F_CONTIGUOUS"),  # bath
        np.ctypeslib.ndpointer(
            dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"
        ),  # dim_bath
        c_int,  # flag_gf
        c_int,
    ]  # flag_mpi
    solve_site.restype = None

    # Define the function signature for the Fortran function `solve_ineq`.
    solve_ineq = self.library.solve_ineq
    solve_ineq.argtypes = [
        np.ctypeslib.ndpointer(dtype=float, ndim=2, flags="F_CONTIGUOUS"),  # bath
        np.ctypeslib.ndpointer(
            dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"
        ),  # dim_bath
        c_int,  # flag_gf
        c_int,
    ]  # flag_mpi
    solve_ineq.restype = None

    dim_bath = np.asarray(np.shape(bath), dtype=np.int64, order="F")

    if len(dim_bath) < 2:
        solve_site(bath, dim_bath, flag_gf, flag_mpi)
    else:
        solve_ineq(bath, dim_bath, flag_gf, mpi_lanc)


# finalize solver
def finalize_solver(self):
    """
       This function cleans up the ED environment, deallocates the relevant \
       arrays and makes a second call to :command:`init_solver` possible
               
       :return: Nothing
       :rtype: None
    """

    finalize_solver_wrapper = self.library.finalize_solver
    finalize_solver_wrapper.argtypes = [c_int]
    finalize_solver_wrapper.restype = None
    if self.Nineq is None:
        print("ED environment is not initialized yet")
        return
    else:
        finalize_solver_wrapper(self.Nineq)
        self.Nineq = None
        self.dim_hloc = 0
        self.Nsym = None

        if hasattr(self, "oldfunc"):
            del self.oldfunc
        if hasattr(self, "gooditer"):
            del self.gooditer
        if hasattr(self, "whichiter"):
            del self.whichiter

        print("ED environment finalized")

        return
