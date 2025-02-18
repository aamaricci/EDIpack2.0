from ctypes import *
import numpy as np
import os, sys
import types


def chi2_fitgf(self, *args, ispin=0, iorb=None, fmpi=True):
    """
      This function fits the Weiss field or Hybridization function (delta) with \
      a discrete set of level. The fit parameters are the bath parameters \
      contained in the user-accessible array. Depending on the type of system \
      we are considering (normal, superconductive, non-SU(2)) \
      a different set of inputs has to be passed. The specifics of the numerical\
      fitting routines are controlled in the input file.
       
        
      :type args: [np.array(dtype=complex,np.array(dtype=complex), \
      np.array(dtype=float)]**or** [np.array(dtype=complex, \
      np.array(dtype=float)]
      :param args: The positional arguments are the function(s) \
      to fit and the bath array. 
       
       If the system is not superconductive ( :f:var:`ed_mode` = :code:`NORMAL` or\
        :f:var:`ed_mode` = :code:`NONSU2`) the argumens are
      
       * :code:`g`: the function to fit
       * :code:`bath`: the bath
      
       If the system is superconductive ( :f:var:`ed_mode` = :code:`SUPERC`) the \
       arguments are

       * :code:`g`: the normal function to fit
       * :code:`f`: the anomalous function to fit
       * :code:`bath`: the bath 
       
       The dimensions of the previous arrays can vary:
       
       The dimension of :code:`bath` can be
     
       * :code:`Nb`: if single-impurity, the output of :func:`get_bath_gimension`
       * :code:`[Nlat ,Nb]`: if real-space DMFT
       
       Accordingly, the dimension of g (and f) can be:
       
       * :code:`3`: in the single-impurity case,  an array of \
         the shape [ :data:`Nspin` :math:`\\cdot` :data:`Norb` ,  :data:`Nspin` \
         :math:`\\cdot` :data:`Norb` , :code:`Lfit` ]. 
       * :code:`3`: in the real-space DMFT case, an array of \
         the shape [ :code:`Nlat` :math:`\\cdot` :data:`Nspin` :math:`\\cdot` \
         :data:`Norb` ,  :code:`Nlat` :math:`\\cdot` :data:`Nspin` :math:`\\cdot` \
         :data:`Norb` , :code:`Lfit` ]
       * :code:`4`: in the real-space DMFT case, an array of \
         the shape [ :code:`Nlat` ,  :data:`Nspin` :math:`\\cdot` \
         :data:`Norb` ,  :data:`Nspin` :math:`\\cdot` :data:`Norb` , :code:`Lfit` ]
       * :code:`5`: in the single-impurity case, an array of \
         the shape [ :data:`Nspin` ,  :data:`Nspin` ,  :data:`Norb` ,  \
         :data:`Norb` , :code:`Lfit` ]
       * :code:`6`: in the real-space DMFT case, an array of \
         the shape [ :code:`Nlat` ,  :data:`Nspin` ,  :data:`Nspin` ,  \
         :data:`Norb` ,  :data:`Norb` , :code:`Lfit` ]
       
       where :code:`Lfit` is a given number of frequencies.

       
     
      :type ispin: int 
      :param ispin: spin species to be fitted. For the normal case, \
      if :data:`Nspin` = :code:`2`, the fitting function \
      needs to be called twice. Only the corresponding elements of :code:`bath` \
      will be updated each time. For the non-SU(2) case, this argument is \
      irrelevant, since all the elements of the Weiss/Delta function need to be \
      fitted. This is also the case if :f:var:`bath_type` = :code:`REPLICA, GENERAL`.
        
      :type iorb: int 
      :param iorb: the orbital to be fitted. If omitted, all orbitals will be fitted
       
      :type fmpi: bool 
      :param fmpi: flag to automatically do and broadcast the fit over MPI, if defined

      :raise ValueError: if the shapes of the positional arguments are incompatible
      :raise ValueError: if a number of positional arguments different \
       from 2 or 3 are passed   
         
      :return: An array of floats that contains the bath parameters \
       for the impurity problem. This is a required input of :func:`solve` \
       and :func:`chi2_fitgf`. Its elements are ordered differently \
       depending on the bath geometry. They are (de)compactified for \
       user interaction via :func:`bath_inspect`. Specific \
       symmetrization operations are implemented and listed in \
       the :ref:`bath` section.
      :rtype: np.array(dtype=float) 
    """

    # single normal
    chi2_fitgf_single_normal_n3 = self.library.chi2_fitgf_single_normal_n3
    chi2_fitgf_single_normal_n3.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=3, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=float, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        c_int,
        c_int,
        c_int,
    ]
    chi2_fitgf_single_normal_n3.restype = None

    chi2_fitgf_single_normal_n5 = self.library.chi2_fitgf_single_normal_n5
    chi2_fitgf_single_normal_n5.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=5, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=float, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        c_int,
        c_int,
        c_int,
    ]
    chi2_fitgf_single_normal_n5.restype = None

    # single superc
    chi2_fitgf_single_superc_n3 = self.library.chi2_fitgf_single_superc_n3
    chi2_fitgf_single_superc_n3.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=3, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=complex, ndim=3, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=float, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        c_int,
        c_int,
        c_int,
    ]
    chi2_fitgf_single_superc_n3.restype = None

    chi2_fitgf_single_superc_n5 = self.library.chi2_fitgf_single_superc_n5
    chi2_fitgf_single_superc_n5.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=5, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=complex, ndim=5, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=float, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        c_int,
        c_int,
        c_int,
    ]
    chi2_fitgf_single_superc_n5.restype = None

    # lattice normal
    chi2_fitgf_lattice_normal_n3 = self.library.chi2_fitgf_lattice_normal_n3
    chi2_fitgf_lattice_normal_n3.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=3, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=float, ndim=2, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        c_int,
    ]
    chi2_fitgf_lattice_normal_n3.restype = None

    chi2_fitgf_lattice_normal_n4 = self.library.chi2_fitgf_lattice_normal_n4
    chi2_fitgf_lattice_normal_n4.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=4, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=float, ndim=2, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        c_int,
    ]
    chi2_fitgf_lattice_normal_n4.restype = None

    chi2_fitgf_lattice_normal_n6 = self.library.chi2_fitgf_lattice_normal_n6
    chi2_fitgf_lattice_normal_n6.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=6, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=float, ndim=2, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        c_int,
    ]
    chi2_fitgf_lattice_normal_n6.restype = None

    # lattice superc
    chi2_fitgf_lattice_superc_n3 = self.library.chi2_fitgf_lattice_superc_n3
    chi2_fitgf_lattice_superc_n3.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=3, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=complex, ndim=3, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=float, ndim=2, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        c_int,
    ]
    chi2_fitgf_lattice_superc_n3.restype = None

    chi2_fitgf_lattice_superc_n4 = self.library.chi2_fitgf_lattice_superc_n4
    chi2_fitgf_lattice_superc_n4.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=4, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=complex, ndim=4, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=float, ndim=2, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        c_int,
    ]
    chi2_fitgf_lattice_superc_n4.restype = None

    chi2_fitgf_lattice_superc_n6 = self.library.chi2_fitgf_lattice_superc_n6
    chi2_fitgf_lattice_superc_n6.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=6, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=complex, ndim=6, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=float, ndim=2, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        c_int,
    ]
    chi2_fitgf_lattice_superc_n6.restype = None

    # main function
    ispin = ispin + 1
    if iorb is None:
        iorb = 0
    else:
        iorb = iorb + 1
    if len(args) == 2:  # normal
        g = np.asarray(args[0], order="F")
        bath = np.asarray(args[1], order="F")
        dim_g = np.asarray(np.shape(g), dtype=np.int64, order="F")
        dim_bath = np.asarray(np.shape(bath), dtype=np.int64, order="F")
        if len(dim_bath) == 1:  # single
            if len(dim_g) == 3:
                chi2_fitgf_single_normal_n3(g, dim_g, bath, dim_bath, ispin, iorb, fmpi)
            elif len(dim_g) == 5:
                chi2_fitgf_single_normal_n5(g, dim_g, bath, dim_bath, ispin, iorb, fmpi)
            else:
                raise ValueError("chi_fitgf_normal: takes dim(g) = 3 or 5")
        elif len(dim_bath) == 2:  # lattice
            if len(dim_g) == 3:
                chi2_fitgf_lattice_normal_n3(g, dim_g, bath, dim_bath, ispin)
            if len(dim_g) == 4:
                chi2_fitgf_lattice_normal_n4(g, dim_g, bath, dim_bath, ispin)
            elif len(dim_g) == 6:
                chi2_fitgf_lattice_normal_n6(g, dim_g, bath, dim_bath, ispin)
            else:
                raise ValueError("chi_fitgf_normal: takes dim(g) = 3 or 5")
        else:
            raise ValueError("chi_fitgf_normal: takes dim(bath) = 1 or 2")
    elif len(args) == 3:  # superc
        g = np.asarray(args[0], order="F")
        f = np.asarray(args[1], order="F")
        bath = np.asarray(args[2], order="F")
        dim_g = np.asarray(np.shape(g), dtype=np.int64, order="F")
        dim_f = np.asarray(np.shape(g), dtype=np.int64, order="F")
        dim_bath = np.asarray(np.shape(bath), dtype=np.int64, order="F")
        if len(dim_bath) == 1:  # single
            if len(dim_g) == 3:
                chi2_fitgf_single_superc_n3(
                    g, dim_g, f, dim_f, bath, dim_bath, ispin, iorb, fmpi
                )
            elif len(dim_g) == 5:
                chi2_fitgf_single_superc_n5(
                    g, dim_g, f, dim_f, bath, dim_bath, ispin, iorb, fmpi
                )
            else:
                raise ValueError("chi_fitgf_superc: takes dim(g,f) = 3 or 5")
        elif len(dim_bath) == 2:  # lattice
            if len(dim_g) == 3:
                chi2_fitgf_lattice_superc_n3(g, dim_g, f, dim_f, bath, dim_bath, ispin)
            if len(dim_g) == 4:
                chi2_fitgf_lattice_superc_n4(g, dim_g, f, dim_f, bath, dim_bath, ispin)
            elif len(dim_g) == 6:
                chi2_fitgf_lattice_superc_n6(g, dim_g, f, dim_f, bath, dim_bath, ispin)
            else:
                raise ValueError("chi_fitgf_superc: takes dim(g,f) = 3 or 5")
        else:
            raise ValueError("chi_fitgf_superc: takes dim(bath) = 1 or 2")
    else:
        raise ValueError("chi_fitgf: takes g,bath or g,f,bath")
    return bath
