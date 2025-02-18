from ctypes import *
import numpy as np
import os, sys
import types


# get_bath_dimension
def get_bath_dimension(self):
    """
    This function returns the correct dimension for the bath to be allocated (for each impurity) given the parameters of the system.

    :return: a number which is the dimension of the bath array for each impurity.
    :rtype: int
    """
    get_bath_dimension_wrap = self.library.get_bath_dimension
    get_bath_dimension_wrap.argtypes = None
    get_bath_dimension_wrap.restype = c_int
    return get_bath_dimension_wrap()


# init_hreplica
def set_hreplica(self, hvec, lambdavec):
    """

       This function is specific to :f:var:`bath_type` = :code:`=replica`. It sets the basis of matrices\
       and scalar parameters that, upon linear combination, make up the bath replica.
        
       :type hvec: np.array(dtype=complex)
       :param hvec: array of bath matrices. They decompose the nonzero part of the replica in a set.\
       Each element of the set correspond to a variational parameter.\
       That way the bath replica matrix is updated while preserving symmetries\
       of the user's choosing. The array can have the following shapes:

        * [ :code:`(Nnambu)` :math:`\\cdot` :data:`Nspin` :math:`\\cdot` :data:`Norb` , :code:`(Nnambu)` :math:`\\cdot` :data:`Nspin` :math:`\\cdot` :data:`Norb` , :code:`Nsym` ]:\
        3-dimensional, where Nnambu refers to the superconducting case and Nsym \
        is the number of matrices that make up the linear combination 
        * [:code:`(Nnambu)` :math:`\\cdot` :data:`Nspin` , :code:`(Nnambu)`  :data:`Nspin` ,  :data:`Norb` ,  :data:`Norb` , :code:`Nsym` ]:\
        5-dimensional, where Nnambu refers to the superconducting case and Nsym is \
        the number of matrices that make up the linear combination 
        
       :type lambdavec: np.array(dtype=float) 
       :param lambdavec: the array of coefficients of the linear combination.\
       This, along with the hybridizations V, are the fitting parameters of \
       the bath. The array has the following shape
       
        * [ :data:`Nbath` , :code:`Nsym` ]: for single-impurity DMFT, 2-dimensional,\
        where Nsym is the number of matrices that make up the linear combination 
        * [ :code:`Nlat`,  :data:`Nbath` , :code:`Nsym` ]: for real-space DMFT, 3-dimensional,\
        where Nlat is the number of inequivalent impurity sites and Nsym is\
        the number of matrices that make up the linear combination 

       :raise ValueError: if the shapes of the arrays are inconsistent
         
       :return: Nothing
       :rtype: None
    """
    init_hreplica_symmetries_d5 = self.library.init_Hreplica_symmetries_d5
    init_hreplica_symmetries_d5.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=5, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=float, ndim=2, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
    ]
    init_hreplica_symmetries_d5.restype = None

    init_hreplica_symmetries_d3 = self.library.init_Hreplica_symmetries_d3
    init_hreplica_symmetries_d3.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=3, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=float, ndim=2, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
    ]
    init_hreplica_symmetries_d3.restype = None

    init_hreplica_symmetries_lattice_d5 = (
        self.library.init_Hreplica_symmetries_lattice_d5
    )
    init_hreplica_symmetries_lattice_d5.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=5, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=float, ndim=3, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
    ]
    init_hreplica_symmetries_lattice_d5.restype = None

    init_hreplica_symmetries_lattice_d3 = (
        self.library.init_Hreplica_symmetries_lattice_d3
    )
    init_hreplica_symmetries_lattice_d3.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=3, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=float, ndim=3, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
    ]
    init_hreplica_symmetries_lattice_d3.restype = None

    aux_norb = c_int.in_dll(self.library, "Norb").value
    aux_nspin = c_int.in_dll(self.library, "Nspin").value
    dim_hvec = np.asarray(np.shape(hvec), dtype=np.int64, order="F")
    dim_lambdavec = np.asarray(np.shape(lambdavec), dtype=np.int64, order="F")

    self.Nsym = dim_lambdavec[1]

    if len(dim_hvec) == 3:
        if len(dim_lambdavec) == 2:
            init_hreplica_symmetries_d3(hvec, dim_hvec, lambdavec, dim_lambdavec)
        elif len(Ddim_lambdavec) == 3:
            init_hreplica_symmetries_lattice_d3(
                hvec, dim_hvec, lambdavec, dim_lambdavec
            )
        else:
            raise ValueError("Shape(lambdavec) != 2 or 3  in set_Hreplica")
    elif len(dim_hvec) == 5:
        if len(dim_lambdavec) == 2:
            init_hreplica_symmetries_d5(hvec, dim_hvec, lambdavec, dim_lambdavec)
        elif len(dim_lambdavec) == 3:
            init_hreplica_symmetries_lattice_d5(
                hvec, dim_hvec, lambdavec, dim_lambdavec
            )
        else:
            raise ValueError("Shape(lambdavec) != 2 or 3  in set_Hreplica")
    else:
        raise ValueError("Shape(Hvec) != 3 or 5  in set_Hreplica")
    return


# init_hgeneral
def set_hgeneral(self, hvec, lambdavec):
    """
       This function is specific to :code:`BATH_TYPE=GENERAL`. It sets the basis of matrices\
       and scalar parameters that, upon linear combination, make up the bath replica. \
       The input is the same as that of :func:`set_hreplica`.
        
       :type hvec: np.array(dtype=complex)
       :param hvec: array of bath matrices. They decompose the nonzero part of the replica in a set.\
       Each element of the set correspond to a variational parameter.\
       That way the bath replica matrix is updated while preserving symmetries\
       of the user's choosing. The array can have the following shapes:

        * [:code:`(Nnambu)` :math:`\\cdot`:data:`Nspin` :math:`\\cdot` :data:`Norb` , :code:`(Nnambu)` :math:`\\cdot` :data:`Nspin` :math:`\\cdot` :data:`Norb` , :code:`Nsym` ]:\
        3-dimensional, where Nnambu refers to the superconducting case and Nsym \
        is the number of matrices that make up the linear combination 
        * [:code:`(Nnambu)` :math:`\\cdot` :data:`Nspin` , :code:`(Nnambu)` :math:`\\cdot` :data:`Nspin` ,  :data:`Norb` ,  :data:`Norb` , :code:`Nsym` ]:\
        5-dimensional, where Nnambu refers to the superconducting case and Nsym is \
        the number of matrices that make up the linear combination 
        
       :type lambdavec: np.array(dtype=float) 
       :param lambdavec: the array of coefficients of the linear combination.\
       This, along with the hybridizations V, are the fitting parameters of \
       the bath. The array has the following shape
       
        * [ :data:`Nbath` , :code:`Nsym` ]: for single-impurity DMFT, 2-dimensional,\
        where Nsym is the number of matrices that make up the linear combination 
        * [ :code:`Nlat`,  :data:`Nbath` , :code:`Nsym` ]: for real-space DMFT, 3-dimensional,\
        where Nlat is the number of inequivalent impurity sites and Nsym is\
        the number of matrices that make up the linear combination 

       :raise ValueError: if the shapes of the arrays are inconsistent
         
       :return: Nothing
       :rtype: None
    """
    init_hgeneral_symmetries_d5 = self.library.init_Hgeneral_symmetries_d5
    init_hgeneral_symmetries_d5.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=5, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=float, ndim=2, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
    ]
    init_hgeneral_symmetries_d5.restype = None

    init_hgeneral_symmetries_d3 = self.library.init_Hgeneral_symmetries_d3
    init_hgeneral_symmetries_d3.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=3, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=float, ndim=2, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
    ]
    init_hgeneral_symmetries_d3.restype = None

    init_hgeneral_symmetries_lattice_d5 = (
        self.library.init_Hgeneral_symmetries_lattice_d5
    )
    init_hgeneral_symmetries_lattice_d5.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=5, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=float, ndim=3, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
    ]
    init_hgeneral_symmetries_lattice_d5.restype = None

    init_hgeneral_symmetries_lattice_d3 = (
        self.library.init_Hgeneral_symmetries_lattice_d3
    )
    init_hgeneral_symmetries_lattice_d3.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=3, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=float, ndim=3, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
    ]
    init_hgeneral_symmetries_lattice_d3.restype = None

    aux_norb = c_int.in_dll(self.library, "Norb").value
    aux_nspin = c_int.in_dll(self.library, "Nspin").value
    dim_hvec = np.asarray(np.shape(hvec), dtype=np.int64, order="F")
    dim_lambdavec = np.asarray(np.shape(lambdavec), dtype=np.int64, order="F")

    self.Nsym = dim_lambdavec[1]

    if len(dim_hvec) == 3:
        if len(dim_lambdavec) == 2:
            init_hgeneral_symmetries_d3(hvec, dim_hvec, lambdavec, dim_lambdavec)
        elif len(Ddim_lambdavec) == 3:
            init_hgeneral_symmetries_lattice_d3(
                hvec, dim_hvec, lambdavec, dim_lambdavec
            )
        else:
            raise ValueError("Shape(lambdavec) != 2 or 3  in set_Hgeneral")
    elif len(dim_hvec) == 5:
        if len(dim_lambdavec) == 2:
            init_hgeneral_symmetries_d5(hvec, dim_hvec, lambdavec, dim_lambdavec)
        elif len(dim_lambdavec) == 3:
            init_hgeneral_symmetries_lattice_d5(
                hvec, dim_hvec, lambdavec, dim_lambdavec
            )
        else:
            raise ValueError("Shape(lambdavec) != 2 or 3  in set_Hgeneral")
    else:
        raise ValueError("Shape(Hvec) != 3 or 5  in set_Hgeneral")
    return


# break_symmetry_bath
def break_symmetry_bath(self, bath, field, sign, save=True):
    """
    
    This function breaks the spin symmetry of the bath, useful \
    for magnetic calculations to incite symmetry breaking.\
    Not compatible with :code:`REPLICA` or :code:`GENERAL` bath types.

    :type bath: np.array(dtype=float)
    :param bath: The user-accessible bath array
   
    :type field: float
    :param field: the magnitude of the symmetry-breaking shift
   
    :type sign: float
    :param sign: the sign of the symmetry-breaking shift
   
    :type save: bool
    :param save: whether to save the symmetry-broken bath for reading
   
    :return: the modified bath array
    :rtype: np.array(dtype=float) 
       
    """

    break_symmetry_bath_site = self.library.break_symmetry_bath_site
    break_symmetry_bath_site.argtypes = [
        np.ctypeslib.ndpointer(dtype=float, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        c_double,
        c_double,
        c_int,
    ]
    break_symmetry_bath_site.restype = None

    break_symmetry_bath_ineq = self.library.break_symmetry_bath_ineq
    break_symmetry_bath_ineq.argtypes = [
        np.ctypeslib.ndpointer(dtype=float, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        c_double,
        c_double,
        c_int,
    ]
    break_symmetry_bath_ineq.restype = None

    if save:
        save_int = 1
    else:
        save_int = 0
    bath_shape = np.asarray(np.shape(bath), dtype=np.int64, order="F")
    if (len(bath_shape)) == 1:
        break_symmetry_bath_site(bath, bath_shape, field, float(sign), save_int)
    else:
        break_symmetry_bath_ineq(bath, bath_shape, field, float(sign), save_int)
    return bath


# spin_symmetrize_bath


def spin_symmetrize_bath(self, bath, save=True):
    """
       This function enforces equality of the opposite-spin components\
       of the bath array. Not compatible with :code:`REPLICA` or :code:`GENERAL` bath types.

       :type bath: np.array(dtype=float)
       :param bath: The user-accessible bath array
          
       :type save: bool
       :param save: whether to save the symmetry-broken bath for reading
       
       :return: the modified bath array
       :rtype: np.array(dtype=float)
    """
    spin_symmetrize_bath_site = self.library.spin_symmetrize_bath_site
    spin_symmetrize_bath_site.argtypes = [
        np.ctypeslib.ndpointer(dtype=float, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        c_int,
    ]
    spin_symmetrize_bath_site.restypes = None

    spin_symmetrize_bath_ineq = self.library.spin_symmetrize_bath_ineq
    spin_symmetrize_bath_ineq.argtypes = [
        np.ctypeslib.ndpointer(dtype=float, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        c_int,
    ]
    spin_symmetrize_bath_ineq.restypes = None
    if save:
        save_int = 1
    else:
        save_int = 0
    bath_shape = np.asarray(np.shape(bath), dtype=np.int64, order="F")
    if (len(bath_shape)) == 1:
        spin_symmetrize_bath_site(bath, bath_shape, save_int)
    else:
        spin_symmetrize_bath_ineq(bath, bath_shape, save_int)
    return bath


# orb_symmetrize_bath
def orb_symmetrize_bath(self, bath, orb1, orb2, save=True):
    """
       This function enforces equality of the different-orbital components \
       of the bath array. Not compatible with :code:`REPLICA` or :code:`GENERAL` bath types.

       :type bath: np.array(dtype=float)
       :param bath: The user-accessible bath array
       
       :type orb1: int
       :param orb1: first orbital index
       
       :type orb2: int
       :param orb2: second orbital index
          
       :type save: bool
       :param save: whether to save the symmetry-broken bath for reading
       
       :return: the modified bath array
       :rtype: np.array(dtype=float)
    """

    orb_symmetrize_bath_site = self.library.orb_symmetrize_bath_site
    orb_symmetrize_bath_site.argtypes = [
        np.ctypeslib.ndpointer(dtype=float, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        c_int,
    ]
    orb_symmetrize_bath_site.restypes = None

    orb_symmetrize_bath_ineq = self.library.orb_symmetrize_bath_ineq
    orb_symmetrize_bath_ineq.argtypes = [
        np.ctypeslib.ndpointer(dtype=float, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        c_int,
    ]
    orb_symmetrize_bath_ineq.restypes = None

    if save:
        save_int = 1
    else:
        save_int = 0
    bath_shape = np.asarray(np.shape(bath), dtype=np.int64, order="F")
    if (len(bath_shape)) == 1:
        orb_symmetrize_bath_site(bath, bath_shape, orb1 + 1, orb2 + 1, save_int)
    else:
        orb_symmetrize_bath_ineq(bath, bath_shape, orb1 + 1, orb2 + 1, save_int)
    return bath


# orb_equality_bath


def orb_equality_bath(self, bath, indx, save=True):
    """
       This function sets every orbital component to be equal to the \
       one of orbital :code:`indx`. Not compatible with :code:`REPLICA` or :code:`GENERAL` bath types.

       :type bath: np.array(dtype=float)
       :param bath: The user-accessible bath array
       
       :type iorb: int 
       :param iorb: the orbital index to which every other will be set as equal
          
       :type save: bool
       :param save: whether to save the symmetry-broken bath for reading
       
       :raise ValueError: if the orbital index is out of bounds
       
       :return: the modified bath array
       :rtype: np.array(dtype=float) 
    """
    orb_equality_bath_site = self.library.orb_equality_bath_site
    orb_equality_bath_site.argtypes = [
        np.ctypeslib.ndpointer(dtype=float, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        c_int,
        c_int,
    ]
    orb_equality_bath_site.restypes = None

    orb_equality_bath_ineq = self.library.orb_equality_bath_ineq
    orb_equality_bath_ineq.argtypes = [
        np.ctypeslib.ndpointer(dtype=float, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        c_int,
        c_int,
    ]
    orb_equality_bath_ineq.restypes = None
    aux_norb = c_int.in_dll(self.library, "Norb").value
    if save:
        save_int = 1
    else:
        save_int = 0
    bath_shape = np.asarray(np.shape(bath), dtype=np.int64, order="F")
    if (indx < 0) or (indx >= aux_norb):
        raise ValueError("orb_equality_bath: orbital index should be in [0,Norb]")
    else:
        indx = indx + 1  # python to fortran convention
        if (len(bath_shape)) == 1:
            orb_equality_bath_site(bath, bath_shape, indx, save_int)
        else:
            orb_equality_bath_ineq(bath, bath_shape, indx, save_int)
    return bath


# ph_symmetrize_bath
def ph_symmetrize_bath(self, bath, save):
    """
       This function enforces particle-hole symmetry of the bath hybridization \
       function. Not compatible with :code:`REPLICA` or :code:`GENERAL` bath types.

       :type bath: np.array(dtype=float)
       :param bath: The user-accessible bath array
          
       :type save: bool
       :param save: whether to save the symmetry-broken bath for reading
       
       :return: the modified bath array
       :rtype: np.array(dtype=float) 
    """
    ph_symmetrize_bath_site = self.library.ph_symmetrize_bath_site
    ph_symmetrize_bath_site.argtypes = [
        np.ctypeslib.ndpointer(dtype=float, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        c_int,
    ]
    ph_symmetrize_bath_site.restypes = None

    ph_symmetrize_bath_ineq = self.library.ph_symmetrize_bath_ineq
    ph_symmetrize_bath_ineq.argtypes = [
        np.ctypeslib.ndpointer(dtype=float, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        c_int,
    ]
    ph_symmetrize_bath_ineq.restypes = None
    if save:
        save_int = 1
    else:
        save_int = 0
    bath_shape = np.asarray(np.shape(bath), dtype=np.int64, order="F")
    if (len(bath_shape)) == 1:
        ph_symmetrize_bath_site(bath, bath_shape, save_int)
    else:
        ph_symmetrize_bath_ineq(bath, bath_shape, save_int)
    return bath


# save array as .restart file
def save_array_as_bath(self, bath):
    """
       This function takes the user-accessible array and saves it in the \
       correct format for every bath type in the file :code:`hamiltonian.restart`

       :type bath: np.array(dtype=float)
       :param bath: The user-accessible bath array
       
       :return: Nothing
       :rtype: None
    """
    save_array_as_bath_site = self.library.save_array_as_bath_site
    save_array_as_bath_site.argtypes = [
        np.ctypeslib.ndpointer(dtype=float, ndim=1, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
    ]
    save_array_as_bath_site.restypes = None

    save_array_as_bath_ineq = self.library.save_array_as_bath_ineq
    save_array_as_bath_ineq.argtypes = [
        np.ctypeslib.ndpointer(dtype=float, ndim=2, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
    ]
    save_array_as_bath_ineq.restypes = None

    bath_shape = np.asarray(np.shape(bath), dtype=np.int64, order="F")
    if (len(bath_shape)) == 1:
        save_array_as_bath_site(bath, bath_shape)
    else:
        save_array_as_bath_ineq(bath, bath_shape)
    return


# auxiliary functions to get/set bath structure. Only works for single-site. User has to do a loop on sites


def bath_inspect(self, bath=None, e=None, v=None, d=None, u=None, l=None):
    """
       This function translates between the user-accessible continuous \
       bath array and the bath components (energy level, hybridization and so on). \
       It functions in both ways, given the array returns the components and \
       vice-versa. It autonomously determines the type of bath and ED mode.

       :type bath: np.array(dtype=float)
       :param bath: The user-accessible bath array
       
       :type e: np.array(dtype=float)
       :param e: an array for the bath levels ( :f:var:`ed_mode` = \
        :code:`NORMAL, NONSU2, SUPERC`). It has dimension [ :data:`Nspin` , \
        :data:`Norb` ,  :data:`Nbath` ] for :code:`NORMAL` bath, \
        [ :data:`Nspin` ,  :data:`Nbath` ] for :code:`HYBRID` bath 
       
       :type v: np.array(dtype=float)
       :param v: an array for the bath hybridizations ( :f:var:`ed_mode` = \
        :code:`NORMAL, NONSU2, SUPERC`). It has dimension [ :data:`Nspin` , :data:`Norb` , \
        :data:`Nbath` ] for :code:`NORMAL` and :code:`HYBRID` bath. \
        For :code:`REPLICA` bath it has dimension [ :data:`Nbath` ] and for \
        :code:`GENERAL` bath it has dimension [ :data:`Nbath` , :data:`Nspin` \
        :math:`\\cdot` :data:`Norb` ]
       
       :type d: np.array(dtype=float)
       :param d: an array for the bath anomalous enery levels( :f:var:`ed_mode` \
        = :code:`SUPERC`). It has dimension [ :data:`Nspin` ,  :data:`Norb` , \
        :data:`Nbath` ] for :code:`NORMAL` bath, [ :data:`Nspin` , :data:`Nbath` ] \
        for :code:`HYBRID` bath
       
       :type u: np.array(dtype=float)
       :param u: an array for the bath spin off-diagonal hybridization \
        ( :f:var:`ed_mode` = :code:`NONSU2`). It has dimension [ :data:`Nspin`, \
        :data:`Norb` ,  :data:`Nbath` ] for :code:`NORMAL` and :code:`HYBRID` bath

       :type l: np.array(dtype=float)
       :param l: an array for the linear coefficients of the Replica matrix \
        linear combination ( :f:var:`bath_type` = :code:`REPLICA,GENERAL`). \
        It has dimension [ :data:`Nbath` , :code:`Nsym` ], the latter being \
        the number of terms on the linear combination

       :raise ValueError: if both :code:`bath` and some among :code:`e,u,v,d,l` \
        are provided, none is provided, the shapes are inconsistent \
        or the inputs are inconsistent with :f:var:`bath_type` and :f:var:`ed_mode` .

       :return: 
         - if :code:`bath` is provided, returns :code:`e,v`, :code:`e,d,v`, \
           :code:`e,v,u` or :code:`l,v` depending on :f:var:`ed_mode`
         - if :code:`e,v`, :code:`e,d,v`, :code:`e,v,u` or :code:`l,v` \
           depending on :f:var:`ed_mode`  are provided, returns :code:`bath` 
       :rtype: np.array(dtype=float) 
    """

    aux_norb = c_int.in_dll(self.library, "Norb").value
    aux_nspin = c_int.in_dll(self.library, "Nspin").value
    aux_nbath = c_int.in_dll(self.library, "Nbath").value

    settings = (self.get_ed_mode(), self.get_bath_type())
    if settings == (1, 1):  # normal ed mode, normal bath
        if bath is None and e is not None and v is not None:
            e = np.asarray(e, order="F")
            v = np.asarray(v, order="F")
            try:
                if np.shape(e) != (aux_nspin, aux_norb, aux_nbath):
                    raise ValueError("e must be (nspin,norb,nbath)")
                if np.shape(v) != (aux_nspin, aux_norb, aux_nbath):
                    raise ValueError("v must be (nspin,norb,nbath)")
            except:
                print(np.shape(e))
                print(np.shape(v))
                raise ValueError("e or v have wrong dimension")

            Nb = self.get_bath_dimension()
            bath = np.zeros(Nb)

            stride = 0
            io = 0
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = (
                            stride
                            + ibath
                            + (iorb) * aux_nbath
                            + (ispin) * aux_nbath * aux_norb
                        )
                        bath[io] = e[ispin, iorb, ibath]
            stride = aux_nspin * aux_norb * aux_nbath
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = (
                            stride
                            + ibath
                            + (iorb) * aux_nbath
                            + (ispin) * aux_nbath * aux_norb
                        )
                        bath[io] = v[ispin, iorb, ibath]
            return bath

        elif bath is not None and e is None and v is None:  # e and v are none
            bath = np.asarray(bath, order="F")
            Nb = self.get_bath_dimension()
            if np.shape(bath)[0] != Nb:
                raise ValueError("bath has the wrong length")

            e = np.zeros((aux_nspin, aux_norb, aux_nbath))
            v = np.zeros((aux_nspin, aux_norb, aux_nbath))

            stride = 0
            io = 0
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = (
                            stride
                            + ibath
                            + (iorb) * aux_nbath
                            + (ispin) * aux_nbath * aux_norb
                        )
                        e[ispin, iorb, ibath] = bath[io]
            stride = aux_nspin * aux_norb * aux_nbath
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = (
                            stride
                            + ibath
                            + (iorb) * aux_nbath
                            + (ispin) * aux_nbath * aux_norb
                        )
                        v[ispin, iorb, ibath] = bath[io]
            return e, v
        else:
            raise ValueError("Wrong input for normal/normal")

    elif settings == (2, 1):  # superc ed mode, normal bath
        if bath is None and e is not None and v is not None and d is not None:
            e = np.asarray(e, order="F")
            v = np.asarray(v, order="F")
            d = np.asarray(u, order="F")
            try:
                if np.shape(e) != (aux_nspin, aux_norb, aux_nbath):
                    raise ValueError("e must be (nspin,norb,nbath)")
                if np.shape(d) != (aux_nspin, aux_norb, aux_nbath):
                    raise ValueError("d must be (nspin,norb,nbath)")
                if np.shape(v) != (aux_nspin, aux_norb, aux_nbath):
                    raise ValueError("v must be (nspin,norb,nbath)")
            except:
                raise ValueError("e,d or v have wrong dimension")

            Nb = self.get_bath_dimension()
            bath = np.zeros(Nb)
            stride = 0
            io = 0
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = (
                            stride
                            + ibath
                            + (iorb) * aux_nbath
                            + (ispin) * aux_nbath * aux_norb
                        )
                        bath[io] = e[ispin, iorb, ibath]
            stride = aux_nspin * aux_norb * aux_nbath
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = (
                            stride
                            + ibath
                            + (iorb) * aux_nbath
                            + (ispin) * aux_nbath * aux_norb
                        )
                        bath[io] = d[ispin, iorb, ibath]
            stride = 2 * aux_nspin * aux_norb * aux_nbath
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = (
                            stride
                            + ibath
                            + (iorb) * aux_nbath
                            + (ispin) * aux_nbath * aux_norb
                        )
                        bath[io] = v[ispin, iorb, ibath]
            return bath
        elif bath is None and e is None and v is None and d is None:
            bath = np.asarray(bath, order="F")
            Nb = self.get_bath_dimension()
            if np.shape(bath)[0] != Nb:
                raise ValueError("bath has the wrong length")

            e = np.zeros((aux_nspin, aux_norb, aux_nbath))
            v = np.zeros((aux_nspin, aux_norb, aux_nbath))
            d = np.zeros((aux_nspin, aux_norb, aux_nbath))

            stride = 0
            io = 0
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = (
                            stride
                            + ibath
                            + (iorb) * aux_nbath
                            + (ispin) * aux_nbath * aux_norb
                        )
                        e[ispin, iorb, ibath] = bath[io]
            stride = aux_nspin * aux_norb * aux_nbath
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = (
                            stride
                            + ibath
                            + (iorb) * aux_nbath
                            + (ispin) * aux_nbath * aux_norb
                        )
                        d[ispin, iorb, ibath] = bath[io]
            stride = 2 * aux_nspin * aux_norb * aux_nbath
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = (
                            stride
                            + ibath
                            + (iorb) * aux_nbath
                            + (ispin) * aux_nbath * aux_norb
                        )
                        v[ispin, iorb, ibath] = bath[io]
            return e, d, v
        else:
            raise ValueError("Wrong input for superc/normal")

    elif settings == (3, 1):  # nonsu2 ed mode, normal bath
        if bath is None and e is not None and v is not None and u is not None:
            try:
                e = np.asarray(e, order="F")
                v = np.asarray(v, order="F")
                u = np.asarray(u, order="F")
                if np.shape(e) != (aux_nspin, aux_norb, aux_nbath):
                    raise ValueError("e must be (nspin,norb,nbath)")
                if np.shape(v) != (aux_nspin, aux_norb, aux_nbath):
                    raise ValueError("v must be (nspin,norb,nbath)")
                if np.shape(u) != (aux_nspin, aux_norb, aux_nbath):
                    raise ValueError("u must be (nspin,norb,nbath)")
            except:
                raise ValueError("e,v or u have wrong dimension")

            Nb = self.get_bath_dimension()
            bath = np.zeros(Nb)
            stride = 0
            io = 0
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = (
                            stride
                            + ibath
                            + (iorb) * aux_nbath
                            + (ispin) * aux_nbath * aux_norb
                        )
                        bath[io] = e[ispin, iorb, ibath]
            stride = aux_nspin * aux_norb * aux_nbath
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = (
                            stride
                            + ibath
                            + (iorb) * aux_nbath
                            + (ispin) * aux_nbath * aux_norb
                        )
                        bath[io] = v[ispin, iorb, ibath]
            stride = 2 * aux_nspin * aux_norb * aux_nbath
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = (
                            stride
                            + ibath
                            + (iorb) * aux_nbath
                            + (ispin) * aux_nbath * aux_norb
                        )
                        bath[io] = u[ispin, iorb, ibath]
            return bath
        elif bath is not None and e is None and v is None and u is None:
            bath = np.asarray(bath, order="F")
            Nb = self.get_bath_dimension()
            if np.shape(bath)[0] != Nb:
                raise ValueError("bath has the wrong length")

            e = np.zeros((aux_nspin, aux_norb, aux_nbath))
            v = np.zeros((aux_nspin, aux_norb, aux_nbath))
            u = np.zeros((aux_nspin, aux_norb, aux_nbath))

            stride = 0
            io = 0
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = (
                            stride
                            + ibath
                            + (iorb) * aux_nbath
                            + (ispin) * aux_nbath * aux_norb
                        )
                        e[ispin, iorb, ibath] = bath[io]
            stride = aux_nspin * aux_norb * aux_nbath
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = (
                            stride
                            + ibath
                            + (iorb) * aux_nbath
                            + (ispin) * aux_nbath * aux_norb
                        )
                        v[ispin, iorb, ibath] = bath[io]
            stride = 2 * aux_nspin * aux_norb * aux_nbath
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = (
                            stride
                            + ibath
                            + (iorb) * aux_nbath
                            + (ispin) * aux_nbath * aux_norb
                        )
                        u[ispin, iorb, ibath] = bath[io]
            return e, v, u
        else:
            raise ValueError("Wrong input for nonsu2/normal")

    elif settings == (1, 2):  # normal ed mode, hybrid bath
        if bath is None and e is not None and v is not None:
            try:
                e = np.asarray(e, order="F")
                v = np.asarray(v, order="F")
                if np.shape(e) != (aux_nspin, aux_nbath):
                    raise ValueError("e must be (nspin,nbath)")
                if np.shape(v) != (aux_nspin, aux_norb, aux_nbath):
                    raise ValueError("v must be (nspin,norb,nbath)")
            except:
                print(np.shape(e))
                print(np.shape(v))
                raise ValueError("e or v have wrong dimension")

            Nb = self.get_bath_dimension()
            bath = np.zeros(Nb)

            stride = 0
            io = 0
            for ispin in range(aux_nspin):
                for ibath in range(aux_nbath):
                    io = stride + ibath + (ispin) * aux_nbath
                    bath[io] = e[ispin, ibath]
            stride = aux_nspin * aux_nbath
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = (
                            stride
                            + ibath
                            + (iorb) * aux_nbath
                            + (ispin) * aux_nbath * aux_norb
                        )
                        bath[io] = v[ispin, iorb, ibath]
            return bath

        elif bath is not None and e is None and v is None:
            bath = np.asarray(bath, order="F")
            Nb = self.get_bath_dimension()
            if np.shape(bath)[0] != Nb:
                raise ValueError("bath has the wrong length")

            e = np.zeros((aux_nspin, aux_nbath))
            v = np.zeros((aux_nspin, aux_norb, aux_nbath))

            stride = 0
            io = 0
            for ispin in range(aux_nspin):
                for ibath in range(aux_nbath):
                    io = stride + ibath + ispin * aux_nbath
                    e[ispin, ibath] = bath[io]
            stride = aux_nspin * aux_nbath
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = (
                            stride
                            + ibath
                            + (iorb) * aux_nbath
                            + (ispin) * aux_nbath * aux_norb
                        )
                        v[ispin, iorb, ibath] = bath[io]
            return e, v
        else:
            raise ValueError("Wrong input for normal/hybrid")

    elif settings == (2, 2):  # superc ed mode, hybrid bath
        if bath is None and e is not None and v is not None and d is not None:
            try:
                e = np.asarray(e, order="F")
                d = np.asarray(d, order="F")
                v = np.asarray(v, order="F")
                if np.shape(e) != (aux_nspin, aux_nbath):
                    raise ValueError("e must be (nspin,nbath)")
                if np.shape(d) != (aux_nspin, aux_nbath):
                    raise ValueError("d must be (nspin,nbath)")
                if np.shape(v) != (aux_nspin, aux_norb, aux_nbath):
                    raise ValueError("v must be (nspin,norb,nbath)")
            except:
                raise ValueError("e,d or v have wrong dimension")

            Nb = self.get_bath_dimension()
            bath = np.zeros(Nb)
            stride = 0
            io = 0
            for ispin in range(aux_nspin):
                for ibath in range(aux_nbath):
                    io = stride + ibath + (ispin) * aux_nbath
                    bath[io] = e[ispin, ibath]
            stride = aux_nspin * aux_nbath
            for ispin in range(aux_nspin):
                for ibath in range(aux_nbath):
                    io = stride + ibath + (ispin) * aux_nbath
                    bath[io] = bath_d[ispin, ibath]
            stride = 2 * aux_nspin * aux_nbath
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = (
                            stride
                            + ibath
                            + (iorb) * aux_nbath
                            + (ispin) * aux_nbath * aux_norb
                        )
                        bath[io] = v[ispin, iorb, ibath]
            return bath
        elif bath is not None and e is None and v is None and d is None:
            bath = np.asarray(bath, order="F")
            Nb = self.get_bath_dimension()
            if np.shape(bath)[0] != Nb:
                raise ValueError("bath has the wrong length")

            e = np.zeros((aux_nspin, aux_nbath))
            d = np.zeros((aux_nspin, aux_nbath))
            v = np.zeros((aux_nspin, aux_norb, aux_nbath))

            stride = 0
            io = 0
            for ispin in range(aux_nspin):
                for ibath in range(aux_nbath):
                    io = stride + ibath + (ispin) * aux_nbath
                    e[ispin, ibath] = bath[io]
            stride = aux_nspin * aux_nbath
            for ispin in range(aux_nspin):
                for ibath in range(aux_nbath):
                    io = stride + ibath + (ispin) * aux_nbath
                    d[ispin, ibath] = bath[io]
            stride = 2 * aux_nspin * aux_nbath
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = (
                            stride
                            + ibath
                            + (iorb) * aux_nbath
                            + (ispin) * aux_nbath * aux_norb
                        )
                        v[ispin, iorb, ibath] = bath[io]
            return e, d, v
        else:
            raise ValueError("Wrong input for superc/hybrid")

    elif settings == (3, 2):  # nonsu2 ed mode, hybrid bath
        if bath is None and e is not None and v is not None and u is not None:
            try:
                e = np.asarray(e, order="F")
                v = np.asarray(v, order="F")
                u = np.asarray(u, order="F")
                if np.shape(e) != (aux_nspin, aux_nbath):
                    raise ValueError("e must be (nspin,norb,nbath)")
                if np.shape(v) != (aux_nspin, aux_norb, aux_nbath):
                    raise ValueError("v must be (nspin,norb,nbath)")
                if np.shape(u) != (aux_nspin, aux_norb, aux_nbath):
                    raise ValueError("u must be (nspin,norb,nbath)")
            except:
                raise ValueError("e,v or u have wrong dimension")

            Nb = self.get_bath_dimension()
            bath = np.zeros(Nb)

            stride = 0
            io = 0
            for ispin in range(aux_nspin):
                for ibath in range(aux_nbath):
                    io = stride + ibath + (ispin) * aux_nbath
                    bat[io] = e[ispin, ibath]
            stride = aux_nspin * aux_nbath
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = (
                            stride
                            + ibath
                            + (iorb) * aux_nbath
                            + (ispin) * aux_nbath * aux_norb
                        )
                        bath[io] = v[ispin, iorb, ibath]
            stride = aux_nspin * aux_nbath + aux_nspin * aux_norb * aux_nbath
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = (
                            stride
                            + ibath
                            + (iorb) * aux_nbath
                            + (ispin) * aux_nbath * aux_norb
                        )
                        bath[io] = u[ispin, iorb, ibath]
            return bath
        elif bath is not None and e is None and v is None and u is None:
            bath = np.asarray(bath, order="F")
            Nb = self.get_bath_dimension()
            if np.shape(bath)[0] != Nb:
                raise ValueError("bath has the wrong length")

            e = np.zeros((aux_nspin, aux_nbath))
            v = np.zeros((aux_nspin, aux_norb, aux_nbath))
            u = np.zeros((aux_nspin, aux_norb, aux_nbath))

            stride = 0
            io = 0
            for ispin in range(aux_nspin):
                for ibath in range(aux_nbath):
                    io = stride + ibath + (ispin) * aux_nbath
                    e[ispin, ibath] = bath[io]
            stride = aux_nspin * aux_nbath
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = (
                            stride
                            + ibath
                            + (iorb) * aux_nbath
                            + (ispin) * aux_nbath * aux_norb
                        )
                        v[ispin, iorb, ibath] = bath[io]
            stride = aux_nspin * aux_nbath + aux_nspin * aux_norb * aux_nbath
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = (
                            stride
                            + ibath
                            + (iorb) * aux_nbath
                            + (ispin) * aux_nbath * aux_norb
                        )
                        u[ispin, iorb, ibath] = bath[io]
            return e, v, u
        else:
            raise ValueError("Wrong input for nonsu2/hybrid")
    elif settings == (1, 3) or settings == (2, 3) or settings == (3, 3):  # replica bath
        if bath is None and l is not None and v is not None:
            l = np.asarray(l, order="F")
            v = np.asarray(v, order="F")
            if self.Nsym is None:
                raise ValueError("Nsym is none, is Hreplica initialized?")
            if np.shape(l) != (aux_nbath, self.Nsym):
                raise ValueError("l must be (nbath,nsym)")
            if np.shape(v) != (aux_nbath):
                raise ValueError("v must be (nbath)")

            Nb = self.get_bath_dimension()
            bath = np.zeros(Nb)

            io = 0
            il = 0

            bath[io] = self.Nsym
            io += 1

            for ibath in range(aux_nbath):
                bath[io] = v[ibath]
                io += 1
                il = 0
                for il in range(Nsym):
                    bath[io] = l[ibath, il]
                    io += 1
                    il += 1
            return bath

        elif bath is not None and l is None and v is None:  # e and v are none
            bath = np.asarray(bath, order="F")
            Nb = self.get_bath_dimension()
            if np.shape(bath)[0] != Nb:
                raise ValueError("bath has the wrong length")
            if bath[0] != self.Nsym:
                raise ValueError("bath[0] is not Nsym")

            l = np.zeros((aux_nbath, self.Nsym))
            v = np.zeros((aux_nbath))

            io = 1
            il = 0

            for ibath in range(aux_nbath):
                v[ibath] = bath[io]
                io += 1
                il = 0
                for il in range(Nsym):
                    l[ibath, il] = bath[io]
                    io += 1
                    il += 1
            return l, v
        else:
            raise ValueError("Wrong input for replica")

    elif settings == (1, 4) or settings == (2, 4) or settings == (3, 4):  # general bath
        if bath is None and l is not None and v is not None:
            l = np.asarray(l, order="F")
            v = np.asarray(v, order="F")
            if self.Nsym is None:
                raise ValueError("Nsym is none, is Hgeneral initialized?")
            if np.shape(l) != (aux_nbath, self.Nsym):
                raise ValueError("l must be (Nbath,Nsym)")
            if np.shape(v) != (aux_nbath, aux_nspin * aux_norb):
                raise ValueError("v must be (Nbath,Nspin*Norb)")

            Nb = self.get_bath_dimension()
            bath = np.zeros(Nb)

            io = 0
            il = 0
            iv = 0

            bath[io] = self.Nsym
            io += 1

            for ibath in range(aux_nbath):
                iv = 0
                for iv in range(aux_nspin * aux_norb):
                    bath[io] = v[ibath, iv]
                    io += 1
                    iv += 1
                il = 0
                for il in range(Nsym):
                    bath[io] = l[ibath, il]
                    io += 1
                    il += 1
            return bath

        elif bath is not None and l is None and v is None:  # e and v are none
            bath = np.asarray(bath, order="F")
            Nb = self.get_bath_dimension()
            if np.shape(bath)[0] != Nb:
                raise ValueError("bath has the wrong length")
            if bath[0] != self.Nsym:
                raise ValueError("bath[0] is not Nsym")

            l = np.zeros((aux_nbath, self.Nsym))
            v = np.zeros((aux_nbath, aux_nspin * aux_norb))

            io = 1
            il = 0
            iv = 0

            for ibath in range(aux_nbath):
                iv = 0
                for iv in range(aux_nspin * aux_norb):
                    v[ibath, iv] = bath[io]
                    io += 1
                    iv += 1
                il = 0
                for il in range(Nsym):
                    l[ibath, il] = bath[io]
                    io += 1
                    il += 1
            return l, v
        else:
            raise ValueError("Wrong input for replica")

    else:
        raise ValueError("EDmode/bath combination not valid or not implemented.")
