from ctypes import *
import numpy as np
import os, sys
import types

# observables


# density
def get_dens(self, ilat=None, iorb=None):
    """
    
    This function returns the value of the charge density
   
    :type ilat: int
    :param ilat: if the case of real-space DMFT, if only the Green's function of \
    a specific inequivalent site is needed, this can be specified.
   
    :type iorb: int
    :param iorb: the orbital index. If none is provided, the whole density \
     vector is returned
   
    :return: the full charge density tensor has dimensions [ :code:`Nlat` ,Norb]. \
     Depending on which keyword arguments are (or not) provided, this is sliced \
     on the corresponding axis.
    :rtype: float **or** np.array(dtype=float) 
    
    """

    aux_norb = c_int.in_dll(self.library, "Norb").value

    ed_get_dens_n1_wrap = self.library.ed_get_dens_n1
    ed_get_dens_n1_wrap.argtypes = [
        np.ctypeslib.ndpointer(dtype=float, ndim=1, flags="F_CONTIGUOUS")
    ]
    ed_get_dens_n1_wrap.restype = None

    ed_get_dens_n2_wrap = self.library.ed_get_dens_n2
    ed_get_dens_n2_wrap.argtypes = [
        np.ctypeslib.ndpointer(dtype=float, ndim=2, flags="F_CONTIGUOUS"),
        c_int,
    ]
    ed_get_dens_n2_wrap.restype = None

    if self.Nineq == 0:
        densvec = np.zeros(aux_norb, dtype=float, order="F")
        ed_get_dens_n1_wrap(densvec)

        if ilat is not None:
            raise ValueError("ilat cannot be none for single-impurity DMFT")
        elif iorb is not None:
            return densvec[iorb]
        else:
            return densvec
    else:
        densvec = np.zeros([self.Nineq, aux_norb], dtype=float, order="F")
        ed_get_dens_n2_wrap(densvec, self.Nineq)
        densvec = np.asarray(densvec)

        if ilat is not None and iorb is not None:
            return densvec[ilat, iorb]
        elif ilat is None and iorb is not None:
            return densvec[:, iorb]
        elif ilat is not None and iorb is None:
            return densvec[ilat, :]
        else:
            return densvec


# magnetization
def get_mag(self, icomp=None, ilat=None, iorb=None):
    """

       
       This function returns the value of the magnetization
      
       :type icomp: str
       :param icomp: the component of the magnetization, :code:`"x"`, \
       :code:`"y"` or :code:`"z"` (default).
       
       :type ilat: int
       :param ilat: if the case of real-space DMFT, if only the Green's function \
       of a specific inequivalent site is needed, this can be specified.
       
       :type iorb: int
       :param iorb: the orbital index. If none is provided, the whole density \
       vector is returned
       
       :return: the full magnetization tensor has dimensions [ :code:`Nlat` ,3,Norb]. \
        Depending on which keyword arguments are (or not) provided, this is \
        sliced on the corresponding axis.
       :rtype: float **or** np.array(dtype=float) 
       
     """

    if icomp == "x" or icomp == "X":
        icomp = 0
    elif icomp == "y" or icomp == "Y":
        icomp = 1
    elif icomp == "z" or icomp == "Z":
        icomp = 2

    aux_norb = c_int.in_dll(self.library, "Norb").value

    ed_get_mag_n2_wrap = self.library.ed_get_mag_n2
    ed_get_mag_n2_wrap.argtypes = [
        np.ctypeslib.ndpointer(dtype=float, ndim=2, flags="F_CONTIGUOUS")
    ]
    ed_get_mag_n2_wrap.restype = None

    ed_get_mag_n3_wrap = self.library.ed_get_mag_n3
    ed_get_mag_n3_wrap.argtypes = [
        np.ctypeslib.ndpointer(dtype=float, ndim=3, flags="F_CONTIGUOUS"),
        c_int,
    ]
    ed_get_mag_n3_wrap.restype = None

    if self.Nineq == 0:
        magvec = np.zeros([3, aux_norb], dtype=float, order="F")
        ed_get_mag_n2_wrap(magvec)

        if ilat is not None:
            raise ValueError("ilat cannot be none for single-impurity DMFT")
        elif iorb is not None and icomp is not None:
            return magvec[icomp, iorb]
        elif iorb is not None and icomp is None:
            return magvec[:, iorb]
        elif iorb is None and icomp is not None:
            return magvec[icomp, :]
        elif iorb is None and icomp is None:
            return magvec
    else:
        magvec = np.zeros([self.Nineq, 3, aux_norb], dtype=float, order="F")
        ed_get_mag_n3_wrap(magvec, self.Nineq)
        magvec = np.asarray(magvec)

        if ilat is not None:
            if iorb is not None and icomp is not None:
                return magvec[ilat, icomp, iorb]
            if iorb is None and icomp is not None:
                return magvec[ilat, icomp, :]
            if iorb is not None and icomp is None:
                return magvec[ilat, :, iorb]
            if iorb is None and icomp is None:
                return magvec[ilat, :, :]
        else:
            if iorb is not None and icomp is not None:
                return magvec[:, icomp, iorb]
            if iorb is None and icomp is not None:
                return magvec[:, icomp, :]
            if iorb is not None and icomp is None:
                return magvec[:, :, iorb]
            if iorb is None and icomp is None:
                return magvec


# double occupation
def get_docc(self, ilat=None, iorb=None):
    """
   This function returns the value of the double occupation
  
   :type ilat: int
   :param ilat: if the case of real-space DMFT, if only the Green's function of a \
    specific inequivalent site is needed, this can be specified.
   
   :type iorb: int
   :param iorb: the orbital index. If none is provided, the whole density vector \
    is returned
   
   :return: the full double-occupation tensor has dimensions [ :code:`Nlat` ,Norb]. \
    Depending on which keyword arguments are (or not) provided, this is sliced \
    on the corresponding axis.
   :rtype: float **or** np.array(dtype=float)
   
   """

    aux_norb = c_int.in_dll(self.library, "Norb").value

    ed_get_docc_n1_wrap = self.library.ed_get_docc_n1
    ed_get_docc_n1_wrap.argtypes = [
        np.ctypeslib.ndpointer(dtype=float, ndim=1, flags="F_CONTIGUOUS")
    ]
    ed_get_docc_n1_wrap.restype = None

    ed_get_docc_n2_wrap = self.library.ed_get_docc_n2
    ed_get_docc_n2_wrap.argtypes = [
        np.ctypeslib.ndpointer(dtype=float, ndim=2, flags="F_CONTIGUOUS"),
        c_int,
    ]
    ed_get_docc_n2_wrap.restype = None

    if self.Nineq == 0:
        doccvec = np.zeros(aux_norb, dtype=float, order="F")
        ed_get_docc_n1_wrap(doccvec)

        if ilat is not None:
            raise ValueError("ilat cannot be none for single-impurity DMFT")
        elif iorb is not None:
            return doccvec[iorb]
        else:
            return doccvec
    else:
        doccvec = np.zeros([self.Nineq, aux_norb], dtype=float, order="F")
        ed_get_docc_n2_wrap(doccvec, self.Nineq)
        doccvec = np.asarray(doccvec)

        if ilat is not None and iorb is not None:
            return doccvec[ilat, iorb]
        elif ilat is None and iorb is not None:
            return doccvec[:, iorb]
        elif ilat is not None and iorb is None:
            return doccvec[ilat, :]
        else:
            return doccvec


# superconductive phi
def get_phi(self, ilat=None, iorb=None, jorb=None):
    """
    This function returns the value of the superconductive order \
    parameter :math:`\\phi = \\langle c_{\\uparrow} c_{\\downarrow} \\rangle`
  
   :type ilat: int
   :param ilat: if the case of real-space DMFT, if only the Green's function of a \
    specific inequivalent site is needed, this can be specified.
   
   :type iorb: int
   :param iorb: the first orbital index
   
   :type jorb: int
   :param jorb: the second orbital index
   
   :return: :math:`\\phi`. The full tensor has dimensions [ :code:`Nlat`, \
   :data:`Norb`, :data:`Norb`]. Depending on which keyword arguments are \
    (or not) provided, this is sliced on the corresponding axis.
   :rtype: float **or** np.array(dtype=float)
   
   """

    aux_norb = c_int.in_dll(self.library, "Norb").value

    ed_get_phisc_n2_wrap = self.library.ed_get_phisc_n2
    ed_get_phisc_n2_wrap.argtypes = [
        np.ctypeslib.ndpointer(dtype=float, ndim=2, flags="F_CONTIGUOUS")  # self
    ]
    ed_get_phisc_n2_wrap.restype = None

    ed_get_phisc_n3_wrap = self.library.ed_get_phisc_n3
    ed_get_phisc_n3_wrap.argtypes = [
        np.ctypeslib.ndpointer(dtype=float, ndim=3, flags="F_CONTIGUOUS"),  # self
        c_int,  # Nlat
    ]
    ed_get_phisc_n3_wrap.restype = None

    if self.Nineq == 0:
        phivec = np.zeros((aux_norb, aux_norb), dtype=float, order="F")
        ed_get_phisc_n2_wrap(phivec)
        phivec = np.asarray(phivec)

        if ilat is not None:
            raise ValueError("ilat cannot be none for single-impurity DMFT")
        elif iorb is not None and jorb is not None:
            return phivec[iorb, jorb]
        elif iorb is not None and jorb is None:
            return phivec[iorb, :]
        elif jorb is not None and iorb is None:
            return phivec[:, jorb]
        else:
            return phivec
    else:
        phivec = np.zeros([self.Nineq, aux_norb, aux_norb], dtype=float, order="F")
        ed_get_docc_n3_wrap(phivec, self.Nineq)
        phivec = np.asarray(phivec)

        if ilat is not None:
            if iorb is not None and jorb is not None:
                return phivec[ilat, iorb, jorb]
            if iorb is not None and jorb is None:
                return phivec[ilat, iorb, :]
            if jorb is not None and iorb is None:
                return phivec[ilat, :, jorb]
            else:
                return phivec[ilat, :, :]
        else:
            if iorb is not None and jorb is not None:
                return phivec[:, iorb, jorb]
            if iorb is not None and jorb is None:
                return phivec[:, iorb, :]
            if jorb is not None and iorb is None:
                return phivec[:, :, jorb]
            else:
                return phivec


# energy
def get_eimp(self, ilat=None, ikind=None):
    """
       This function returns the value of the local energy components
         
       :type ilat: int
       :param ilat: if the case of real-space DMFT, if only the Green's function of \
       a specific inequivalent site is needed, this can be specified.
       
       :type ikind: int
       :param ikind: index of the component. It is
        
        * :code:`1`: ed_Epot: the potential energy from interaction
        * :code:`2`: ed_Eint: ed-Epot - ed_Ehartree (? it is not assigned) 
        * :code:`3`: ed_Ehartree: Hartree part of interaction energy
        * :code:`4`: ed_Eknot: on-site part of the kinetic term
       
       :return: the full local energy tensor has dimensions [ :code:`Nlat` ,4]. Depending on \
        which keyword arguments are (or not) provided, this is sliced on the corresponding axis.
       :rtype: float **or** np.array(dtype=float)    
    """

    ed_get_eimp_n1_wrap = self.library.ed_get_eimp_n1
    ed_get_eimp_n1_wrap.argtypes = [
        np.ctypeslib.ndpointer(dtype=float, ndim=1, flags="F_CONTIGUOUS")
    ]
    ed_get_eimp_n1_wrap.restype = None

    ed_get_eimp_n2_wrap = self.library.ed_get_eimp_n2
    ed_get_eimp_n2_wrap.argtypes = [
        np.ctypeslib.ndpointer(dtype=float, ndim=2, flags="F_CONTIGUOUS"),
        c_int,
    ]
    ed_get_eimp_n2_wrap.restype = None

    if self.Nineq == 0:
        eimp_vec = np.zeros(4, dtype=float, order="F")
        ed_get_eimp_n1_wrap(eimp_vec)

        if ilat is not None:
            raise ValueError("ilat cannot be none for single-impurity DMFT")
        elif ikind is not None:
            return eimp_vec[ikind]
        else:
            return eimp_vec
    else:
        eimp_vec = np.zeros([self.Nineq, 4], dtype=float, order="F")
        ed_get_eimp_n2_wrap(eimp_vec, self.Nineq)
        eimp_vec = np.asarray(eimp_vec)

        if ilat is not None and ikind is not None:
            return eimp_vec[ilat, ikind]
        elif ilat is None and ikind is not None:
            return eimp_vec[:, ikind]
        elif ilat is not None and ikind is None:
            return eimp_vec[ilat, :]
        else:
            return eimp_vec


########################
#   SIGMA              #
########################


# backcompatibility, undocumented
def build_sigma(self, zeta, ilat=None, ishape=None, typ="n"):
    return self.get_sigma(zeta=zeta, ilat=ilat, ishape=ishape, typ=typ)


# get Sigma
def get_sigma(self, ilat=None, ishape=None, axis="m", typ="n", zeta=None):
    """
    This function generates the self-energy for a user-chosen set of frequencies \
    in the complex plane

    :type ilat: int
    :param ilat: if the case of real-space DMFT, if only the self-energy of \
    a specific inequivalent site is needed, this can be specified.
        
    :type ishape: int 
    :param ishape: this variable determines the shape of the returned array. \
     Possible values:
   
     * :code:`None`: the same shape as :code:`Hloc` plus one axis for frequency 
     * :code:`3`: in the single-impurity case, it will return an array of the shape \
       [ :data:`Nspin` :math:`\\cdot`  :data:`Norb` ,  :data:`Nspin` :math:`\\cdot` \
       :data:`Norb` , :code:`len(zeta)` ]. In the real-space DMFT case, \
       it will return an array of the shape \
       [ :code:`Nlat` :math:`\\cdot`  :data:`Nspin` :math:`\\cdot`  :data:`Norb` , \
       :code:`Nlat` :math:`\\cdot` :data:`Nspin` :math:`\\cdot`  :data:`Norb` , \
       :code:`len(zeta)` ]. :code:`Nlat` will be determined from the module 
       by assessing the shape of Hloc. If :code:`ilat` is set, ValueError is returned.
     * :code:`4`: in the real-space DMFT case, it will return an array of the shape \
       [ :code:`Nlat` ,  :data:`Nspin` :math:`\\cdot`  :data:`Norb` ,  :data:`Nspin` \
       :math:`\\cdot`  :data:`Norb` , :code:`len(zeta)` `. :code:`Nlat` will \
       be determined from the module by assessing the shape of Hloc. If :code:`ilat` is \
       set, the output will have one dimension less.
     * :code:`5`: in the single-impurity case, it will return an array of the \
       shape [ :data:`Nspin` ,  :data:`Nspin` ,  :data:`Norb` ,  :data:`Norb` , \
       :code:`len(zeta)` ].
     * :code:`6`: in the real-space DMFT case, it will return an array of the \
       shape [ :code:`Nlat` ,  :data:`Nspin` ,  :data:`Nspin` ,  :data:`Norb` , \
       :data:`Norb` , :code:`len(zeta)` ]. :code:`Nlat` will be determined from \
       the module by assessing the shape of Hloc. If :code:`ilat` is set, \
       the output will have one dimension less.

    :type axis: str 
    :param axis: if :code:`zeta` is not provided, return the self-energy on the \
     Matsubara or Real axis with parameters set in the input file. \
     Can be :code:`m` for Matsubara(default) or :code:`r` for real.
        
    :type typ: str 
    :param typ: whether to return the normal or anomalous self-energy \
     (for the superconducting case). Can be :code:`n` for normal (default)  \
     or :code:`a` for anomalous.
    

    :type zeta: complex **or** [complex] **or** np.array(dtype=complex)
    :param zeta: user-defined array of frequencies in the whole complex plane. \
     If none is provided, according to :code:`axis` the Matsubara or real axis is chosen

   
    :raise ValueError: If :code:`ishape` is incompatible woth :code:`ilat` or \
     not in the previous list.
    :raise ValueError: If :code:`axis` is not in the previous list.
     
    :return: An array of floats that contains the self-energy along the \
     specific axis, with dimension set by :code:`ishape` and :code:`zeta`, if present.  
    :rtype: np.array(dtype=float) 
    """

    ed_get_sigma_site_n3 = self.library.get_sigma_site_n3
    ed_get_sigma_site_n3.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=3, flags="F_CONTIGUOUS"),  # self
        c_int,  # axis
        c_int,  # typ
        np.ctypeslib.ndpointer(dtype=complex, ndim=1, flags="F_CONTIGUOUS"),  # zeta
        c_int,  # dz
        c_int,  # zflag
    ]
    ed_get_sigma_site_n3.restype = None

    ed_get_sigma_site_n5 = self.library.get_sigma_site_n5
    ed_get_sigma_site_n5.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=5, flags="F_CONTIGUOUS"),  # self
        c_int,  # axis
        c_int,  # typ
        np.ctypeslib.ndpointer(dtype=complex, ndim=1, flags="F_CONTIGUOUS"),  # zeta
        c_int,  # dz
        c_int,  # zflag
    ]
    ed_get_sigma_site_n5.restype = None

    ed_get_sigma_lattice_n3 = self.library.get_sigma_lattice_n3
    ed_get_sigma_lattice_n3.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=3, flags="F_CONTIGUOUS"),  # self
        c_int,  # nineq
        c_int,  # axis
        c_int,  # typ
        np.ctypeslib.ndpointer(dtype=complex, ndim=1, flags="F_CONTIGUOUS"),  # zeta
        c_int,  # dz
        c_int,  # zflag
    ]
    ed_get_sigma_lattice_n3.restype = None

    ed_get_sigma_lattice_n4 = self.library.get_sigma_lattice_n4
    ed_get_sigma_lattice_n4.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=4, flags="F_CONTIGUOUS"),  # self
        c_int,  # nineq
        c_int,  # axis
        c_int,  # typ
        np.ctypeslib.ndpointer(dtype=complex, ndim=1, flags="F_CONTIGUOUS"),  # zeta
        c_int,  # dz
        c_int,  # zflag
    ]
    ed_get_sigma_lattice_n4.restype = None

    ed_get_sigma_lattice_n6 = self.library.get_sigma_lattice_n6
    ed_get_sigma_lattice_n6.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=4, flags="F_CONTIGUOUS"),  # self
        c_int,  # nineq
        c_int,  # axis
        c_int,  # typ
        np.ctypeslib.ndpointer(dtype=complex, ndim=1, flags="F_CONTIGUOUS"),  # zeta
        c_int,  # dz
        c_int,  # zflag
    ]
    ed_get_sigma_lattice_n6.restype = None

    # Global vars
    norb_aux = c_int.in_dll(self.library, "Norb").value
    nspin_aux = c_int.in_dll(self.library, "Nspin").value

    # zeta
    if zeta is not None:
        if np.isscalar(zeta):
            zeta = [zeta]
        zeta = np.asarray(zeta, dtype=complex, order="F")
        nfreq = np.shape(zeta)[0]
        zflag = 1
        if any(abs(np.real(zeta)) > 1e-10):
            axis = "r"
    else:
        zeta = np.asarray([0.0], dtype=complex, order="F")
        if axis == "m":
            nfreq = c_int.in_dll(self.library, "Lmats").value
        else:
            nfreq = c_int.in_dll(self.library, "Lreal").value
        zflag = 0

    # ishape
    if ishape is None:
        ishape = self.dim_hloc + 1

    # axis
    if axis == "m":
        axisint = 0
    elif axis == "r":
        axisint = 1
    else:
        raise ValueError("get_sigma: axis can only be 'm' or 'r'")

    # typ
    if typ == "n":
        typint = 0
    elif typ == "a":
        typint = 1
    else:
        raise ValueError("get_sigma: typ can only be 'n' or 'a'")

    if self.Nineq == 0:
        if ilat is not None:
            raise ValueError("ilat is not defined in single-impurity DMFT")
        if ishape == 3:
            Sigma = np.zeros(
                [nspin_aux * norb_aux, nspin_aux * norb_aux, nfreq],
                dtype=complex,
                order="F",
            )
            ed_get_sigma_site_n3(Sigma, axisint, typint, zeta, nfreq, zflag)
        elif ishape == 5:
            Sigma = np.zeros(
                [nspin_aux, nspin_aux, norb_aux, norb_aux, nfreq],
                dtype=complex,
                order="F",
            )
            ed_get_sigma_site_n5(Sigma, axisint, typint, zeta, nfreq, zflag)
        else:
            raise ValueError("Shape(array) != 3,5 in get_sigma_site")
        return Sigma
    else:
        if ishape == 3:
            Sigma = np.zeros(
                [
                    self.Nineq * nspin_aux * norb_aux,
                    self.Nineq * nspin_aux * norb_aux,
                    nfreq,
                ],
                dtype=complex,
                order="F",
            )
            ed_get_sigma_site_n3(Sigma, self.Nineq, axisint, typint, zeta, nfreq, zflag)
        elif ishape == 4:
            Sigma = np.zeros(
                [self.Nineq, nspin_aux * norb_aux, nspin_aux * norb_aux, nfreq],
                dtype=complex,
                order="F",
            )
            ed_get_sigma_site_n4(Sigma, self.Nineq, axisint, typint, zeta, nfreq, zflag)
        elif ishape == 6:
            Sigma = np.zeros(
                [self.Nineq, nspin_aux, nspin_aux, norb_aux, norb_aux, nfreq],
                dtype=complex,
                order="F",
            )
            ed_get_sigma_site_n6(Sigma, self.Nineq, axisint, typint, zeta, nfreq, zflag)
        else:
            raise ValueError("Shape(array) != 3,4,6 in get_sigma_lattice")
        if ilat is not None and ishape != 3:
            return Sigma[ilat]
        else:
            return Sigma


#######################
#   GIMP              #
#######################


# backcompatibility, undocumented
def build_gimp(self, zeta, ilat=None, ishape=None, typ="n"):
    return self.get_gimp(zeta=zeta, ilat=ilat, ishape=ishape, typ=typ)


# get gimp
def get_gimp(self, ilat=None, ishape=None, axis="m", typ="n", zeta=None):
    """
    This function generates the impurity Green's function for a user-chosen set \
     of frequencies in the complex plane

    :type ilat: int
    :param ilat: if the case of real-space DMFT, if only the self-energy of \
    a specific inequivalent site is needed, this can be specified.
        
    :type ishape: int 
    :param ishape: this variable determines the shape of the returned array. \
     Possible values:
   
     * :code:`None`: the same shape as :code:`Hloc` plus one axis for frequency 
     * :code:`3`: in the single-impurity case, it will return an array of the shape \
       [ :data:`Nspin` :math:`\\cdot`  :data:`Norb` ,  :data:`Nspin` :math:`\\cdot` \
       :data:`Norb` , :code:`len(zeta)` ]. In the real-space DMFT case, \
       it will return an array of the shape \
       [ :code:`Nlat` :math:`\\cdot`  :data:`Nspin` :math:`\\cdot`  :data:`Norb` , \
       :code:`Nlat` :math:`\\cdot` :data:`Nspin` :math:`\\cdot`  :data:`Norb` , \
       :code:`len(zeta)` ]. :code:`Nlat` will be determined from the module by \
       assessing the shape of Hloc. If :code:`ilat` is set, ValueError is returned.
     * :code:`4`: in the real-space DMFT case, it will return an array of the shape \
       [ :code:`Nlat` ,  :data:`Nspin` :math:`\\cdot`  :data:`Norb` ,  :data:`Nspin` \
       :math:`\\cdot`  :data:`Norb` , :code:`len(zeta)` `. :code:`Nlat` will \
       be determined from the module by assessing the shape of Hloc. If :code:`ilat`\
       is set, the output will have one dimension less.
     * :code:`5`: in the single-impurity case, it will return an array of the \
       shape [ :data:`Nspin` ,  :data:`Nspin` ,  :data:`Norb` ,  :data:`Norb` , \
       :code:`len(zeta)` ].
     * :code:`6`: in the real-space DMFT case, it will return an array of the \
       shape [ :code:`Nlat` ,  :data:`Nspin` ,  :data:`Nspin` ,  :data:`Norb` , \
       :data:`Norb` , :code:`len(zeta)` ]. \
       :code:`Nlat` will be determined from the module by assessing the shape of Hloc. \
       If :code:`ilat` is set, the output will have one dimension less.

    :type axis: str 
    :param axis: if :code:`zeta` is not provided, return the self-energy on the \
     Matsubara or Real axis with parameters set in the input file. \
     Can be :code:`m` for Matsubara(default) or :code:`r` for real.
        
    :type typ: str 
    :param typ: whether to return the normal or anomalous self-energy \
     (for the superconducting case). Can be :code:`n` for normal (default) or \
    :code:`a` for anomalous.
    

    :type zeta: complex **or** [complex] **or** np.array(dtype=complex)
    :param zeta: user-defined array of frequencies in the whole complex plane. \
     If none is provided, according to :code:`axis` the Matsubara or real axis is chosen

   
    :raise ValueError: If :code:`ishape` is incompatible woth :code:`ilat` or \
     not in the previous list.
    :raise ValueError: If :code:`axis` is not in the previous list.
     
    :return: An array of floats that contains the impurity Green's function along the \
     specific axis, with dimension set by :code:`ishape` and :code:`zeta`, if present.  
    :rtype: np.array(dtype=float) 
    """

    ed_get_gimp_site_n3 = self.library.get_gimp_site_n3
    ed_get_gimp_site_n3.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=3, flags="F_CONTIGUOUS"),  # self
        c_int,  # axis
        c_int,  # typ
        np.ctypeslib.ndpointer(dtype=complex, ndim=1, flags="F_CONTIGUOUS"),  # zeta
        c_int,  # dz
        c_int,  # zflag
    ]
    ed_get_gimp_site_n3.restype = None

    ed_get_gimp_site_n5 = self.library.get_gimp_site_n5
    ed_get_gimp_site_n5.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=5, flags="F_CONTIGUOUS"),  # self
        c_int,  # axis
        c_int,  # typ
        np.ctypeslib.ndpointer(dtype=complex, ndim=1, flags="F_CONTIGUOUS"),  # zeta
        c_int,  # dz
        c_int,  # zflag
    ]
    ed_get_gimp_site_n5.restype = None

    ed_get_gimp_lattice_n3 = self.library.get_gimp_lattice_n3
    ed_get_gimp_lattice_n3.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=3, flags="F_CONTIGUOUS"),  # self
        c_int,  # nineq
        c_int,  # axis
        c_int,  # typ
        np.ctypeslib.ndpointer(dtype=complex, ndim=1, flags="F_CONTIGUOUS"),  # zeta
        c_int,  # dz
        c_int,  # zflag
    ]
    ed_get_gimp_lattice_n3.restype = None

    ed_get_gimp_lattice_n4 = self.library.get_gimp_lattice_n4
    ed_get_gimp_lattice_n4.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=4, flags="F_CONTIGUOUS"),  # self
        c_int,  # nineq
        c_int,  # axis
        c_int,  # typ
        np.ctypeslib.ndpointer(dtype=complex, ndim=1, flags="F_CONTIGUOUS"),  # zeta
        c_int,  # dz
        c_int,  # zflag
    ]
    ed_get_gimp_lattice_n4.restype = None

    ed_get_gimp_lattice_n6 = self.library.get_gimp_lattice_n6
    ed_get_gimp_lattice_n6.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=4, flags="F_CONTIGUOUS"),  # self
        c_int,  # nineq
        c_int,  # axis
        c_int,  # typ
        np.ctypeslib.ndpointer(dtype=complex, ndim=1, flags="F_CONTIGUOUS"),  # zeta
        c_int,  # dz
        c_int,  # zflag
    ]
    ed_get_gimp_lattice_n6.restype = None

    # Global vars
    norb_aux = c_int.in_dll(self.library, "Norb").value
    nspin_aux = c_int.in_dll(self.library, "Nspin").value

    # zeta
    if zeta is not None:
        if np.isscalar(zeta):
            zeta = [zeta]
        zeta = np.asarray(zeta, dtype=complex, order="F")
        nfreq = np.shape(zeta)[0]
        zflag = 1
        if any(
            abs(np.real(zeta)) > 1e-10
        ):  # if the provided frequency array is not Matsubara, set axis="r"
            axis = "r"
    else:
        zeta = np.asarray([0.0], dtype=complex, order="F")
        if axis == "m":
            nfreq = c_int.in_dll(self.library, "Lmats").value
        else:
            nfreq = c_int.in_dll(self.library, "Lreal").value
        zflag = 0

    # ishape
    if ishape is None:
        ishape = self.dim_hloc + 1

    # axis
    if axis == "m":
        axisint = 0
    elif axis == "r":
        axisint = 1
    else:
        raise ValueError("get_gimp: axis can only be 'm' or 'r'")

    # typ
    if typ == "n":
        typint = 0
    elif typ == "a":
        typint = 1
    else:
        raise ValueError("get_gimp: typ can only be 'n' or 'a'")

    if self.Nineq == 0:
        if ilat is not None:
            raise ValueError("ilat is not defined in single-impurity DMFT")
        if ishape == 3:
            gimp = np.zeros(
                [nspin_aux * norb_aux, nspin_aux * norb_aux, nfreq],
                dtype=complex,
                order="F",
            )
            ed_get_gimp_site_n3(gimp, axisint, typint, zeta, nfreq, zflag)
        elif ishape == 5:
            gimp = np.zeros(
                [nspin_aux, nspin_aux, norb_aux, norb_aux, nfreq],
                dtype=complex,
                order="F",
            )
            ed_get_gimp_site_n5(gimp, axisint, typint, zeta, nfreq, zflag)
        else:
            raise ValueError("Shape(array) != 3,5 in get_gimp_site")
        return gimp
    else:
        if ishape == 3:
            gimp = np.zeros(
                [
                    self.Nineq * nspin_aux * norb_aux,
                    self.Nineq * nspin_aux * norb_aux,
                    nfreq,
                ],
                dtype=complex,
                order="F",
            )
            ed_get_gimp_site_n3(gimp, self.Nineq, axisint, typint, zeta, nfreq, zflag)
        elif ishape == 4:
            gimp = np.zeros(
                [self.Nineq, nspin_aux * norb_aux, nspin_aux * norb_aux, nfreq],
                dtype=complex,
                order="F",
            )
            ed_get_gimp_site_n4(gimp, self.Nineq, axisint, typint, zeta, nfreq, zflag)
        elif ishape == 6:
            gimp = np.zeros(
                [self.Nineq, nspin_aux, nspin_aux, norb_aux, norb_aux, nfreq],
                dtype=complex,
                order="F",
            )
            ed_get_gimp_site_n6(gimp, self.Nineq, axisint, typint, zeta, nfreq, zflag)
        else:
            raise ValueError("Shape(array) != 3,4,6 in get_gimp_lattice")
        if ilat is not None and ishape != 3:
            return gimp[ilat]
        else:
            return gimp


# Anderson Impurity Model functions


# gimp
def get_g0and(self, zeta, bath, ishape=None, typ="n"):
    """

       This function calculates the value of the Anderson Impurity Model's \
       noninteracting Green's function on a given frequency array.

       :type zeta: complex 
       :param zeta: the array of frequencies (only frequencies on the real and \
        imaginary axes are supported)
       
       :type bath: float
       :param bath: the user-accessible bath array  
            
       :type ishape: int 
       :param ishape: this variable determines the shape of the returned array.\
        Possible values:
       
        * :code:`None`: the same shape as :code:`Hloc` plus one axis for frequency 
        * :code:`3`: the output array will have shape [ :data:`Nspin` \
          :math:`\\cdot` :data:`Norb` ,  :data:`Nspin` :math:`\\cdot`  \
          :data:`Norb` , :code:`len(zeta)` ]
        * :code:`5`: the output array will have shape [ :data:`Nspin` ,  \
          :data:`Nspin` ,  :data:`Norb` ,  :data:`Norb` ,  :code:`len(zeta)` ]
          
       :type typ: str 
       :param typ: whether to return the normal or anomalous Green's function \
        (for the superconducting case). Can be :code:`n` for normal or :code:`a`\
        for anomalous.
       
       :raise ValueError: If :code:`zeta` is not completely real or completely \
        imaginary
       :raise ValueError: If :code:`ishape` is not 3 or 5.
               
       :return: An array of complex that contains :math:`G^{And}_{0}(z)` function \
         along the specific frequency array, with dimension set by \
        :code:`ishape` and :code:`zeta`.  
       :rtype: np.array(dtype=complex)
       
    """

    ed_get_g0and_n3 = self.library.get_g0and_n3
    ed_get_g0and_n3.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=1, flags="F_CONTIGUOUS"),
        c_int,
        np.ctypeslib.ndpointer(dtype=float, ndim=1, flags="F_CONTIGUOUS"),
        c_int,
        np.ctypeslib.ndpointer(dtype=complex, ndim=3, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        c_char_p,
        c_char_p,
    ]
    ed_get_g0and_n3.restype = None

    ed_get_g0and_n5 = self.library.get_g0and_n5
    ed_get_g0and_n5.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=1, flags="F_CONTIGUOUS"),
        c_int,
        np.ctypeslib.ndpointer(dtype=float, ndim=1, flags="F_CONTIGUOUS"),
        c_int,
        np.ctypeslib.ndpointer(dtype=complex, ndim=5, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        c_char_p,
        c_char_p,
    ]
    ed_get_g0and_n5.restype = None

    norb_aux = c_int.in_dll(self.library, "Norb").value
    nspin_aux = c_int.in_dll(self.library, "Nspin").value

    zeta = zeta.astype(complex)

    nfreq = np.shape(zeta)[0]
    dimbath = np.shape(bath)[0]

    if any(abs(np.real(zeta)) > 1e-10):
        axis = "r"
    elif any(abs(np.imag(zeta)) > 1e-10):
        axis = "m"
    else:
        raise ValueError(
            "get_g0and: frequencies can only be purely real or purely imaginary"
        )
    if ishape is None:
        ishape = self.dim_hloc + 1

    if ishape == 3:
        G0and = np.zeros(
            [nspin_aux * norb_aux, nspin_aux * norb_aux, nfreq],
            dtype=complex,
            order="F",
        )
        DimG0and = np.asarray(
            [nspin_aux * norb_aux, nspin_aux * norb_aux, nfreq],
            dtype=np.int64,
            order="F",
        )
        ed_get_g0and_n3(
            zeta,
            nfreq,
            bath,
            dimbath,
            G0and,
            DimG0and,
            c_char_p(axis.encode()),
            c_char_p(typ.encode()),
        )
    elif ishape == 5:
        G0and = np.zeros(
            [nspin_aux, nspin_aux, norb_aux, norb_aux, nfreq], dtype=complex, order="F"
        )
        DimG0and = np.asarray(
            [nspin_aux, nspin_aux, norb_aux, norb_aux, nfreq], dtype=np.int64, order="F"
        )
        ed_get_g0and_n5(
            zeta,
            nfreq,
            bath,
            dimbath,
            G0and,
            DimG0and,
            c_char_p(axis.encode()),
            c_char_p(typ.encode()),
        )
    else:
        raise ValueError("Shape(array) != 3,5 in get_g0and")
    return G0and


# Delta
def get_delta(self, zeta, bath, ishape=None, typ="n"):
    """

       This function calculates the value of the Anderson Impurity Model's \
       hybridization function on a given frequency array.

       :type zeta: complex 
       :param zeta: the array of frequencies (only frequencies on the real and \
        imaginary axes are supported)   

       :type bath: float
       :param bath: the user-accessible bath array  
            
       :type ishape: int 
       :param ishape: this variable determines the shape of the returned array. \
        Possible values:
       
        * :code:`None`: the same shape as :code:`Hloc` plus one axis for frequency 
        * :code:`3`: the output array will have shape [ :data:`Nspin` \
          :math:`\\cdot` :data:`Norb` ,  :data:`Nspin` :math:`\\cdot` \
          :data:`Norb` , :code:`len(zeta)` ]
        * :code:`5`: the output array will have shape [ :data:`Nspin` ,  \
          :data:`Nspin` ,  :data:`Norb` ,  :data:`Norb` ,  :code:`len(zeta)` ]
          
       :type typ: str 
       :param typ: whether to return the normal or anomalous Green's function \
        (for the superconducting case). Can be :code:`n` for normal or :code:`a`\
        for anomalous.
       
       :raise ValueError: If :code:`zeta` is not completely real or completely \
        imaginary
       :raise ValueError: If :code:`ishape` is not 3 or 5.
               
       :return: An array of complex that contains :math:`\\Delta(z)` along the \
        specific frequency array, with dimension set by :code:`ishape` and \
        :code:`zeta`.  
       :rtype: np.array(dtype=complex) 
       
    """

    ed_get_delta_n3 = self.library.get_delta_n3
    ed_get_delta_n3.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=1, flags="F_CONTIGUOUS"),
        c_int,
        np.ctypeslib.ndpointer(dtype=float, ndim=1, flags="F_CONTIGUOUS"),
        c_int,
        np.ctypeslib.ndpointer(dtype=complex, ndim=3, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        c_char_p,
        c_char_p,
    ]
    ed_get_delta_n3.restype = None

    ed_get_delta_n5 = self.library.get_delta_n5
    ed_get_delta_n5.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=1, flags="F_CONTIGUOUS"),
        c_int,
        np.ctypeslib.ndpointer(dtype=float, ndim=1, flags="F_CONTIGUOUS"),
        c_int,
        np.ctypeslib.ndpointer(dtype=complex, ndim=5, flags="F_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags="F_CONTIGUOUS"),
        c_char_p,
        c_char_p,
    ]
    ed_get_delta_n5.restype = None

    norb_aux = c_int.in_dll(self.library, "Norb").value
    nspin_aux = c_int.in_dll(self.library, "Nspin").value

    zeta = zeta.astype(complex)

    nfreq = np.shape(zeta)[0]
    dimbath = np.shape(bath)[0]

    if any(abs(np.real(zeta)) > 1e-10):
        axis = "r"
    elif any(abs(np.imag(zeta)) > 1e-10):
        axis = "m"
    else:
        raise ValueError(
            "get_delta: frequencies can only be purely real or purely imaginary"
        )
    if ishape is None:
        ishape = self.dim_hloc + 1

    if ishape == 3:
        Delta = np.zeros(
            [nspin_aux * norb_aux, nspin_aux * norb_aux, nfreq],
            dtype=complex,
            order="F",
        )
        DimDelta = np.asarray(
            [nspin_aux * norb_aux, nspin_aux * norb_aux, nfreq],
            dtype=np.int64,
            order="F",
        )
        ed_get_delta_n3(
            zeta,
            nfreq,
            bath,
            dimbath,
            Delta,
            DimDelta,
            c_char_p(axis.encode()),
            c_char_p(typ.encode()),
        )
    elif ishape == 5:
        Delta = np.zeros(
            [nspin_aux, nspin_aux, norb_aux, norb_aux, nfreq], dtype=complex, order="F"
        )
        DimDelta = np.asarray(
            [nspin_aux, nspin_aux, norb_aux, norb_aux, nfreq], dtype=np.int64, order="F"
        )
        ed_get_delta_n5(
            zeta,
            nfreq,
            bath,
            dimbath,
            Delta,
            DimDelta,
            c_char_p(axis.encode()),
            c_char_p(typ.encode()),
        )
    else:
        raise ValueError("Shape(array) != 3,5 in get_delta")
    return Delta


###################
# Susceptibilities #
###################


def get_chi(self, chan="spin", zeta=None, axis=None, ilat=None):
    """

       This function calculates the value of the Anderson Impurity Model's \
       response functions :math:`\\chi` in the spin, density, pairing and \
       exciton channels.

       :type chan: str 
       :param chan: the channel of the response function. Possible values are \
        :code:`spin` ( :code:`s` ), :code:`dens` ( :code:`c` ), :code:`pair` \
        ( :code:`p` ), :code:`exct` ( :code:`e` ). Default is :code:`spin`.
       
       :type zeta: complex 
       :param zeta: the array of frequencies or times (only frequencies on \
        the real imaginary axes and imaginary times are supported). If no array
        is provided, depending on the value of :data:`axis`, an array will be \
        generated on the positive Matsubara axis ( :data:`Lmats` values ), on \
        the real axis ( :data:`Lreal` values between :code:`0` and  \
        :data:`wfin`) or on the imaginary time axis (:data:`Ltau` values \
        between 0 and :data:`beta`)

       :type axis: str
       :param axis: the axis on which to calculate :math:`\\chi`. Possible values \
       :code:`matsubara` ( :code:`m`), :code:`real` ( :code:`r`), :code:`time` \
        ( :code:`t`). This parameter has to be specificed.
       
       :type ilat: int
       :param ilat: for real-space DMFT, if this flag is specified only the \
        susceptibility for the relative inequivalent atom is returned
        specified.
       
       
       :raise ValueError: If :f:var:`ed_mode` is not :code:`normal`
       :raise ValueError: If :code:`axis` not provided or invalid
       :raise ValueError: If :code:`chan` is invalid
               
       :return: An array of complex that contains :math:`\\chi` along the chosen \
        axis. The dimension of the array depends on the chosen channel:

        * [ :data:`Norb`, :data:`Norb` , :code:`len(zeta)` ] for \
          channel :code:`spin`, :code:`dens`, :code:`pair`
        * [ :code:`3`, :data:`Norb`, :data:`Norb` , :code:`len(zeta)` ] for channel \
          :code:`exct`, corresponding to singlet, triplet(xy) and triplet(z)
          
        One dimension corresponding to the number of inequivalent sites is added \
        at the beginning for the case of real-space DMFT and if ilat is None.
       :rtype: np.array(dtype=complex) 
               
       
    """
    ed_get_spinchi = self.library.ed_get_spinchi
    ed_get_spinchi.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=4, flags="F_CONTIGUOUS"),  # self
        np.ctypeslib.ndpointer(dtype=complex, ndim=1, flags="F_CONTIGUOUS"),  # zeta
        c_int,  # dim_zeta
        c_int,  # zetaflag
        c_int,  # axis
        c_int,  # nsites
        c_int,  # latticeflag
    ]
    ed_get_spinchi.restype = None

    ed_get_denschi = self.library.ed_get_denschi
    ed_get_denschi.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=4, flags="F_CONTIGUOUS"),  # self
        np.ctypeslib.ndpointer(dtype=complex, ndim=1, flags="F_CONTIGUOUS"),  # zeta
        c_int,  # dim_zeta
        c_int,  # zetaflag
        c_int,  # axis
        c_int,  # nsites
        c_int,  # latticeflag
    ]
    ed_get_denschi.restype = None

    ed_get_pairchi = self.library.ed_get_pairchi
    ed_get_pairchi.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=4, flags="F_CONTIGUOUS"),  # self
        np.ctypeslib.ndpointer(dtype=complex, ndim=1, flags="F_CONTIGUOUS"),  # zeta
        c_int,  # dim_zeta
        c_int,  # zetaflag
        c_int,  # axis
        c_int,  # nsites
        c_int,  # latticeflag
    ]
    ed_get_pairchi.restype = None

    ed_get_exctchi = self.library.ed_get_exctchi
    ed_get_exctchi.argtypes = [
        np.ctypeslib.ndpointer(dtype=complex, ndim=5, flags="F_CONTIGUOUS"),  # self
        np.ctypeslib.ndpointer(dtype=complex, ndim=1, flags="F_CONTIGUOUS"),  # zeta
        c_int,  # dim_zeta
        c_int,  # zetaflag
        c_int,  # axis
        c_int,  # nsites
        c_int,  # latticeflag
    ]
    ed_get_exctchi.restype = None

    if self.Nineq == 0:
        Nsites = 1
        latticeflag = 0
    else:
        Nsites = self.Nineq
        latticeflag = 1

    aux_norb = c_int.in_dll(self.library, "Norb").value
    aux_Lmats = c_int.in_dll(self.library, "Lmats").value
    aux_Lreal = c_int.in_dll(self.library, "Lreal").value
    aux_Ltau = c_int.in_dll(self.library, "Ltau").value
    edmode = self.get_ed_mode()

    if edmode != 1:
        raise ValueError(
            "Susceptibility calculation not supported for ed_mode " "not = normal"
        )

    zetaflag = 1

    if axis is None:
        raise ValueError("Axis is required")
    else:
        if zeta is None:
            if axis == "m":
                zeta = np.array([0.0], dtype=complex)
                zetaflag = 0
                axisflag = 0
                nfreq = aux_Lmats
            if axis == "r":
                zeta = np.array([0.0], dtype=complex)
                zetaflag = 0
                axisflag = 1
                nfreq = aux_Lreal
            if axis == "t":
                zeta = np.array([0.0], dtype=complex)
                zetaflag = 0
                axisflag = 2
                nfreq = aux_Ltau
        else:
            if axis == "m":
                axisflag = 0
            elif axis == "r":
                axisflag = 1
            elif axis == "t":
                axisflag = 2
            else:
                raise ValueError("axis can only be m, r or t")
            nfreq = np.shape(zeta)[0]

    if chan.lower() == "spin" or chan.lower() == "s":
        chi = np.zeros([Nsites, aux_norb, aux_norb, nfreq], dtype=complex, order="F")
        ed_get_spinchi(chi, zeta, nfreq, zetaflag, axisflag, Nsites, latticeflag)
    if chan.lower() == "dens" or chan.lower() == "d":
        chi = np.zeros([Nsites, aux_norb, aux_norb, nfreq], dtype=complex, order="F")
        ed_get_denschi(chi, zeta, nfreq, zetaflag, axisflag, Nsites, latticeflag)
    if chan.lower() == "pair" or chan.lower() == "p":
        chi = np.zeros([Nsites, aux_norb, aux_norb, nfreq], dtype=complex, order="F")
        ed_get_pairchi(chi, zeta, nfreq, zetaflag, axisflag, Nsites, latticeflag)
    if chan.lower() == "exct" or chan.lower() == "e":
        chi = np.zeros([Nsites, 3, aux_norb, aux_norb, nfreq], dtype=complex, order="F")
        ed_get_exctchi(chi, zeta, nfreq, zetaflag, axisflag, Nsites, latticeflag)

    if self.Nineq == 0:
        return chi[0]
    else:
        if ilat == None:
            return chi
        else:
            return chi[ilat]
