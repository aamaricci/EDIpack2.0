from ctypes import *
import numpy as np
import os,sys
import types

#get_bath_dimension
def get_bath_dimension(self):
    get_bath_dimension_wrap = self.library.get_bath_dimension
    get_bath_dimension_wrap.argtypes = None  
    get_bath_dimension_wrap.restype = c_int
    return get_bath_dimension_wrap()

#init_hreplica
   
def set_Hreplica(self,hvec,lambdavec):
    init_hreplica_symmetries_d5 = self.library.init_Hreplica_symmetries_d5
    init_hreplica_symmetries_d5.argtypes =[np.ctypeslib.ndpointer(dtype=complex,ndim=5, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS')]
    init_hreplica_symmetries_d5.restype = None
    
    init_hreplica_symmetries_d3 = self.library.init_Hreplica_symmetries_d3
    init_hreplica_symmetries_d3.argtypes =[np.ctypeslib.ndpointer(dtype=complex,ndim=3, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS')]
    init_hreplica_symmetries_d3.restype = None
    
    init_hreplica_symmetries_lattice_d5 = self.library.init_Hreplica_symmetries_lattice_d5
    init_hreplica_symmetries_lattice_d5.argtypes =[np.ctypeslib.ndpointer(dtype=complex,ndim=5, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=float,ndim=3, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS')]
    init_hreplica_symmetries_lattice_d5.restype = None
    
    init_hreplica_symmetries_lattice_d3 = self.library.init_Hreplica_symmetries_lattice_d3
    init_hreplica_symmetries_lattice_d3.argtypes =[np.ctypeslib.ndpointer(dtype=complex,ndim=3, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=float,ndim=3, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS')]
    init_hreplica_symmetries_lattice_d3.restype = None

    
    aux_norb=c_int.in_dll(self.library, "Norb").value
    aux_nspin=c_int.in_dll(self.library, "Nspin").value
    dim_hvec = np.asarray(np.shape(hvec),dtype=np.int64,order="F")
    dim_lambdavec = np.asarray(np.shape(lambdavec),dtype=np.int64,order="F")


    if(len(dim_hvec) == 3):
        if(len(dim_lambdavec)==2):
            init_hreplica_symmetries_d3(hvec,dim_hvec,lambdavec,dim_lambdavec)
        elif(len(Ddim_lambdavec)==3):
            init_hreplica_symmetries_lattice_d3(hvec,dim_hvec,lambdavec,dim_lambdavec)
        else:
             raise ValueError('Shape(lambdavec) != 2 or 3  in set_Hreplica')
    elif(len(dim_hvec) == 5):
        if(len(dim_lambdavec)==2):
            init_hreplica_symmetries_d5(hvec,dim_hvec,lambdavec,dim_lambdavec)
        elif(len(Ddim_lambdavec)==3):
            init_hreplica_symmetries_lattice_d5(hvec,dim_hvec,lambdavec,dim_lambdavec)
        else:
             raise ValueError('Shape(lambdavec) != 2 or 3  in set_Hreplica')
    else:
         raise ValueError('Shape(Hvec) != 3 or 5  in set_Hreplica')
    return ;
    
#break_symmetry_bath
def break_symmetry_bath(self, bath, field, sign, save=True):
    break_symmetry_bath_site = self.library.break_symmetry_bath_site
    break_symmetry_bath_site.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),
                                         np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                         c_double,
                                         c_double,
                                         c_int]
    break_symmetry_bath_site.restype = None

    break_symmetry_bath_ineq = self.library.break_symmetry_bath_ineq
    break_symmetry_bath_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),
                                         np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                         c_double,
                                         c_double,
                                         c_int]
    break_symmetry_bath_ineq.restype = None

    if save:
        save_int=1
    else:
        save_int=0
    bath_shape = np.asarray(np.shape(bath),dtype=np.int64,order="F")
    if (len(bath_shape)) == 1:
        break_symmetry_bath_site(bath,bath_shape,field,float(sign),save_int)
    else:
        break_symmetry_bath_ineq(bath,bath_shape,field,float(sign),save_int)
    return bath
    
#spin_symmetrize_bath
    
def spin_symmetrize_bath(self, bath, save=True):
    spin_symmetrize_bath_site = self.library.spin_symmetrize_bath_site
    spin_symmetrize_bath_site.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),
                                          np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                          c_int]
    spin_symmetrize_bath_site.restypes = None

    spin_symmetrize_bath_ineq = self.library.spin_symmetrize_bath_ineq
    spin_symmetrize_bath_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),
                                          np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'), 
                                          c_int]
    spin_symmetrize_bath_ineq.restypes = None
    if save:
        save_int=1
    else:
        save_int=0
    bath_shape = np.asarray(np.shape(bath),dtype=np.int64,order="F")
    if (len(bath_shape)) == 1:
        spin_symmetrize_bath_site(bath,bath_shape,save_int)
    else:
        spin_symmetrize_bath_ineq(bath,bath_shape,save_int)
    return bath
    
#orb_symmetrize_bath
def orb_symmetrize_bath(self, bath, save=True):
    orb_symmetrize_bath_site = self.library.orb_symmetrize_bath_site
    orb_symmetrize_bath_site.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),
                                          np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                          c_int]
    orb_symmetrize_bath_site.restypes = None

    orb_symmetrize_bath_ineq = self.library.orb_symmetrize_bath_ineq
    orb_symmetrize_bath_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),
                                         np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'), 
                                         c_int]
    orb_symmetrize_bath_ineq.restypes = None
    
    if save:
        save_int=1
    else:
        save_int=0
    bath_shape = np.asarray(np.shape(bath),dtype=np.int64,order="F")
    if (len(bath_shape)) == 1:
        orb_symmetrize_bath_site(bath,bath_shape,save_int)
    else:
        orb_symmetrize_bath_ineq(bath,bath_shape,save_int)
    return bath
    
#orb_equality_bath
    
def orb_equality_bath(self, bath, indx, save=True):
    orb_equality_bath_site = self.library.orb_equality_bath_site
    orb_equality_bath_site.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),
                                       np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'), 
                                       c_int,
                                       c_int]
    orb_equality_bath_site.restypes = None

    orb_equality_bath_ineq = self.library.orb_equality_bath_ineq
    orb_equality_bath_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),
                                       np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                       c_int,
                                       c_int]
    orb_equality_bath_ineq.restypes = None
    aux_norb=c_int.in_dll(self.library, "Norb").value
    if save:
        save_int=1
    else:
        save_int=0
    bath_shape = np.asarray(np.shape(bath),dtype=np.int64,order="F")
    if (indx < 0) or (indx >= aux_norb):
        raise ValueError("orb_equality_bath: orbital index should be in [0,Norb]")
    else:
        indx = indx + 1 #python to fortran convention 
        if (len(bath_shape)) == 1:
            orb_equality_bath_site(bath,bath_shape,indx,save_int)
        else:
            orb_equality_bath_ineq(bath,bath_shape,indx,save_int)
    return bath
    
    
#ph_symmetrize_bath
def ph_symmetrize_bath(self, bath, save):
    ph_symmetrize_bath_site = self.library.ph_symmetrize_bath_site
    ph_symmetrize_bath_site.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),
                                        np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                        c_int]
    ph_symmetrize_bath_site.restypes = None

    ph_symmetrize_bath_ineq = self.library.ph_symmetrize_bath_ineq
    ph_symmetrize_bath_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),
                                        np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                        c_int]
    ph_symmetrize_bath_ineq.restypes = None
    if save:
        save_int=1
    else:
        save_int=0
    bath_shape = np.asarray(np.shape(bath),dtype=np.int64,order="F")
    if (len(bath_shape)) == 1:
        ph_symmetrize_bath_site(bath,bath_shape,save_int)
    else:
        ph_symmetrize_bath_ineq(bath,bath_shape,save_int)
    return bath

