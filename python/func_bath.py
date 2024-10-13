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
def set_hreplica(self,hvec,lambdavec):
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
        elif(len(dim_lambdavec)==3):
            init_hreplica_symmetries_lattice_d5(hvec,dim_hvec,lambdavec,dim_lambdavec)
        else:
             raise ValueError('Shape(lambdavec) != 2 or 3  in set_Hreplica')
    else:
         raise ValueError('Shape(Hvec) != 3 or 5  in set_Hreplica')
    return ;
    
#init_hgeneral
def set_hgeneral(self,hvec,lambdavec):
    init_hgeneral_symmetries_d5 = self.library.init_Hgeneral_symmetries_d5
    init_hgeneral_symmetries_d5.argtypes =[np.ctypeslib.ndpointer(dtype=complex,ndim=5, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS')]
    init_hgeneral_symmetries_d5.restype = None
    
    init_hgeneral_symmetries_d3 = self.library.init_Hgeneral_symmetries_d3
    init_hgeneral_symmetries_d3.argtypes =[np.ctypeslib.ndpointer(dtype=complex,ndim=3, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS')]
    init_hgeneral_symmetries_d3.restype = None
    
    init_hgeneral_symmetries_lattice_d5 = self.library.init_Hgeneral_symmetries_lattice_d5
    init_hgeneral_symmetries_lattice_d5.argtypes =[np.ctypeslib.ndpointer(dtype=complex,ndim=5, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=float,ndim=3, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS')]
    init_hgeneral_symmetries_lattice_d5.restype = None
    
    init_hgeneral_symmetries_lattice_d3 = self.library.init_Hgeneral_symmetries_lattice_d3
    init_hgeneral_symmetries_lattice_d3.argtypes =[np.ctypeslib.ndpointer(dtype=complex,ndim=3, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=float,ndim=3, flags='F_CONTIGUOUS'),
                                           np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS')]
    init_hgeneral_symmetries_lattice_d3.restype = None

    
    aux_norb=c_int.in_dll(self.library, "Norb").value
    aux_nspin=c_int.in_dll(self.library, "Nspin").value
    dim_hvec = np.asarray(np.shape(hvec),dtype=np.int64,order="F")
    dim_lambdavec = np.asarray(np.shape(lambdavec),dtype=np.int64,order="F")


    if(len(dim_hvec) == 3):
        if(len(dim_lambdavec)==2):
            init_hgeneral_symmetries_d3(hvec,dim_hvec,lambdavec,dim_lambdavec)
        elif(len(Ddim_lambdavec)==3):
            init_hgeneral_symmetries_lattice_d3(hvec,dim_hvec,lambdavec,dim_lambdavec)
        else:
             raise ValueError('Shape(lambdavec) != 2 or 3  in set_Hgeneral')
    elif(len(dim_hvec) == 5):
        if(len(dim_lambdavec)==2):
            init_hgeneral_symmetries_d5(hvec,dim_hvec,lambdavec,dim_lambdavec)
        elif(len(dim_lambdavec)==3):
            init_hgeneral_symmetries_lattice_d5(hvec,dim_hvec,lambdavec,dim_lambdavec)
        else:
             raise ValueError('Shape(lambdavec) != 2 or 3  in set_Hgeneral')
    else:
         raise ValueError('Shape(Hvec) != 3 or 5  in set_Hgeneral')
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
def orb_symmetrize_bath(self, bath,orb1,orb2, save=True):
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
        orb_symmetrize_bath_site(bath,bath_shape,orb1+1,orb2+1,save_int)
    else:
        orb_symmetrize_bath_ineq(bath,bath_shape,orb1+1,orb2+1,save_int)
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
    
#save array as .restart file
def save_array_as_bath(self, bath):
    save_array_as_bath_site = self.library.save_array_as_bath_site
    save_array_as_bath_site.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),
                                         np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS')]
    save_array_as_bath_site.restypes = None

    save_array_as_bath_ineq = self.library.save_array_as_bath_ineq
    save_array_as_bath_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS'),
                                        np.ctypeslib.ndpointer(dtype=np.int64,ndim=1, flags='F_CONTIGUOUS')]
    save_array_as_bath_ineq.restypes = None
    
    bath_shape = np.asarray(np.shape(bath),dtype=np.int64,order="F")
    if (len(bath_shape)) == 1:
        save_array_as_bath_site(bath,bath_shape)
    else:
        save_array_as_bath_ineq(bath,bath_shape)
    return ;


#auxiliary functions to get/set bath structure. Only works for single-site. User has to do a loop on sites


def bath_inspect(self,bath=None,e=None,v=None,d=None,u=None):
    
    aux_norb=c_int.in_dll(self.library, "Norb").value
    aux_nspin=c_int.in_dll(self.library, "Nspin").value
    aux_nbath=c_int.in_dll(self.library, "Nbath").value
    
    settings=(self.get_ed_mode(),self.get_bath_type())
    if settings == (1,1): #normal ed mode, normal bath
        if bath is None and e is not None and v is not None:
            e=np.asarray(e,order="F")
            v=np.asarray(v,order="F")
            try:
                if np.shape(e) != (aux_nspin,aux_norb,aux_nbath):
                    raise ValueError("e must be (nspin,norb,nbath)")
                if np.shape(v) != (aux_nspin,aux_norb,aux_nbath):
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
                        io = stride + ibath + (iorb)*aux_nbath + (ispin)*aux_nbath*aux_norb
                        bath[io]=e[ispin,iorb,ibath]
            stride = aux_nspin * aux_norb * aux_nbath
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = stride + ibath + (iorb)*aux_nbath + (ispin)*aux_nbath*aux_norb
                        bath[io]=v[ispin,iorb,ibath]
            return bath
                
        elif bath is not None and e is None and v is None: #e and v are none
            bath = np.asarray(bath,order="F")
            Nb = self.get_bath_dimension()
            if np.shape(bath)[0] != Nb:
                raise ValueError("bath has the wrong length")
                
            e=np.zeros((aux_nspin,aux_norb,aux_nbath))
            v=np.zeros((aux_nspin,aux_norb,aux_nbath))
            
            stride = 0
            io = 0
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = stride + ibath + (iorb)*aux_nbath + (ispin)*aux_nbath*aux_norb
                        e[ispin,iorb,ibath] = bath[io]
            stride = aux_nspin * aux_norb * aux_nbath            
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = stride + ibath + (iorb)*aux_nbath + (ispin)*aux_nbath*aux_norb
                        v[ispin,iorb,ibath] = bath[io]
            return e,v
        else:
            raise ValueError("Wrong input for normal/normal")
            
    elif settings == (2,1): #superc ed mode, normal bath
        if bath is None and e is not None and v is not None and d is not None:
            e=np.asarray(e,order="F")
            v=np.asarray(v,order="F")
            d=np.asarray(u,order="F")
            try:
                if np.shape(e) != (aux_nspin,aux_norb,aux_nbath):
                    raise ValueError("e must be (nspin,norb,nbath)")
                if np.shape(d) != (aux_nspin,aux_norb,aux_nbath):
                    raise ValueError("d must be (nspin,norb,nbath)")
                if np.shape(v) != (aux_nspin,aux_norb,aux_nbath):
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
                        io = stride + ibath + (iorb)*aux_nbath + (ispin)*aux_nbath*aux_norb
                        bath[io] = e[ispin,iorb,ibath]
            stride = aux_nspin * aux_norb * aux_nbath 
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = stride + ibath + (iorb)*aux_nbath + (ispin)*aux_nbath*aux_norb
                        bath[io] = d[ispin,iorb,ibath]
            stride = 2 * aux_nspin * aux_norb * aux_nbath 
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = stride + ibath + (iorb)*aux_nbath + (ispin)*aux_nbath*aux_norb
                        bath[io] = v[ispin,iorb,ibath]
            return bath
        elif bath is None and e is None and v is None and d is None:
            bath = np.asarray(bath,order="F")
            Nb = self.get_bath_dimension()
            if np.shape(bath)[0] != Nb:
                raise ValueError("bath has the wrong length")
                
            e=np.zeros((aux_nspin,aux_norb,aux_nbath))
            v=np.zeros((aux_nspin,aux_norb,aux_nbath))
            d=np.zeros((aux_nspin,aux_norb,aux_nbath))
            
            stride = 0
            io = 0
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = stride + ibath + (iorb)*aux_nbath + (ispin)*aux_nbath*aux_norb
                        e[ispin,iorb,ibath] = bath[io]
            stride = aux_nspin * aux_norb * aux_nbath 
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = stride + ibath + (iorb)*aux_nbath + (ispin)*aux_nbath*aux_norb
                        d[ispin,iorb,ibath] = bath[io]
            stride = 2 * aux_nspin * aux_norb * aux_nbath 
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = stride + ibath + (iorb)*aux_nbath + (ispin)*aux_nbath*aux_norb
                        v[ispin,iorb,ibath] = bath[io]
            return e,d,v
        else:
            raise ValueError("Wrong input for superc/normal")
            
    elif settings == (3,1): #nonsu2 ed mode, normal bath
        if bath is None and e is not None and v is not None and u is not None:
            try:
                e=np.asarray(e,order="F")
                v=np.asarray(v,order="F")
                u=np.asarray(u,order="F")
                if np.shape(e) != (aux_nspin,aux_norb,aux_nbath):
                    raise ValueError("e must be (nspin,norb,nbath)")
                if np.shape(v) != (aux_nspin,aux_norb,aux_nbath):
                    raise ValueError("v must be (nspin,norb,nbath)")
                if np.shape(u) != (aux_nspin,aux_norb,aux_nbath):
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
                        io = stride + ibath + (iorb)*aux_nbath + (ispin)*aux_nbath*aux_norb
                        bath[io] = e[ispin,iorb,ibath]
            stride = aux_nspin * aux_norb * aux_nbath 
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = stride + ibath + (iorb)*aux_nbath + (ispin)*aux_nbath*aux_norb
                        bath[io] = v[ispin,iorb,ibath]
            stride = 2 * aux_nspin * aux_norb * aux_nbath 
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = stride + ibath + (iorb)*aux_nbath + (ispin)*aux_nbath*aux_norb
                        bath[io] = u[ispin,iorb,ibath]
            return bath
        elif bath is not None and e is None and v is None and u is None:
            bath = np.asarray(bath,order="F")
            Nb = self.get_bath_dimension()
            if np.shape(bath)[0] != Nb:
                raise ValueError("bath has the wrong length")
                
            e=np.zeros((aux_nspin,aux_norb,aux_nbath))
            v=np.zeros((aux_nspin,aux_norb,aux_nbath))
            u=np.zeros((aux_nspin,aux_norb,aux_nbath))
            
            stride = 0
            io = 0
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = stride + ibath + (iorb)*aux_nbath + (ispin)*aux_nbath*aux_norb
                        e[ispin,iorb,ibath] = bath[io]
            stride = aux_nspin * aux_norb * aux_nbath 
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = stride + ibath + (iorb)*aux_nbath + (ispin)*aux_nbath*aux_norb
                        v[ispin,iorb,ibath] = bath[io]
            stride = 2 * aux_nspin * aux_norb * aux_nbath 
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = stride + ibath + (iorb)*aux_nbath + (ispin)*aux_nbath*aux_norb
                        u[ispin,iorb,ibath] = bath[io]
            return e,v,u
        else:
            raise ValueError("Wrong input for nonsu2/normal")
            
            
    elif settings == (1,2): #normal ed mode, hybrid bath
        if bath is None and e is not None and v is not None:
            try:
                e=np.asarray(e,order="F")
                v=np.asarray(v,order="F")
                if np.shape(e) != (aux_nspin,aux_nbath):
                    raise ValueError("e must be (nspin,nbath)")
                if np.shape(v) != (aux_nspin,aux_norb,aux_nbath):
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
                    io = stride + ibath + (ispin)*aux_nbath
                    bath[io] = e[ispin,ibath]
            stride = aux_nspin * aux_nbath
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = stride + ibath + (iorb)*aux_nbath + (ispin)*aux_nbath*aux_norb
                        bath[io] = v[ispin,iorb,ibath]
            return bath
                
        elif bath is not None and e is None and v is None:
            bath = np.asarray(bath,order="F")
            Nb = self.get_bath_dimension()
            if np.shape(bath)[0] != Nb:
                raise ValueError("bath has the wrong length")
                
            e = np.zeros((aux_nspin,aux_nbath))
            v = np.zeros((aux_nspin,aux_norb,aux_nbath))
            
            stride = 0
            io = 0
            for ispin in range(aux_nspin):
                for ibath in range(aux_nbath):
                    io = stride + ibath + ispin*aux_nbath
                    e[ispin,ibath] = bath[io]
            stride = aux_nspin * aux_nbath            
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = stride + ibath + (iorb)*aux_nbath + (ispin)*aux_nbath*aux_norb
                        v[ispin,iorb,ibath] = bath[io]
            return e,v
        else:
            raise ValueError("Wrong input for normal/hybrid")
            
    elif settings == (2,2): #superc ed mode, hybrid bath
        if bath is None and e is not None and v is not None and d is not None:
            try:
                e=np.asarray(e,order="F")
                d=np.asarray(d,order="F")
                v=np.asarray(v,order="F")
                if np.shape(e) != (aux_nspin,aux_nbath):
                    raise ValueError("e must be (nspin,nbath)")
                if np.shape(d) != (aux_nspin,aux_nbath):
                    raise ValueError("d must be (nspin,nbath)")
                if np.shape(v) != (aux_nspin,aux_norb,aux_nbath):
                    raise ValueError("v must be (nspin,norb,nbath)")
            except:
                raise ValueError("e,d or v have wrong dimension")

            Nb = self.get_bath_dimension()
            bath = np.zeros(Nb)
            stride = 0
            io = 0
            for ispin in range(aux_nspin):
                for ibath in range(aux_nbath):
                    io = stride + ibath + (ispin)*aux_nbath
                    bath[io] = e[ispin,ibath]
            stride = aux_nspin * aux_nbath 
            for ispin in range(aux_nspin):
                for ibath in range(aux_nbath):
                    io = stride + ibath + (ispin)*aux_nbath
                    bath[io] = bath_d[ispin,ibath]
            stride = 2 * aux_nspin * aux_nbath 
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = stride + ibath + (iorb)*aux_nbath + (ispin)*aux_nbath*aux_norb
                        bath[io] = v[ispin,iorb,ibath]
            return bath
        elif bath is not None and e is None and v is None and d is None:
            bath = np.asarray(bath,order="F")
            Nb = self.get_bath_dimension()
            if np.shape(bath)[0] != Nb:
                raise ValueError("bath has the wrong length")
                
            e = np.zeros((aux_nspin,aux_nbath))
            d = np.zeros((aux_nspin,aux_nbath))
            v = np.zeros((aux_nspin,aux_norb,aux_nbath))
            
            stride = 0
            io = 0
            for ispin in range(aux_nspin):
                for ibath in range(aux_nbath):
                    io = stride + ibath + (ispin)*aux_nbath
                    e[ispin,ibath] = bath[io]
            stride = aux_nspin * aux_nbath 
            for ispin in range(aux_nspin):
                for ibath in range(aux_nbath):
                    io = stride + ibath + (ispin)*aux_nbath
                    d[ispin,ibath] = bath[io]
            stride = 2 * aux_nspin * aux_nbath 
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = stride + ibath + (iorb)*aux_nbath + (ispin)*aux_nbath*aux_norb
                        v[ispin,iorb,ibath] = bath[io]
            return e,d,v
        else:
            raise ValueError("Wrong input for superc/hybrid")
            
    elif settings == (3,2): #nonsu2 ed mode, hybrid bath
        if bath is None and e is not None and v is not None and u is not None:
            try:
                e=np.asarray(e,order="F")
                v=np.asarray(v,order="F")
                u=np.asarray(u,order="F")
                if np.shape(e) != (aux_nspin,aux_nbath):
                    raise ValueError("e must be (nspin,norb,nbath)")
                if np.shape(v) != (aux_nspin,aux_norb,aux_nbath):
                    raise ValueError("v must be (nspin,norb,nbath)")
                if np.shape(u) != (aux_nspin,aux_norb,aux_nbath):
                    raise ValueError("u must be (nspin,norb,nbath)")
            except:
                raise ValueError("e,v or u have wrong dimension")

            Nb = self.get_bath_dimension()
            bath = np.zeros(Nb)
            
            stride = 0
            io = 0
            for ispin in range(aux_nspin):
                for ibath in range(aux_nbath):
                    io = stride + ibath + (ispin)*aux_nbath
                    bat[io] = e[ispin,ibath]
            stride = aux_nspin * aux_nbath 
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = stride + ibath + (iorb)*aux_nbath + (ispin)*aux_nbath*aux_norb
                        bath[io] = v[ispin,iorb,ibath]
            stride = aux_nspin * aux_nbath + aux_nspin * aux_norb * aux_nbath 
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = stride + ibath + (iorb)*aux_nbath + (ispin)*aux_nbath*aux_norb
                        bath[io] = u[ispin,iorb,ibath]
            return bath
        elif bath is not None and e is None and v is None and u is None:
            bath = np.asarray(bath,order="F")
            Nb = self.get_bath_dimension()
            if np.shape(bath)[0] != Nb:
                raise ValueError("bath has the wrong length")
                
            e = np.zeros((aux_nspin,aux_nbath))
            v = np.zeros((aux_nspin,aux_norb,aux_nbath))
            u = np.zeros((aux_nspin,aux_norb,aux_nbath))
            
            stride = 0
            io = 0
            for ispin in range(aux_nspin):
                for ibath in range(aux_nbath):
                    io = stride + ibath + (ispin)*aux_nbath
                    e[ispin,ibath] = bath[io]
            stride = aux_nspin * aux_nbath 
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = stride + ibath + (iorb)*aux_nbath + (ispin)*aux_nbath*aux_norb
                        v[ispin,iorb,ibath] = bath[io]
            stride = aux_nspin * aux_nbath + aux_nspin * aux_norb * aux_nbath 
            for ispin in range(aux_nspin):
                for iorb in range(aux_norb):
                    for ibath in range(aux_nbath):
                        io = stride + ibath + (iorb)*aux_nbath + (ispin)*aux_nbath*aux_norb
                        u[ispin,iorb,ibath] = bath[io]
            return e,v,u
        else:
            raise ValueError("Wrong input for nonsu2/hybrid")
    else:
        raise ValueError("EDmode/bath combination not valid or not implemented.")
        
    

