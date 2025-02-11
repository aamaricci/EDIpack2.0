from ctypes import *
import numpy as np
import os,sys
import types

#################################
#Classes
#################################

#dummy classes, to be filled
class LinkSolver:
    def __init__(self,library):
        self.library = library
        self.Nineq = None
        self.dim_hloc = 0
        self.Nsym = None
#utils: colors and bold text
        self.PURPLE = '\033[95m'
        self.CYAN = '\033[96m'
        self.DARKCYAN = '\033[36m'
        self.BLUE = '\033[94m'
        self.GREEN = '\033[92m'
        self.YELLOW = '\033[93m'
        self.RED = '\033[91m'
        self.BOLD = '\033[1m'
        self.UNDERLINE = '\033[4m'
        self.COLOREND = '\033[0m'
        
class LinkFit:
    def __init__(self,library):
        self.library = library
        self.library = library
        self.Nineq = None
        self.dim_hloc = 0
        self.Nsym = None
#utils: colors and bold text
        self.PURPLE = '\033[95m'
        self.CYAN = '\033[96m'
        self.DARKCYAN = '\033[36m'
        self.BLUE = '\033[94m'
        self.GREEN = '\033[92m'
        self.YELLOW = '\033[93m'
        self.RED = '\033[91m'
        self.BOLD = '\033[1m'
        self.UNDERLINE = '\033[4m'
        self.COLOREND = '\033[0m'

#################################
#AUXILIARY FUNCTIONS
#################################

#function that will add a variable to the dummy class, will be called in variable definition
def add_global_variable(obj, dynamic_name, target_object, target_attribute):
    @property
    def getter(self):
        try:
            attrib = getattr(target_object, target_attribute)
            try: #this is for strings
                attrib = attrib.decode()
            except:
                pass
        except: #this is for arrays
            if(len(target_object)>1):
                return [target_object[x] for x in range(len(target_object))]
        return attrib

    @getter.setter
    def setter(self, new_value):
        try: #this is for arrays
            if(len(target_object)>1):
                minlength=min(len(target_object),len(new_value))
                target_object[0:minlength]=new_value[0:minlength]
        except:
            try:
                new_value = new_value.encode()
            except:
                pass
            setattr(target_object, target_attribute, new_value)

    # Dynamically add the property to the class
    setattr(obj.__class__, dynamic_name, getter)
    setattr(obj.__class__, dynamic_name, setter)


#get bath type
def get_bath_type(self):
    """
    
     This function returns an integer number related to the value of  :f:var:`bath_type` in the input file
    
      - :code:`1` for **normal** bath
      - :code:`2` for **hybrid** bath
      - :code:`3` for **replica** bath
      - :code:`4` for **general** bath
   
    :return: the integer index
    :rtype: int
    
   """

    get_bath_type_wrap = self.library.get_bath_type
    get_bath_type_wrap.argtypes = None  
    get_bath_type_wrap.restype = c_int
    return get_bath_type_wrap()
    
#get ed mode
def get_ed_mode(self):
    """
    
     This function returns an integer number related to the value of  :f:var:`ed_mode` in the input file
     
      - :code:`1` for **normal** mode
      - :code:`2` for **superc** mode
      - :code:`3` for **nonsu2** mode
   
     :return: the integer index
     :rtype: int
    
    """
    
    get_ed_mode_wrap = self.library.get_ed_mode
    get_ed_mode_wrap.argtypes = None  
    get_ed_mode_wrap.restype = c_int
    return get_ed_mode_wrap()

######################################
# Load shared library with C-bindings
######################################

system = sys.platform
libext = '.so' 
if(system=='darwin'):
    libext = '.dylib'
libpath = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, libpath)

libfile_solver = os.path.join(libpath, 'libedi2pysolver'+libext)
libfile_fit     = os.path.join(libpath, 'libedi2pyfit'+libext)

try:
    libedi2py_solver = CDLL(libfile_solver)
    libedi2py_fit     = CDLL(libfile_fit)
except:
    libedi2py_solver = None
    libedi2py_fit    = None

####################################################################
# Create the solver. class (this is what the python module sees)
####################################################################

solver = LinkSolver(libedi2py_solver)
fit    = LinkFit(libedi2py_fit)

######################################
# GLOBAL VARIABLES FOR SOLVER
######################################

try:
    add_global_variable(solver, "Nbath", c_int.in_dll(libedi2py_solver, "Nbath"), "value")
    add_global_variable(solver, "Norb", c_int.in_dll(libedi2py_solver, "Norb"), "value")
    add_global_variable(solver, "Nspin", c_int.in_dll(libedi2py_solver, "Nspin"), "value")
    add_global_variable(solver, "Nloop", c_int.in_dll(libedi2py_solver, "Nloop"), "value")
    add_global_variable(solver, "Nph", c_int.in_dll(libedi2py_solver, "Nph"), "value")
    add_global_variable(solver, "Nsuccess", c_int.in_dll(libedi2py_solver, "Nsuccess"), "value")
    add_global_variable(solver, "Lmats", c_int.in_dll(libedi2py_solver, "Lmats"), "value")
    add_global_variable(solver, "Lreal", c_int.in_dll(libedi2py_solver, "Lreal"), "value")
    add_global_variable(solver, "Ltau", c_int.in_dll(libedi2py_solver, "Ltau"), "value")
    add_global_variable(solver, "Lpos", c_int.in_dll(libedi2py_solver, "Lpos"), "value")
    add_global_variable(solver, "LOGfile", c_int.in_dll(libedi2py_solver, "LOGfile"), "value")

    add_global_variable(solver, "Uloc", ARRAY(c_double, 5).in_dll(libedi2py_solver, "Uloc"), "value")
    add_global_variable(solver, "Ust", c_double.in_dll(libedi2py_solver, "Ust"), "value")
    add_global_variable(solver, "Jh", c_double.in_dll(libedi2py_solver, "Jh"), "value")
    add_global_variable(solver, "Jx", c_double.in_dll(libedi2py_solver, "Jx"), "value")
    add_global_variable(solver, "Jp", c_double.in_dll(libedi2py_solver, "Jp"), "value")
    add_global_variable(solver, "xmu", c_double.in_dll(libedi2py_solver, "xmu"), "value")
    add_global_variable(solver, "beta", c_double.in_dll(libedi2py_solver, "beta"), "value")
    add_global_variable(solver, "dmft_error", c_double.in_dll(libedi2py_solver, "dmft_error"), "value")
    add_global_variable(solver, "eps", c_double.in_dll(libedi2py_solver, "eps"), "value")
    add_global_variable(solver, "wini", c_double.in_dll(libedi2py_solver, "wini"), "value")
    add_global_variable(solver, "wfin", c_double.in_dll(libedi2py_solver, "wfin"), "value")
    add_global_variable(solver, "xmin", c_double.in_dll(libedi2py_solver, "xmin"), "value")
    add_global_variable(solver, "xmax", c_double.in_dll(libedi2py_solver, "xmax"), "value")
    add_global_variable(solver, "sb_field", c_double.in_dll(libedi2py_solver, "sb_field"), "value")
    add_global_variable(solver, "nread", c_double.in_dll(libedi2py_solver, "nread"), "value")

    add_global_variable(solver, "ed_total_ud", c_bool.in_dll(libedi2py_solver, "ed_total_ud"), "value")
    add_global_variable(solver, "ed_twin", c_bool.in_dll(libedi2py_solver, "ed_twin"), "value")
except:
    print("Could not load edipy2solver library. Is the python module correctly installed?")
    
######################################
# GLOBAL VARIABLES FOR FIT
######################################

try:
    add_global_variable(fit, "Nbath", c_int.in_dll(libedi2py_fit, "Nbath"), "value")
    add_global_variable(fit, "Norb", c_int.in_dll(libedi2py_fit, "Norb"), "value")
    add_global_variable(fit, "Nspin", c_int.in_dll(libedi2py_fit, "Nspin"), "value")
    add_global_variable(fit, "Nloop", c_int.in_dll(libedi2py_fit, "Nloop"), "value")
    add_global_variable(fit, "Nph", c_int.in_dll(libedi2py_fit, "Nph"), "value")
    add_global_variable(fit, "Nsuccess", c_int.in_dll(libedi2py_fit, "Nsuccess"), "value")
    add_global_variable(fit, "Lmats", c_int.in_dll(libedi2py_fit, "Lmats"), "value")
    add_global_variable(fit, "Lreal", c_int.in_dll(libedi2py_fit, "Lreal"), "value")
    add_global_variable(fit, "Ltau", c_int.in_dll(libedi2py_fit, "Ltau"), "value")
    add_global_variable(fit, "Lfit", c_int.in_dll(libedi2py_fit, "Lfit"), "value")
    add_global_variable(fit, "Lpos", c_int.in_dll(libedi2py_fit, "Lpos"), "value")
    add_global_variable(fit, "LOGfile", c_int.in_dll(libedi2py_fit, "LOGfile"), "value")

    add_global_variable(fit, "xmu", c_double.in_dll(libedi2py_fit, "xmu"), "value")
    add_global_variable(fit, "beta", c_double.in_dll(libedi2py_fit, "beta"), "value")
    add_global_variable(fit, "dmft_error", c_double.in_dll(libedi2py_fit, "dmft_error"), "value")
    add_global_variable(fit, "eps", c_double.in_dll(libedi2py_fit, "eps"), "value")
    add_global_variable(fit, "wini", c_double.in_dll(libedi2py_fit, "wini"), "value")
    add_global_variable(fit, "wfin", c_double.in_dll(libedi2py_fit, "wfin"), "value")
    add_global_variable(fit, "xmin", c_double.in_dll(libedi2py_fit, "xmin"), "value")
    add_global_variable(fit, "xmax", c_double.in_dll(libedi2py_fit, "xmax"), "value")
except:
    print("Could not load edipy2fit library. Is the python module correctly installed?")


######################################
# GLOBAL FUNCTIONS FOR SOLVER
######################################

#from here
solver.get_bath_type = types.MethodType(get_bath_type, solver)
solver.get_ed_mode = types.MethodType(get_ed_mode, solver)

#read_input
import func_read_input
solver.read_input = types.MethodType(func_read_input.read_input, solver)

#aux_funx
import func_aux_funx
solver.set_hloc = types.MethodType(func_aux_funx.set_hloc, solver)
solver.search_variable = types.MethodType(func_aux_funx.search_variable, solver)
solver.check_convergence = types.MethodType(func_aux_funx.check_convergence, solver)

#bath
import func_bath
solver.get_bath_dimension = types.MethodType(func_bath.get_bath_dimension, solver)
solver.set_hreplica = types.MethodType(func_bath.set_hreplica, solver)
solver.set_hgeneral = types.MethodType(func_bath.set_hgeneral, solver)
solver.break_symmetry_bath = types.MethodType(func_bath.break_symmetry_bath, solver)
solver.spin_symmetrize_bath = types.MethodType(func_bath.spin_symmetrize_bath, solver)
solver.orb_symmetrize_bath = types.MethodType(func_bath.orb_symmetrize_bath, solver)
solver.orb_equality_bath = types.MethodType(func_bath.orb_equality_bath, solver)
solver.ph_symmetrize_bath = types.MethodType(func_bath.ph_symmetrize_bath, solver)
solver.save_array_as_bath = types.MethodType(func_bath.save_array_as_bath, solver)

solver.bath_inspect = types.MethodType(func_bath.bath_inspect, solver)


#main
import func_main
solver.init_solver = types.MethodType(func_main.init_solver, solver)
solver.solve = types.MethodType(func_main.solve, solver)
solver.finalize_solver = types.MethodType(func_main.finalize_solver, solver)

#io
import func_io
solver.get_sigma = types.MethodType(func_io.get_sigma, solver)
solver.get_gimp = types.MethodType(func_io.get_gimp, solver)
solver.get_g0and= types.MethodType(func_io.get_g0and, solver)
solver.get_delta= types.MethodType(func_io.get_delta, solver)
solver.get_dens = types.MethodType(func_io.get_dens, solver)
solver.get_mag = types.MethodType(func_io.get_mag, solver)
solver.get_docc = types.MethodType(func_io.get_docc, solver)
solver.get_eimp = types.MethodType(func_io.get_eimp, solver)
solver.build_sigma = types.MethodType(func_io.build_sigma, solver)
solver.build_gimp = types.MethodType(func_io.build_gimp, solver)


######################################
# GLOBAL FUNCTIONS FOR FIT
######################################

#read_input
import func_read_input
fit.read_input = types.MethodType(func_read_input.read_input, fit)

#bath_fit
import func_bath_fit
fit.chi2_fitgf = types.MethodType(func_bath_fit.chi2_fitgf, fit)

#aux_funx
import func_aux_funx
fit.check_convergence = types.MethodType(func_aux_funx.check_convergence, fit)

#set hreplica, hgeneral
import func_bath
fit.set_hreplica = types.MethodType(func_bath.set_hreplica, fit)
fit.set_hgeneral = types.MethodType(func_bath.set_hgeneral, fit)

import func_aux_funx
fit.set_hloc = types.MethodType(func_aux_funx.set_hloc, fit)
