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
libfile_fit     = os.path.join(libpath, 'libedi2pysolver'+libext)

try:
    libedi2py_solver = CDLL(libfile_solver)
    libedi2py_fit     = CDLL(libfile_fit)
except:
    libedi2py_solver = None
    libedi2py_fit    = None

####################################################################
# Create the global_env class (this is what the python module sees)
####################################################################

solver = LinkSolver(libedi2py_solver)
fit    = LinkFit(libedi2py_fit)

######################################
# GLOBAL VARIABLES FOR SOLVER
######################################

try:
    add_global_variable(global_env, "Nbath", c_int.in_dll(libedi2py_solver, "Nbath"), "value")
    add_global_variable(global_env, "Norb", c_int.in_dll(libedi2py_solver, "Norb"), "value")
    add_global_variable(global_env, "Nspin", c_int.in_dll(libedi2py_solver, "Nspin"), "value")
    add_global_variable(global_env, "Nloop", c_int.in_dll(libedi2py_solver, "Nloop"), "value")
    add_global_variable(global_env, "Nph", c_int.in_dll(libedi2py_solver, "Nph"), "value")
    add_global_variable(global_env, "Nsuccess", c_int.in_dll(libedi2py_solver, "Nsuccess"), "value")
    add_global_variable(global_env, "Lmats", c_int.in_dll(libedi2py_solver, "Lmats"), "value")
    add_global_variable(global_env, "Lreal", c_int.in_dll(libedi2py_solver, "Lreal"), "value")
    add_global_variable(global_env, "Ltau", c_int.in_dll(libedi2py_solver, "Ltau"), "value")
    add_global_variable(global_env, "Lpos", c_int.in_dll(libedi2py_solver, "Lpos"), "value")
    add_global_variable(global_env, "LOGfile", c_int.in_dll(libedi2py_solver, "LOGfile"), "value")

    add_global_variable(global_env, "Uloc", ARRAY(c_double, 5).in_dll(libedi2py_solver, "Uloc"), "value")
    add_global_variable(global_env, "Ust", c_double.in_dll(libedi2py_solver, "Ust"), "value")
    add_global_variable(global_env, "Jh", c_double.in_dll(libedi2py_solver, "Jh"), "value")
    add_global_variable(global_env, "Jx", c_double.in_dll(libedi2py_solver, "Jx"), "value")
    add_global_variable(global_env, "Jp", c_double.in_dll(libedi2py_solver, "Jp"), "value")
    add_global_variable(global_env, "xmu", c_double.in_dll(libedi2py_solver, "xmu"), "value")
    add_global_variable(global_env, "beta", c_double.in_dll(libedi2py_solver, "beta"), "value")
    add_global_variable(global_env, "dmft_error", c_double.in_dll(libedi2py_solver, "dmft_error"), "value")
    add_global_variable(global_env, "eps", c_double.in_dll(libedi2py_solver, "eps"), "value")
    add_global_variable(global_env, "wini", c_double.in_dll(libedi2py_solver, "wini"), "value")
    add_global_variable(global_env, "wfin", c_double.in_dll(libedi2py_solver, "wfin"), "value")
    add_global_variable(global_env, "xmin", c_double.in_dll(libedi2py_solver, "xmin"), "value")
    add_global_variable(global_env, "xmax", c_double.in_dll(libedi2py_solver, "xmax"), "value")
    add_global_variable(global_env, "sb_field", c_double.in_dll(libedi2py_solver, "sb_field"), "value")
    add_global_variable(global_env, "nread", c_double.in_dll(libedi2py_solver, "nread"), "value")

    add_global_variable(global_env, "ed_total_ud", c_bool.in_dll(libedi2py_solver, "ed_total_ud"), "value")
    add_global_variable(global_env, "ed_twin", c_bool.in_dll(libedi2py_solver, "ed_twin"), "value")
except:
    print("Could not load edipy2solver library. Is the python module correctly installed?")
    
######################################
# GLOBAL VARIABLES FOR FIT
######################################

try:
    add_global_variable(global_env, "Nbath", c_int.in_dll(libedi2py_fit, "Nbath"), "value")
    add_global_variable(global_env, "Norb", c_int.in_dll(libedi2py_fit, "Norb"), "value")
    add_global_variable(global_env, "Nspin", c_int.in_dll(libedi2py_fit, "Nspin"), "value")
    add_global_variable(global_env, "Nloop", c_int.in_dll(libedi2py_fit, "Nloop"), "value")
    add_global_variable(global_env, "Nph", c_int.in_dll(libedi2py_fit, "Nph"), "value")
    add_global_variable(global_env, "Nsuccess", c_int.in_dll(libedi2py_fit, "Nsuccess"), "value")
    add_global_variable(global_env, "Lmats", c_int.in_dll(libedi2py_fit, "Lmats"), "value")
    add_global_variable(global_env, "Lreal", c_int.in_dll(libedi2py_fit, "Lreal"), "value")
    add_global_variable(global_env, "Ltau", c_int.in_dll(libedi2py_fit, "Ltau"), "value")
    add_global_variable(global_env, "Lfit", c_int.in_dll(libedi2py_fit, "Lfit"), "value")
    add_global_variable(global_env, "Lpos", c_int.in_dll(libedi2py_fit, "Lpos"), "value")
    add_global_variable(global_env, "LOGfile", c_int.in_dll(libedi2py_fit, "LOGfile"), "value")

    add_global_variable(global_env, "xmu", c_double.in_dll(libedi2py_fit, "xmu"), "value")
    add_global_variable(global_env, "beta", c_double.in_dll(libedi2py_fit, "beta"), "value")
    add_global_variable(global_env, "dmft_error", c_double.in_dll(libedi2py_fit, "dmft_error"), "value")
    add_global_variable(global_env, "eps", c_double.in_dll(libedi2py_fit, "eps"), "value")
    add_global_variable(global_env, "wini", c_double.in_dll(libedi2py_fit, "wini"), "value")
    add_global_variable(global_env, "wfin", c_double.in_dll(libedi2py_fit, "wfin"), "value")
    add_global_variable(global_env, "xmin", c_double.in_dll(libedi2py_fit, "xmin"), "value")
    add_global_variable(global_env, "xmax", c_double.in_dll(libedi2py_fit, "xmax"), "value")
except:
    print("Could not load edipy2fit library. Is the python module correctly installed?")


######################################
# GLOBAL FUNCTIONS FOR SOLVER
######################################

#from here
global_env.get_bath_type = types.MethodType(get_bath_type, solver)
global_env.get_ed_mode = types.MethodType(get_ed_mode, solver)

#read_input
import func_read_input
global_env.read_input = types.MethodType(func_read_input.read_input, solver)

#aux_funx
import func_aux_funx
global_env.set_hloc = types.MethodType(func_aux_funx.set_hloc, solver)
global_env.search_variable = types.MethodType(func_aux_funx.search_variable, solver)
global_env.check_convergence = types.MethodType(func_aux_funx.check_convergence, solver)

#bath
import func_bath
global_env.get_bath_dimension = types.MethodType(func_bath.get_bath_dimension, solver)
global_env.set_hreplica = types.MethodType(func_bath.set_hreplica, solver)
global_env.set_hgeneral = types.MethodType(func_bath.set_hgeneral, solver)
global_env.break_symmetry_bath = types.MethodType(func_bath.break_symmetry_bath, solver)
global_env.spin_symmetrize_bath = types.MethodType(func_bath.spin_symmetrize_bath, solver)
global_env.orb_symmetrize_bath = types.MethodType(func_bath.orb_symmetrize_bath, solver)
global_env.orb_equality_bath = types.MethodType(func_bath.orb_equality_bath, solver)
global_env.ph_symmetrize_bath = types.MethodType(func_bath.ph_symmetrize_bath, solver)
global_env.save_array_as_bath = types.MethodType(func_bath.save_array_as_bath, solver)

global_env.bath_inspect = types.MethodType(func_bath.bath_inspect, solver)


#main
import func_main
global_env.init_solver = types.MethodType(func_main.init_solver, solver)
global_env.solve = types.MethodType(func_main.solve, solver)
global_env.finalize_solver = types.MethodType(func_main.finalize_solver, solver)

#io
import func_io
global_env.get_sigma = types.MethodType(func_io.get_sigma, solver)
global_env.get_gimp = types.MethodType(func_io.get_gimp, solver)
global_env.get_g0and= types.MethodType(func_io.get_g0and, solver)
global_env.get_delta= types.MethodType(func_io.get_delta, solver)
global_env.get_dens = types.MethodType(func_io.get_dens, solver)
global_env.get_mag = types.MethodType(func_io.get_mag, solver)
global_env.get_docc = types.MethodType(func_io.get_docc, solver)
global_env.get_eimp = types.MethodType(func_io.get_eimp, solver)
global_env.build_sigma = types.MethodType(func_io.build_sigma, solver)
global_env.build_gimp = types.MethodType(func_io.build_gimp, solver)


######################################
# GLOBAL FUNCTIONS FOR FIT
######################################

#read_input
import func_read_input
global_env.read_input = types.MethodType(func_read_input.read_input, fit)

#bath_fit
import func_bath_fit
global_env.chi2_fitgf = types.MethodType(func_bath_fit.chi2_fitgf, fit)

#aux_funx
import func_aux_funx
global_env.check_convergence = types.MethodType(func_aux_funx.check_convergence, fit)
