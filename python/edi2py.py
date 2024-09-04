from ctypes import *
import numpy as np
import os,sys
import types

#################################
#AUXILIARY FUNCTIONS
#################################

#dummy class, to be filled
class Link:
    def __init__(self,library):
        self.library = library

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


######################################
# Load shared library with C-bindings
######################################

system = sys.platform
libext = '.so' 
if(system=='darwin'):
    libext = '.dylib'
libpath = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, libpath)
libfile = os.path.join(libpath, 'libedi2py'+libext)
libedi2py = CDLL(libfile)

####################################################################
# Create the global_env class (this is what the python module sees)
####################################################################

global_env=Link(libedi2py)

######################################
# GLOBAL VARIABLES
######################################

add_global_variable(global_env, "Nbath", c_int.in_dll(libedi2py, "Nbath"), "value")
add_global_variable(global_env, "Norb", c_int.in_dll(libedi2py, "Norb"), "value")
add_global_variable(global_env, "Nspin", c_int.in_dll(libedi2py, "Nspin"), "value")
add_global_variable(global_env, "Nloop", c_int.in_dll(libedi2py, "Nloop"), "value")
add_global_variable(global_env, "Nph", c_int.in_dll(libedi2py, "Nph"), "value")
add_global_variable(global_env, "Nsuccess", c_int.in_dll(libedi2py, "Nsuccess"), "value")
add_global_variable(global_env, "Lmats", c_int.in_dll(libedi2py, "Lmats"), "value")
add_global_variable(global_env, "Lreal", c_int.in_dll(libedi2py, "Lreal"), "value")
add_global_variable(global_env, "Ltau", c_int.in_dll(libedi2py, "Ltau"), "value")
add_global_variable(global_env, "Lpos", c_int.in_dll(libedi2py, "Lpos"), "value")
add_global_variable(global_env, "LOGfile", c_int.in_dll(libedi2py, "LOGfile"), "value")

add_global_variable(global_env, "Uloc", ARRAY(c_double, 5).in_dll(libedi2py, "Uloc"), "value")
add_global_variable(global_env, "Ust", c_double.in_dll(libedi2py, "Ust"), "value")
add_global_variable(global_env, "Jh", c_double.in_dll(libedi2py, "Jh"), "value")
add_global_variable(global_env, "Jx", c_double.in_dll(libedi2py, "Jx"), "value")
add_global_variable(global_env, "Jp", c_double.in_dll(libedi2py, "Jp"), "value")
add_global_variable(global_env, "xmu", c_double.in_dll(libedi2py, "xmu"), "value")
add_global_variable(global_env, "beta", c_double.in_dll(libedi2py, "beta"), "value")
add_global_variable(global_env, "dmft_error", c_double.in_dll(libedi2py, "dmft_error"), "value")
add_global_variable(global_env, "eps", c_double.in_dll(libedi2py, "eps"), "value")
add_global_variable(global_env, "wini", c_double.in_dll(libedi2py, "wini"), "value")
add_global_variable(global_env, "wfin", c_double.in_dll(libedi2py, "wfin"), "value")
add_global_variable(global_env, "xmin", c_double.in_dll(libedi2py, "xmin"), "value")
add_global_variable(global_env, "xmax", c_double.in_dll(libedi2py, "xmax"), "value")
add_global_variable(global_env, "sb_field", c_double.in_dll(libedi2py, "sb_field"), "value")
add_global_variable(global_env, "nread", c_double.in_dll(libedi2py, "nread"), "value")

######################################
# GLOBAL FUNCTIONS
######################################

#read_input
import func_read_input
global_env.read_input = types.MethodType(func_read_input.read_input, global_env)

#aux_funx
import func_aux_funx
global_env.set_hloc = types.MethodType(func_aux_funx.set_hloc, global_env)
global_env.search_variable = types.MethodType(func_aux_funx.search_variable, global_env)
global_env.check_convergence = types.MethodType(func_aux_funx.check_convergence, global_env)

#bath
import func_bath
global_env.get_bath_dimension = types.MethodType(func_bath.get_bath_dimension, global_env)
global_env.set_Hreplica = types.MethodType(func_bath.set_Hreplica, global_env)
global_env.break_symmetry_bath = types.MethodType(func_bath.break_symmetry_bath, global_env)
global_env.spin_symmetrize_bath = types.MethodType(func_bath.spin_symmetrize_bath, global_env)
global_env.orb_symmetrize_bath = types.MethodType(func_bath.orb_symmetrize_bath, global_env)
global_env.orb_equality_bath = types.MethodType(func_bath.orb_equality_bath, global_env)
global_env.ph_symmetrize_bath = types.MethodType(func_bath.ph_symmetrize_bath, global_env)

#main
import func_main
global_env.init_solver = types.MethodType(func_main.init_solver, global_env)
global_env.solve = types.MethodType(func_main.solve, global_env)

#io
import func_io
global_env.get_sigma = types.MethodType(func_io.get_sigma, global_env)
global_env.get_gimp = types.MethodType(func_io.get_gimp, global_env)

#bath_fit
import func_bath_fit
global_env.chi2_fitgf = types.MethodType(func_bath_fit.chi2_fitgf, global_env)