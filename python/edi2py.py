from ctypes import *
import numpy as np
import os,sys
import types

#################################
#link to access global variables
#################################

#dummy class, to be filled
class Link:
    pass

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
libfile = os.path.join(libpath, 'libedi2py'+libext)
libedi2py = CDLL(libfile)


######################################
# READ_INPUT
######################################

read_input_wrap = libedi2py.read_input
read_input_wrap.argtypes = [c_char_p]  
read_input_wrap.restype = None

def read_input(self,input_string):
    c_string = c_char_p(input_string.encode())
    read_input_wrap(c_string)
    


####################################################################
# Create the global_env class (this is what the python module sees)
####################################################################

global_env=Link()

#variables
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


#functions
global_env.read_input = types.MethodType(read_input, global_env)

