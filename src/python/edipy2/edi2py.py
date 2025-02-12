from ctypes import *
import ctypes
import numpy as np
import os,sys
import types

#################################
#Classes
#################################

#dummy classes, to be filled
class LinkSolver:
    def __init__(self,library):
        self.library = CDLL(library,mode=RTLD_LOCAL)
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
     
    @property
    def nspin(self):
        targetobj = c_int.in_dll(self.library, "Nspin")
        return getattr(targetobj, "value")

    @nspin.setter
    def nspin(self,value):
        targetobj = c_int.in_dll(self.library, "Nspin")
        setattr(targetobj, "value", value)

#    def add_global_variable(self, dynamic_name, target_object_type, target_attribute="value"):
#        target_object = target_object_type.in_dll(self.library, dynamic_name)
#        
#        def _get_from_dll(self):
#            try:
#                attrib = getattr(target_object, target_attribute)
#                try: #this is for strings
#                    attrib = attrib.decode()
#                except:
#                    pass
#            except: #this is for arrays
#                if(len(target_object)>1):
#                    return [target_object[x] for x in range(len(target_object))]
#            return attrib

#        def _set_to_dll(self, new_value):
#            try: #this is for arrays
#                if(len(target_object)>1):
#                    minlength=min(len(target_object),len(new_value))
#                    target_object[0:minlength]=new_value[0:minlength]
#            except:
#                try:
#                    new_value = new_value.encode()
#                except:
#                    pass
#                setattr(target_object, target_attribute, new_value)

#################################
#AUXILIARY FUNCTIONS
#################################




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


####################################################################
# Create the solver. class (this is what the python module sees)
####################################################################

solver1 = LinkSolver(libfile_solver)
solver2 = LinkSolver(libfile_solver)

#solver1.add_global_variable("Nbath", c_int)
#solver2.add_global_variable("Nbath", c_int)


