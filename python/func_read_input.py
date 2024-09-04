from ctypes import *
import numpy as np
import os,sys
import types

def read_input(self,input_string):
    read_input_wrap = self.library.read_input
    read_input_wrap.argtypes = [c_char_p]  
    read_input_wrap.restype = None
    c_string = c_char_p(input_string.encode())
    read_input_wrap(c_string)
    return "lol"