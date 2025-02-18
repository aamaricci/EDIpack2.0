from ctypes import *
import numpy as np
import os, sys
import types


def read_input(self, input_string):
    """

       This function reads from the input file of EDIpack2. If the file does not \
       exist, a template file is generated with default parameters.
       This is generated with the prefix "used." which will need to be \
       removed for it to be read. "used.${input_string}" will be updated within
       the DMFT loop with the current value of the input variables.

       :param input_string: The name of the input file to be read, including 
        the extension
       :type input_string: str
       :return: Nothing
       :rtype: None

    """

    read_input_wrap = self.library.read_input
    read_input_wrap.argtypes = [c_char_p]
    read_input_wrap.restype = None
    c_string = c_char_p(input_string.encode())
    read_input_wrap(c_string)
