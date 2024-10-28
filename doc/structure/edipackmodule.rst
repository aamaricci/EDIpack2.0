EDIpack2
==========================


:f:mod:`edipack2` is the top module of the `EDIpack2.0` code. It provides access to the all the relevant procedures of the library, realizing the Fortran API. The user needs to invoke use of this module to get access to `EDIpack2.0` as:

   .. code-block:: fortran

      program test
          USE EDIPACK2
	  ...

   		   
The module also contains a subset of the global and input variables that can be accessed in the userspace. 

.. f:automodule::   edipack2
