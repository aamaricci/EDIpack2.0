EDIpack 2.0: A massively parallel Exact Diagonalization solver for quantum Impurity problems
============================================================================================================

(under construction)

EDIpack2.0 is a suitable extension of `EDIpack
<https://github.com/aamaricci/EDIpack>`_ : a  Lanczos based method 
for the solution of generic Quantum Impurity problems,  exploiting distributed memory MPI parallelisation.
This 2.0 version, aims to solve single-site, multi-orbital models, in either  *normal*, *superconducting* (s-wave) or *Spin-non-conserving* (e.g. with Spin-Orbit Coupling or in-plane magnetization) phases, including electron-phonons coupling. The code works at zero and low temperatures.   
 
See
`j.cpc.2021.108261 <https://doi.org/10.1016/j.cpc.2021.108261>`_
for further information about the underlying algorithm in the *normal*
channel. 
Substantial modifications have been developed to address the Superconducting and non-SU(2) channels.  
An updated manuscript will be released soon. 


.. toctree::
   :maxdepth: 2

   fortran_src/edipack2
   python_api/edipy2



   

