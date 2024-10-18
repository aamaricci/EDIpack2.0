EDIpack2.0
################

A massively parallel Exact Diagonalization solver for quantum Impurity problems.
***************************************************************************************************************

(under construction)

Introduction
===============

EDIpack2.0 is a Lanczos based method 
for the solution of generic Quantum Impurity problems,  exploiting distributed memory MPI parallelisation.
The 2.0 version extends the `EDIpack <https://github.com/aamaricci/EDIpack>`_ by enabling the solution of
single-site, multi-orbital models with different conserved
quantum numbers :math:`\vec{Q}` corresponding to separate operational modes:

* :math:`\vec{Q}=[\vec{N}_\uparrow,\vec{N}_\downarrow]` for which
  either the number of total or orbital spin up and down electrons is
  conserved: **normal**: 
* :math:`\vec{Q}=S_z`  with *s*-wave pairing with conserved total
  magnetization:  **superconducting**  
* :math:`\vec{Q}=N_{\rm tot}`  where spin degrees
  freedom is not fully conserved:  **non-SU(2)**
  (used, for instance, in presence of local Spin-Orbit coupling, in-plane magnetization or in-plane triplet excitonic condensation)

All modes include electron-phonon coupling (local or Holstein
phonons). EDIpack2.0  is designed to obtain the lowest part of the
spectrum of the problem, thus it naturally works at zero temperature
but can also be used to explore low temperature properties.  
 
The EDIpack2.0 diagonalization algorithm is based on a massively
parallel execution of matrix-vector products, required in the context
of Lanczos-Arnoldi linear procedures.  See `j.cpc.2021.108261
<https://doi.org/10.1016/j.cpc.2021.108261>`_ for a detailed
descriptions of these algorithms.
However, substantial modifications have been introduced in the 2.0
version to address the *Superconducting* and *non-SU(2)* channels.  
An updated manuscript will be released soon. 




.. toctree::
   :maxdepth: 2

   fortran_src/edipack2
   python_api/edipy2



   

