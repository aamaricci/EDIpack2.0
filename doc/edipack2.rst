Here we give an overview of the structure of the EDIpack2.0 library,
with a detailed description of the relevant modules and procedures.


General 
###########################

This include a set of modules which contains input variables
:f:mod:`ed_input_vars`, global ones and shared classes/data structure in
:f:mod:`ed_vars_global`, global memory (de)allocation and evaluation
of symmetry sectors dimensions :f:mod:`ed_setup` as well as
generic functions which are used throughout the code :f:mod:`ed_aux_funx`.

:f:mod:`ed_vars_global`  contains the definition of simple data
structures, such as  :f:type:`gfmatrix`  storing all weights and poles of the Green's functions or
:f:type:`effective_bath` gathering the different bath components
according to value of :f:var:`bath_type` and :f:var:`ed_mode`.
   

.. toctree::
   :maxdepth: 2
   :glob:

   structure/general/*


Sparse Matrix
###########################

The module  :f:mod:`ed_sparse_matrix` contains the implementation of a
suitable Compact Row Stored sparse matrix, used to store sparse Hamiltonians in each symmetry
sector. The data structure essentially corresponds to an array of
rows as dynamical vectors, one for the values and one for the columns index of
the non-zero elements of the matrix. This enables :math:`O(1)`
access. If MPI parallelization is enabled the structure includes an
automatic rows split balanced among the different threads and a
further subdivision of the columns/values in local and non-local
blocks.   



.. toctree::
   :maxdepth: 2

   structure/classes/01_ed_sparse_matrix


   

EigenSpace
###########################

In  :f:mod:`ed_eigenspace` we implement an ordered
single linked list to efficiently store the lower part of the energy
spectrum, including eigenvectors, eigenvalues, symmetry sector index and
twin states, i.e. states in sectors with exchanged quantum numbers which do
have same energy. 


.. toctree::
   :maxdepth: 2

   structure/classes/02_ed_eigenspace

   
Sectors
###########################

The direct diagonalization of the Hamiltonian for a finite system
becomes quickly impossible due to the exponential increase of its
size :math:`D^{Ns}`, where :math:`D` is the dimension of the local
Hilbert space (:math:`D=4` for spin-:math:`1/2` electrons).

A way to partially soften such severe behavior is to take into account
conserved quantities. For any quantity :math:`{\cal Q}` such that
:math:`[H,{\cal Q}]=0`,  we can separate the structure of the
Hamiltonian matrix into different blocks, each corresponding to a
discrete value :math:`q` of the spectrum of :math:`Q`.
Analyzing one-by-one such symmetry sector it becomes possible to
construct at least part of the energy spectrum.

In :f:mod:`ed_sector` we implement the construction of the symmetry
sectors for the three cases considered in `EDIpack2.0`, that we
recall are:

* :math:`\vec{Q}=[\vec{N}_\uparrow,\vec{N}_\downarrow]` for which
  either the number of total or orbital spin up and down electrons is
  conserved: **NORMAL** 

* :math:`\vec{Q}=S_z`  with conserved total magnetization:  **SUPERConducting**  

* :math:`\vec{Q}=N_{\rm tot}`  where spin degrees freedom is not fully
  conserved:  **NON-SU(2)**
  

.. toctree::
   :maxdepth: 1
   :glob:

   structure/ed_sector



Bath 
###########################

The construction and the handling of the bath is a crucial part of the
description of the generic single impurity Anderson problem.
In `EDIpack2.0` we implemented different bath topologies, which can be
selected using the variable :f:var:`bath_type` = 
:code:`normal, hybrid, replica, general` according to the nature of the problem at hand.

.. toctree::
   :maxdepth: 1
   :glob:

   structure/bath
   






   


Exact Diagonalization
###########################

This part of the `EDIpack2.0` code implements the exact
diagonalization of the general, single-site, multi-orbital quantum
impurity problem using  different operational modes,  corresponding to the
choice of the specific symmetry implemented in the code, i.e. which
quantum numbers are to be conserved. The operational modes are selected
by the variable :f:var:`ed_mode` =  :code:`normal, superc, nosu2`. See
:f:mod:`ed_sector` for more info about the symmetries implemented in
the code.

To keep the code simple we implemented the three different channels in
distinct class of modules, essentially performing all the main
operations required for the construction of the sector Hamiltonian,
their diagonalization, the evaluation of the impurity Green's functions,
the impurity susceptibilities and observables.  

.. toctree::
   :maxdepth: 1
   :glob:

   structure/diag
   



Green's functions 
###########################
This part of the `EDIpack2.0` code implements the calculation of the
impurity interacting Green's functions, self-energy functions and
impurity susceptibilities. Calculations are performed in  each operational mode,  corresponding to the
choice of the specific symmetry implemented in the code, i.e. which
quantum numbers are to be conserved. The operational modes are selected
by the variable :f:var:`ed_mode` =  :code:`normal, superc, nosu2`. See
:f:mod:`ed_sector` for more info about the symmetries implemented in
the code.


.. toctree::
   :maxdepth: 1
   :glob:

   structure/greensfunctions


Observables
###########################

This part of the `EDIpack2.0` code implements the calculation of the
impurity observables and static correlations, such as density,
internal energy or double occupation. Calculations are performed in
each operational mode,  corresponding to the choice of the specific
symmetry implemented in the code, i.e. which quantum numbers are to be
conserved. The operational modes are selected by the variable
:f:var:`ed_mode` =  :code:`normal, superc, nosu2`. See
:f:mod:`ed_sector` for more info about the symmetries implemented
in the code.


.. toctree::
   :maxdepth: 1
   :glob:

   structure/observables




Input/Output
###########################

This module provides access to the results of the exact
diagonalization. All quantities such as dynamical response functions,
self-energy components  or impurity observables  can be retrieved  
using specific functions. Additionally we provide procedure to perform
on-the-fly re-calculation of the impurity Green's functions and
self-energy on a given arbitrary set of points in the complex
frequency domain.    

.. toctree::
   :maxdepth: 1

   structure/io/ed_io


   
EDIpack2
###########################

.. toctree::
   :maxdepth: 1
   :glob:

   structure/edipackmodule


..
   EDIpack2.0 Library
   #############################

   An overview of the structure of the `EDIpack2.0` library with a detailed description of the relevant modules.


   .. toctree::
      :maxdepth: 1
      :glob:

      structure
      




