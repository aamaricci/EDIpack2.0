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
The bath is described by two set of parameters: the local
hamiltonian :math:`\hat{h}^p` and the hybridization :math:`\hat{V}^p`,
for :math:`p=1,\dots,N_{bath}`. The first describes the local
properties of each bath element (be that a single electronic level or
a more complex structure made of few levels), the second describes the
coupling with the impurity.

In `EDIpack2.0` the bath is handled using a Reverse Communication Strategy.
The user accesses the bath as a double precision rank-1 array
containing in a given order all the parameters. This array is passed
as input to `EDIpack2.0` procedures and dumped into an internal data
structure, implemented in :f:mod:`ed_bath_dmft`.

We implemented different bath topologies, which can be
selected using the variable :f:var:`bath_type` = 
:code:`normal, hybrid, replica, general` .


For :f:var:`bath_type` = :code:`normal` (:code:`hybrid`) a number :f:var:`nbath` of electronic
levels are coupled to each orbital level (to any orbital level) of the
impurity site. The bath local Hamiltonian is diagonal  
:math:`\hat{h}^p\equiv\epsilon^p_a\delta_{ab}` while the
hybridizations are:  :math:`\hat{V}^p=V^p_{a}\delta_{ab}`
(:math:`\hat{V}^p=V^p_{ab}`).  
If  :f:var:`ed_mode` = :code:`superc` the bath includes a set of parameters
:math:`\Delta_p` describing the superconductive amplitude on each bath
level. 


For :f:var:`bath_type` = :code:`replica` (:code:`general`) a number :f:var:`nbath` of copies of the
impurity structure are  coupled to the impurity itself.  Each bath
element is made of a number :math:`N_{orb}` of electronic levels,
i.e. the number of orbitals in the impurity site.
The hybridization to the impurity site is
:math:`\hat{V}^p=V^p_{a}\delta_{ab}` (:math:`\hat{V}^p=V^p_{ab}`).  
The local bath  Hamiltonian is :math:`\hat{h}^p = \sum_{m=1}^{M} \lambda^p_m O_m`.
The set  :math:`\{O\}_m` is a user defined matrix basis for the impurity
Hamiltonian or, equally, for the local Hamiltonian of the lattice
problem. The numbers  :math:`\lambda^p_m\in{\mathbb R}` are
variational parameters.  

.. note::
   The enumeration of the total bath electronic levels is different
   among the different cases. This number is automatically evaluated
   upon calling :f:func:`get_bath_dimension`, see
   :f:mod:`ed_bath_dim`.


.. note::
   The :code:`replica, general` bath topologies are available also for
   the superconductive case :f:var:`ed_mode` = :code:`superc`. In this case the
   structure of the matrix basis should be set to the proper
   multi-orbital Nambu basis, so that off-diagonal blocks corresponds
   to anomalous components.


      
.. toctree::
   :maxdepth: 1
   :glob:

   structure/bath
   



Exact Diagonalization
###########################

..
 .. raw:: html
    :file:  structure/graphs/diag.html

 |

 .. raw:: html
    :file:  structure/graphs/hamiltonian.html

 |


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
   

Browse Source Code
###########################

.. toctree::
   :maxdepth: 1
   :glob:

   browsesource/module/*



..
   EDIpack2.0 Library
   #############################

   An overview of the structure of the `EDIpack2.0` library with a detailed description of the relevant modules.


   .. toctree::
      :maxdepth: 1
      :glob:

      structure
      




