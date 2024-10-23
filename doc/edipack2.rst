Here we give an overview of the structure of the EDIpack2.0 library,
with a detailed description of the relevant modules and procedures.


General 
###########################

This include a set of modules which contains input variables
:f:mod:`ed_input_vars`, global ones and shared classes/data structure in
:f:mod:`ed_vars_global`, global memory (de)allocation :f:mod:`ed_setup` as well as
generic functions which are used throughout the code :f:mod:`ed_aux_funx`.

.. toctree::
   :maxdepth: 2
   :glob:

   structure/general/*


Classes
###########################

Some simple data structures contained in :f:mod:`ed_vars_global` have
been discussed in :doc:`structure/general/ed_vars_global`, such as  :f:type:`gfmatrix` 
type storing all weights and poles of the Green's functions or
:f:type:`effective_bath` gathering the different bath components
according to value of :f:var:`bath_type` and :f:var:`ed_mode`.
However other tasks in the code requires more elaborated classes which
are implemented in two main modules.

The first is :f:mod:`ed_sparse_matrix` describing Compact Row Stored
sparse matrices used to store sparse Hamiltonians in each symmetry
sector. The data structure essentially corresponds to an array of
rows as dynamical vectors, one for the values and one for the columns index of
the non-zero elements of the matrix. This enables :math:`O(1)`
access. If MPI parallelization is enabled the structure includes an
automatic rows split balanced among the different threads and a
further subdivision of the columns/values in local and non-local
blocks.   

The second class is :mod:`ed_eigenspace` which describes an ordered
single linked list efficiently storing the lower part of the energy
spectrum, including eigenvectors, eigenvalues, symmetry sector and
twin states (states in sectors with exchanged quantum numbers which do
have same energy).


.. toctree::
   :maxdepth: 2
   :glob:

   structure/classes/*


   
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
The bath is described in terms of two set of parameters: the local
hamiltonian :math:`\hat{h}^p` and the hybridization :math:`\hat{V}^p`,
for :math:`p=1,\dots,N_{bath}`. The first describes the local
properties of each bath element (be that a single electronic level or
a more complex structure made of few levels), the second describes the
coupling with the impurity.

Depending on the nature of such two parameters we can distinguish
different cases corresponding to the value of :f:var:`bath_type`. 


For :f:var:`bath_type` = :code:`normal` a number :f:var:`nbath` of electronic
levels are coupled directly to each orbital level of the impurity
site.  Thus, :math:`\hat{V}^p=V^p_{a}\delta_{ab}` and
:math:`\hat{h}^p\equiv \epsilon^p_a\delta_{ab}` are both diagonal in
the orbital index.
For :f:var:`ed_mode` = :code:`superc` the bath includes a set of parameters
:math:`\Delta_p` describing the superconductive amplitude on each bath
level. 


For :f:var:`bath_type` = :code:`hybrid` a number :code:`nbath` of electronic
levels are all coupled to any orbital level of the impurity
site.  Thus, :math:`\hat{V}^p=V^p_{ab}` is a matrix in the orbital
index whereas :math:`\hat{h}^p\equiv \epsilon^p_a\delta_{ab}` remains
diagonal.
For :f:var:`ed_mode` = :code:`superc` the bath includes a set of parameters
:math:`\Delta_p` describing the superconductive amplitude on each bath
level. 


For :code:`bath_type=replica` a number :code:`Nbath` of copies of the
impurity structure are  coupled to the impurity itself.  Each bath
element is made of a number :math:`N_{orb}` of electronic levels,
i.e. the number of orbitals in the impurity site. Each bath element is
coupled to the impurity with an amplitude
:math:`\hat{V}^p=V^p_{a}\delta_{ab}`, while it is described by a local
Hamiltonian :math:`\hat{h}^p = \sum_{m=1}^{M} \lambda^p_m O_m`.
The set  :math:`\{O\}_m` is a user defined matrix basis for the impurity
Hamiltonian or, equally, for the local Hamiltonian of the lattice
problem. The numbers  :math:`\lambda^p_m\in{\mathbb R}` are
variational parameters.  

For :f:var:`bath_type` = :code:`general` a number :f:var:`nbath` of copies of the
impurity structure are  coupled to the impurity itself.  Each bath
element is made of a number :math:`N_{orb}` of electronic levels,
i.e. the number of orbitals in the impurity site. However, contrary to
the previous case, each bath element is
coupled to the impurity with a set of  amplitudes values
:math:`\hat{V}^p=V^p_{ab}`, e.g. one value for each degree of freedom
contained in the bath element.  

.. note::
   The enumeration of the total bath electronic levels is different
   among the different cases. This number is automatically evaluated
   upon calling :f:func:`get_bath_dimension`, see
   :doc:`structure/bath/ed_bath_dim`.


.. note::
   The :code:`replica, general` bath topologies are available also for
   the superconductive case :f:var:`ed_mode` = :code:`superc`. In this case the
   structure of the matrix basis should be set to the proper
   multi-orbital Nambu basis, so that off-diagonal blocks corresponds
   to anomalous components.

   

In :doc:`structure/bath/ed_bath_functions` we implement on-the-fly construction of
the Hybridization functions :math:`\Delta(z) = \sum_p
\hat{V}^p\left[z-\hat{h}^p \right]^{-1}\hat{V}^p` as well as the
non-interacting Anderson Green's functions
:math:`G_0(z) = \left[z +\mu - H_{loc} - \Delta(z) \right]^{-1}` for
all different cases.

Finally, in :doc:`structure/bath/ed_bath_fit` we provide to the
user a generic function :f:func:`ed_chi2_fitgf` performing the
minimization of a user provided Weiss field against the corresponding
model of non-interacting Anderson Green's function with the aim of
updating the bath parameters.




.. toctree::
   :maxdepth: 1
   :glob:

   structure/bath/*

Input/Output
###########################

.. toctree::
   :maxdepth: 1
   :glob:

   structure/io/*



   
Main
###########################


.. toctree::
   :maxdepth: 2
   :glob:

   structure/hamiltonian
   structure/diag
   structure/greensfunctions
   structure/chifunctions
   structure/observables


EDIpack2 FORTRAN module
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
      




