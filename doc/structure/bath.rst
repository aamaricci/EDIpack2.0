Bath Auxiliary
###########################

In this set of modules we implement a number of auxiliary procedures
which are required to enumerate the bath levels, performs all the
relevant checks on the user input bath, build up of the
:code:`replica/general` matrix basis or apply given symmetry operation
on the user side. 

.. toctree::
   :maxdepth: 1

   bath_aux/00_ed_vars_global
   bath_aux/01_ed_bath_aux
   bath_aux/02_ed_bath_dim
   bath_aux/03_ed_bath_user
   bath_aux/04_ed_bath_replica
	      
Bath DMFT
###########################

In :f:mod:`ed_bath_dmft` we implement operations on the  :f:type:`effective_bath` data
structure: a suitable representation of the effective bath used
internally in the code. We refer to the generic shared instance of this bath
as :f:var:`dmft_bath`.  Depending on the value :f:var:`bath_type` this quantity collect
different bath parameters which can be directly accessed in the
construction of symmetry sectors Hamiltonian.  

.. toctree::
   :maxdepth: 1

   bath_dmft/01_ed_bath_dmft


Bath Functions
###########################

In :f:mod:`ed_bath_functions` we implement on-the-fly construction of
the hybridization functions :math:`\Delta(z) = \sum_p
\hat{V}^p\left[z-\hat{h}^p \right]^{-1}\hat{V}^p`, as well as of the
non-interacting Anderson Green's functions
:math:`G_0(z) = \left[z +\mu - H_{loc} - \Delta(z) \right]^{-1}` for
all different cases selected by :f:var:`ed_mode` and :f:var:`bath_type`.

.. toctree::
   :maxdepth: 1

   bath_functions/ed_bath_functions

   
