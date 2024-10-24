Normal mode
###########################

This set of modules implements the exact diagonalization of the single
impurity problems assuming :math:`\vec{Q}=\left[\vec{N}_\uparrow,\vec{N}_\downarrow \right]`.
Where :math:`\vec{N}_\sigma=N_\sigma` if the total number of electrons
with spin :math:`\sigma` is conserved (:f:var:`ed_total_ud` = T ) or
:math:`\vec{N}_\sigma=[ N_{1\sigma},\dots,N_{N_{orb}\sigma} ]` if the
number of electrons in the orbital :math:`\alpha=1,\dots,N_{orb}` and
spin :math:`\sigma` is conserved (:f:var:`ed_total_ud` = F). 

This case corresponds to the normal phase in presence of spin
conservation, possibly reduced to :math:`U(1)` in presence of long
range magnetic order along :math:`z` quantization axis of the spin
operator.   

.. toctree::
   :maxdepth: 2
   :glob:

   normal/01_diag
   normal/02_hamiltonian


Superconductive mode
###############################

This set of modules implements the exact diagonalization of the single
impurity problems assuming 
:math:`\vec{Q}\equiv S_z=N_\uparrow-N_\downarrow`.

This case corresponds to the superconductive phase with :math:`s-`
wave pairing.


.. toctree::
   :maxdepth: 2
   :glob:

   superc/01_diag
   superc/02_hamiltonian



Non-SU(2) mode
###############################

This set of modules implements the exact diagonalization of the single
impurity problems assuming 
:math:`\vec{Q}\equiv N_{tot}=N_\uparrow+N_\downarrow`.

This case corresponds to the normal phase in the absence of spin
conservation, as for instance in presence of Spin-Orbit coupling.  


.. toctree::
   :maxdepth: 2
   :glob:

   nonsu2/01_diag
   nonsu2/02_hamiltonian
