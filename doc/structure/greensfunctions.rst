Normal mode
###########################

This set of modules implements the calculations of impurity dynamical
response functions, e.g. the Green's functions and different
susceptibilities,  assuming :math:`\vec{Q}=\left[\vec{N}_\uparrow,\vec{N}_\downarrow \right]`.
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

   normal/03_greensfunctions
   normal/04_chifunctions


Superconductive mode
###############################

This set of modules implements the calculations of impurity dynamical
response functions, e.g. the Green's functions and different
susceptibilities,  assuming 
:math:`\vec{Q}\equiv S_z=N_\uparrow-N_\downarrow`.

This case corresponds to the superconductive phase with :math:`s-`
wave pairing.


.. toctree::
   :maxdepth: 2
   :glob:

   superc/03_greensfunctions
   superc/04_chifunctions



Non-SU(2) mode
###############################

This set of modules implements the calculations of impurity dynamical
response functions, e.g. the Green's functions and different
susceptibilities,  assuming 
:math:`\vec{Q}\equiv N_{tot}=N_\uparrow+N_\downarrow`.

This case corresponds to the normal phase in the absence of spin
conservation, as for instance in presence of Spin-Orbit coupling.  


.. toctree::
   :maxdepth: 2
   :glob:

   nonsu2/03_greensfunctions
   nonsu2/04_chifunctions
