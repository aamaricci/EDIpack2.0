Pair Susceptibility
============================



In :f:mod:`ed_chi_dens` we evaluate the impurity pair 
susceptibility, defined as:

.. math::

   \chi^{\Delta}_{ab}(\omega) = \langle \Delta_a(\omega) \Delta_b(\omega) \rangle = \frac{1}{\cal
   Z}\sum_m e^{-\beta E_m} \langle m | \Delta_a [\omega-H]^{-1} \Delta_b  | m \rangle

where :math:`\Delta_a = c_{a\uparrow} c_{a\downarrow}` is the fermion
singlet pair operator of the orbital :math:`a` and :math:`\omega \in {\mathbb C}`. As for the
Green's functions, the susceptibility is evaluated using the dynamical
Lanczos method: a) the partial tridiagonalization of the 
sector Hamiltonian :math:`H` with quantum numbers
:math:`\vec{Q}=[\vec{N}_\uparrow,\vec{N}_\downarrow]` on the Krylov
basis of :math:`n_a|m\rangle` is obtained; b) the resulting
tridiagonal matrix is further diagonalized to obtained excitations
amplitudes or **weights**  :math:`\langle p | \Delta_a | m \rangle` for
any state :math:`| p \rangle` in the spectrum (*without knowing the
state itself* ) and the excitations energies :math:`\delta E = E_p -
E_m` or **poles**; c) an controlled approximation to the
Kallen-Lehmann sum is constructed for  :math:`a,b=1,\dots,N_{\rm
orb}`. 




.. f:automodule::  ed_chi_pair



.. |Nbath| replace:: :f:var:`nbath`
.. |Nlat| replace:: :f:var:`nlat`
.. |Nspin| replace:: :f:var:`nspin`
.. |Norb| replace:: :f:var:`norb`
.. |Nso| replace:: :f:var:`nspin` . :f:var:`norb`
.. |Nlso| replace:: :f:var:`nlat`. :f:var:`nspin` . :f:var:`norb`
