Stored Hamiltonian :math:`H\times\vec{v}`  
==============================================

The module :f:mod:`ed_hamiltonian_superc_stored_hxv` constructs the
terms of the sector Hamiltonian storing them into different
:f:var:`sparse_matrix` instances: :f:var:`sph0` for the electronic
part and three others for the phononic and electron-phonon terms. 

The main output of this module are the matrix vector products performed
using the stored sparse matrices :math:`\vec{w} = H\times \vec{v}`.



.. f:automodule::  ed_hamiltonian_superc_stored_hxv



.. |Nbath| replace:: :f:var:`nbath`
.. |Nlat| replace:: :f:var:`nlat`
.. |Nspin| replace:: :f:var:`nspin`
.. |Norb| replace:: :f:var:`norb`
.. |Nso| replace:: :f:var:`nspin` . :f:var:`norb`
.. |Nlso| replace:: :f:var:`nlat`. :f:var:`nspin` . :f:var:`norb`
