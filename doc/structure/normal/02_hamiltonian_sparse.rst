Stored Hamiltonian :math:`H\times\vec{v}`  
==============================================

The module :f:mod:`ed_hamiltonian_normal_stored_hxv` constructs the
terms of the sector Hamiltonian storing them into different
:f:var:`sparse_matrix` instances.
The main output of this module are the matrix vector products performed
using the stored sparse matrices :math:`\vec{w} = H\times \vec{v}`.

Different constructions and matrix-vector products are implemented
according to the serial or paralle mode and the value of the variable
:f:var:`ed_total_ud` which distinguish between conservation of total
or orbital resolved spin occupations. 


.. f:automodule::  ed_hamiltonian_normal_stored_hxv

.. |Nbath| replace:: :f:var:`nbath`
.. |Nlat| replace:: :f:var:`nlat`
.. |Nspin| replace:: :f:var:`nspin`
.. |Norb| replace:: :f:var:`norb`
.. |Nso| replace:: :f:var:`nspin` . :f:var:`norb`
.. |Nlso| replace:: :f:var:`nlat`. :f:var:`nspin` . :f:var:`norb`
