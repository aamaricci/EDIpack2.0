Hamiltonian Direct :math:`H\times\vec{v}`
==============================================

The module :f:mod:`ed_hamiltonian_normal_direct_hxv` constructs and
applied on-the-fly each term of the sector Hamiltonian to the input
vector as :math:`\vec{w} = H\times \vec{v}` in a Arpack/Lanczos framework.

Different functions are implemented according to the value of the variable
:f:var:`ed_total_ud`, which distinguish between conservation of total
or orbital resolved spin occupations and the serial or parallel case.


.. f:automodule::  ed_hamiltonian_normal_direct_hxv



.. |Nbath| replace:: :f:var:`nbath`
.. |Nlat| replace:: :f:var:`nlat`
.. |Nspin| replace:: :f:var:`nspin`
.. |Norb| replace:: :f:var:`norb`
.. |Nso| replace:: :f:var:`nspin` . :f:var:`norb`
.. |Nlso| replace:: :f:var:`nlat`. :f:var:`nspin` . :f:var:`norb`
