Bath Replica routines
============================

..
 .. raw:: html
    :file:  ../graphs/bath_dmft/03_ed_bath_replica.html

 |

In these module we implement the functions to set the matrix basis
:math:`\{ \hat{O}_i \}_{i=1,\dots,N_{sym}}` and the initial
variational parameters :math:`\vec{\lambda}` used to decompose each
local bath hamiltonian for the  :f:var:`replica` and :f:var:`general`
bath types.  

.. f:automodule::   ed_bath_replica


.. |Nbath| replace:: :f:var:`nbath`
.. |Nsym| replace:: :f:var:`nsym`
.. |Nnambu| replace:: :f:var:`nnambu`
.. |Nlat| replace:: :f:var:`nlat`
.. |Nspin| replace:: :f:var:`nspin`
.. |Norb| replace:: :f:var:`norb`
.. |Nso| replace:: :f:var:`nspin` . :f:var:`norb`
.. |Nns| replace:: :f:var:`nnambu` . :f:var:`nspin`
.. |Nlso| replace:: :f:var:`nlat`. :f:var:`nspin` . :f:var:`norb`
.. |Nnso| replace:: :f:var:`nnambu` . :f:var:`nspin`. :f:var:`norb`
