Global variables
=================

These are global variables of the edipy2 module. They form a subset of the global variables of the EDIpack2 library. Along with all the other global variables, they can be set in the input file, and are read when calling the :func:`read_input` function.

The exposed global variables can be accessed as methods of the :class:`global_env` class.

.. code-block:: python

    import numpy as np
    from edipy2 import global_env as ed
   
    ed.Nspin = 1            # set a global variable
    mylocalvar = ed.Nspin   # assing to a local variable (the value of mylocalvar will not change if ed.Nspin changes)
    print(ed.Nspin)         # all functions can have global variables as arguments
    np.arange(ed.Nspin)


.. data:: beta

   Value of the inverse temperature, at T=0 is used as a IR cut-off
   
   :type: float
   :default: 1000.0

.. data:: Jh

   Value of the Hund's coupling
   
   :type: float
   :default: 0.0
   
.. data:: dmft_error

   Error threshold for DMFT convergence
   
   :type: float
   :default: 1e-05
   
.. data:: ed_total_ud

   Flag to select which type of quantum numbers have to be considered: T (default) total Nup-Ndw, F orbital based Nup-Ndw
   
   :type: bool
   :default: True
   
.. data:: ed_twin

   Flag to reduce (T) or not (F,default) the number of visited sector using twin symmetry
   
   :type: bool
   :default: False
   
.. data:: eps

   Broadening on the float-axis
   
   :type: float
   :default: 1e-02

.. data:: Jx

   Value of the spin exchange coupling
   
   :type: float
   :default: 0.0

.. data:: Jp

   Value of the pair hopping coupling
   
   :type: float
   :default: 0.0
   
.. data:: Lfloat

   Number of frequencies, float frequency axis
   
   :type: int
   :default: 5000

.. data:: Lmats

   Number of frequencies, Matsubara axis
   
   :type: int
   :default: 4096
  
.. data:: LOGfile

   Log unit
   
   :type: int
   :default: 6
   
.. data:: Lpos

   Number of points for the lattice PDF
   
   :type: int
   :default: 100
   
.. data:: Ltau

   Number of imaginary time points
   
   :type: int
   :default: 1024

.. data:: Nbath

   Number of bath levels. See the specifics of the bath geometries
   
   :type: int
   :default: 6
   
.. data:: Nloop

   Maximum number of DMFT loops
   
   :type: int
   :default: 100

.. data:: Norb

   Number of correlated orbitals. Maximum 5 orbitals are supported
   
   :type: int
   :default: 1

.. data:: Nph

   Max number of phonons allowed (cut off)
   
   :type: int
   :default: 0
   
.. data:: nread

   Value of the target density for fixed density calculations.
   If valued 0, it is discarded.
   
   :type: float
   :default: 0.0

.. data:: Nspin

   Number of explicitly defined spin degrees of freedom. If Nspin=1, the two spin block of the Hamiltonian, Green's function, Self-Energy and so on are assumed equal.
   If Nspin=2 they may differ (e.g. for non-SU(2) or magnetic systems).
   The superconductive variant of the code requires Nspin=1
   
   :type: int
   :default: 1
   
.. data:: Nsuccess

   Number of successive iterations below threshold for convergence
   
   :type: int
   :default: 1
   
.. data:: sb_field

   Value of a symmetry breaking field for magnetic solutions
   
   :type: float
   :default: 0.1


.. data:: Uloc

   Values of the local interaction per orbital (max 5). If less values are provided, the array is filled in increasing order
   
   :type: float
   :default: [2.0, 0.0, 0.0, 0.0, 0.0]
   
.. data:: Ust

   Value of the inter-orbital interaction term
   
   :type: float
   :default: 0.0
   
.. data:: wini

   Value of the smallest float-axis frequency
   
   :type: float
   :default: -5.0
   
.. data:: wfin

   Value of the largest float-axis frequency
   
   :type: float
   :default: -5.0
   
.. data:: xmin

   Value for the smallest position for the lattice PDF
   
   :type: float
   :default: -3.0

.. data:: xmax

   Value for the largest position for the lattice PDF
   
   :type: float
   :default: 3.0

   
.. data:: xmu

   Value of the chemical potential. If HFMODE = T, xmu=0 satisfied the half-filling condition
   
   :type: float
   :default: 0.0

