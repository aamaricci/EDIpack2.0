EDIpack2.0
################

A massively parallel Exact Diagonalization solver for quantum Impurity problems.
***************************************************************************************************************

*The documentation is under construction*


EDIpack2.0 is a Lanczos based method 
for the solution of generic Quantum Impurity problems,  exploiting distributed memory MPI parallelisation.
The 2.0 version extends the EDIpack_ by enabling the solution of
single-site, multi-orbital models with different conserved
quantum numbers :math:`\vec{Q}` corresponding to separate operational
modes which, in `EDIpack2.0` software, are selected by the input
variable `ed_mode=normal,superc,nonsu2` as follow: 

* :math:`\vec{Q}=[\vec{N}_\uparrow,\vec{N}_\downarrow]` for which
  either the number of total or orbital spin up and down electrons is
  conserved: **NORMAL** 

* :math:`\vec{Q}=S_z`  with conserved total magnetization:  **SUPERConducting**  

* :math:`\vec{Q}=N_{\rm tot}`  where spin degrees freedom is not fully conserved:  **NON-SU(2)**


.. note::
   The `superc` mode deals with local *s*-wave pairing although in 
   diagonal and off-diagonal orbital channels. The actual
   implementation does not support long-range magnetic ordering.
   
.. note::
   The `nonsu2` operational mode deals with any situation in which
   spin symmetry group is not fully conserved, for instance in
   presence of local Spin-Orbit Coupling :math:`\vec{L} \cdot \vec{S}`,
   in-plane magnetization :math:`\langle S_x\rangle\gt0`  or in-plane
   triplet excitonic condensation, see `PhysRevB.107.115117`_. 

.. _PhysRevB.107.115117: https://journals.aps.org/prb/abstract/10.1103/PhysRevB.107.115117


All modes include electron-phonon coupling (local or Holstein
phonons). EDIpack2.0  is designed to obtain the lowest part of the
spectrum of the problem, thus it naturally works at zero temperature
but can also be used to explore low temperature properties.  
 
The EDIpack2.0 diagonalization algorithm is based on a massively
parallel execution of matrix-vector products, required in the context
of Lanczos-Arnoldi linear procedures.  See `j.cpc.2021.108261`_  for a detailed
descriptions of these algorithms.
However, substantial modifications have been introduced in the 2.0
version to address the *Superconducting* and *non-SU(2)* channels.  
An updated manuscript will be released soon. 

.. _EDIPACK: https://github.com/aamaricci/EDIpack
.. _j.cpc.2021.108261: https://doi.org/10.1016/j.cpc.2021.108261


Authors
=================

The `EDIpack` libraries (1.0 and 2.0) have been developed as a
collective effort by different authors, each contributing to diverse
aspects of the library. The following list does not follow any
particular order:  

* `Adriano Amaricci`_ (leading author)
  
* `Lorenzo Crippa`_
  
* `Samuele Giuli`_

* `Gabriele Bellomia`_

* `Giacomo Mazza`_

* Alberto Scazzola

* Luca de Medici
  
* Massimo Capone

.. _Adriano Amaricci: https://github.com/aamaricci
.. _Lorenzo Crippa: https://github.com/lcrippa    
.. _Samuele Giuli: https://github.com/SamueleGiuli
.. _Gabriele Bellomia: https://github.com/beddalumia
.. _Giacomo Mazza: https://github.com/GiacMazza


Installation
=================

:doc:`dependencies`
     Software requirements to install `EDIpack2.0`
     
:doc:`installation`
     Build, install and configure `EDIpack2.0`
     

Quick Start 
=================

:doc:`quickstart`
     A quick start guide to  `EDIpack2.0` usage

:doc:`examples`
     Some examples illustrating the use of `EDIpack2.0` for simple test problems
     

EDIpack2.0
======================

:doc:`edipack2`
     An overview of the structure of the library and a detailed
     description of the relevant modules.


EDIpy2
======================

:doc:`edipy2`
     Installation and basic use of `EDIpy2`: the `python` API of `EDIpack2.0`

.. Hidden TOCs

.. toctree::
   :caption: Installation
   :maxdepth: 2
   :hidden:

   dependencies
   installation

.. toctree::
   :caption: Quick Start
   :maxdepth: 2
   :hidden:

   quickstart
   examples
   
.. toctree::
   :caption: EDIpack2.0
   :maxdepth: 2
   :hidden:
      
   edipack2


.. toctree::
   :caption: EDIpy2
   :maxdepth: 2
   :hidden:

   edipy2




   

