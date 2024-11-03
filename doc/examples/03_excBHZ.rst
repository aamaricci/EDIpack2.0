In-plane excitons in quantum spin Hall insulators
#####################################################

In this finale example we discuss an interesting symmetry breaking
case in the same QSHI model introduced in :doc:`01_bhz`.
Specifically, we consider the case of excitonic condensation described
by the order  parameters :math:`\vec{E}=[E_0,E_1,E_2,E_3]`  which, in terms of the Gamma
matrices, read :math:`E_a = \langle \psi_i^\dagger \Gamma_{a1} \psi_i
\rangle`.

The analysis of the impurity susceptibilities in the previous section,
suggest a possible instability towards the in-plane triplet exciton
state :math:`E_1` or :math:`E_2`. Interestingly, the onset of this
state breaks several symmetries, e.g. time-reversal and spin SU(2).
As thus, this is a physical case to investigate the :f:var:`ed_mode` =
**nonsu2** mode in `EDIpack2.0`. 


The considerations about the model Hamiltonian and the non-interacting
solution remain identical to the previous case :doc:`01_bhz`. Here we
discuss how to change the program to tackle the specific issue of
in-plane exciton condensation in QSHI.

The general structure of the code is unchanged but for one important
part: the bath construction. We get:

.. code-block:: fortran
   :linenos:

      
   !> Get local Hamiltonian summing over k (one can do better)
   allocate(Hloc(Nso,Nso))
   Hloc = sum(Hk,dim=3)/Lk
   where(abs(dreal(Hloc))<1d-6)Hloc=zero
   !> Set H_{loc} in EDIpack2
   call ed_set_hloc(Hloc)
   !> Get bath dimension and allocate user bath to this size
   ! ~removed~[Nb=ed_get_bath_dimension()]~
   !> Setup the replica bath basis for the case E0EzEx(singlet,tripletZ,tripletX)
   allocate(lambdasym_vector(Nbath,4))
   allocate(Hsym_basis(Nso,Nso,4))
   Hsym_basis(:,:,1)=Gamma5  ;lambdasym_vector(:,1)= Mh
   Hsym_basis(:,:,2)=GammaE0 ;lambdasym_vector(:,2)= sb_field
   Hsym_basis(:,:,3)=GammaEz ;lambdasym_vector(:,3)= sb_field
   Hsym_basis(:,:,4)=GammaEx ;lambdasym_vector(:,4)=-sb_field
   !> Set the replica bath
   call ed_set_Hreplica(Hsym_basis,lambdasym_vector)
   !> Get bath dimension and allocate user bath to this size
   Nb=ed_get_bath_dimension(4)   !(Hsym_basis)
   allocate(Bath(Nb))
   !
   !> Initialize the ED solver (bath is guessed or read from file) 
   call ed_init_solver(bath)


We first generate a basis of 4 matrices
:math:`O_i=[\Gamma_{03},\Gamma_{01},\Gamma_{31},\Gamma_{11}]` and a
set of parameters :math:`\vec{\lambda}^p=[ \lambda^p_1,\lambda^p_2,\lambda^p_3,\lambda^p_4]`
which will be used to parametrize any bath Hamiltonian as :math:`h^p
=\sum_i \lambda_i^p O_i`.
The rest of the implementation is unaltered, except for a couple of
printing flags. 


.. raw:: html

   <hr>


We can now discuss some results obtained with this code for the
exciton condensation. More results can be found in `PhysRevB.107.115117`_. 

.. _PhysRevB.107.115117: https://journals.aps.org/prb/abstract/10.1103/PhysRevB.107.115117







.. raw:: html

   <hr>


The program to solve the main model can be found here:
:download:`Exciton BHZ Code <03_excBHZ.f90>`
