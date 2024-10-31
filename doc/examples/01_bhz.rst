Interacting Bernevig-Hughes-Zhang model	
###################################################################

In this section we discuss a more thorough example, using the
`EDIpack2.0` ED algorithm as a solver for DMFT in a paradigmatic
multi-orbital model of interacting electrons.

This is a Fermi-Hubbard model describing two-orbital electrons on a
square lattice. Under suitable conditions the system describes an
interacting quantum spin Hall insulator. The model is an extension of
the Bernevig-Hughes-Zhang model supplemented with Hubbard-Kanamori
interaction.

We first introduce a suitable basis of Dirac matrices
$\Gamma_{a\alpha}=\sigma_a\otimes \tau_\alpha$, where $\sigma_a$ and
$\tau_\alpha$ are Pauli matrices, respectively, in the spin and
orbital pseudo-spin space.

In terms of these matrices the model Hamiltonian reads
:math:`H=\sum_{k}\psi_{k}^\dagger H(k)\psi_{k} + H_{\rm int}`, where the
spinor :math:`\psi_{k}=[c_{1\uparrow k}, c_{2\uparrow k},
c_{1\downarrow k}, c_{2\downarrow k} ]` collect the two orbital and
two spins second quantization annihilation operators. We have then:

.. math::

   H(k) = \left[M-2t(\cos{k_x}+\cos{k_y} \right]\Gamma_{03} +
   \lambda\sin{k_x}\Gamma_{01} +   \lambda\sin{k_y}\Gamma_{02}

where :math:`M` is the mass term, which plays the role of a crystal
field splitting among the orbitals. The presence of this term breaks
the symmetry in the orbital pseudo-spin channel.

The  interaction term describes the density-density part of the
Kanamori interaction. We neglect the pair-hopping and spin-flip purely
for numerical reasons in this context.

.. math::

   H_{\rm int} = (U-J)\frac{N(N-1)}{2} - J\left( \frac{1}{4}N^2 +
   S_z^2 - 2 T_z^2\right)
   
where :math:`N=\tfrac{1}{2}\psi_i^\dagger \Gamma_{00}\psi_i` is the
total density operator,
:math:`S_z=\tfrac{1}{2}\psi_i^\dagger \Gamma_{30}\psi_i` is the total
spin polarization operator and :math:`T_z=\tfrac{1}{2}\psi_i^\dagger
\Gamma_{03}\psi_i` is the orbital pseudo-spin polarization operator.


We are now ready to start discussing the solution of this model,
starting from the non-interacting regime :math:`U=J=0`.
In this regime the model describes a quantum spin Hall insulating
phase for :math:`M<4t` and a trivial band insulator for
:math:`M>4t`. As expected a gapless Dirac state is realized at the
transition point :math:`M=4t`. 

In presence of interaction we can solve the model using DMFT. We
review here the setup of the program using `EDIpack2.0`. 


.. code-block::
   :linenos:

   program ed_bhz
      USE EDIPACK2
      USE SCIFOR
      USE DMFT_TOOLS
      USE MPI
      implicit none

      integer :: Nso,iloop,Lk
      logical :: converged
      !Bath:
      integer                                 :: Nb
      real(8),allocatable                     :: Bath(:)
      !The local hybridization function:
      complex(8),dimension(:,:,:),allocatable :: Weiss, Weiss_
      complex(8),dimension(:,:,:),allocatable :: Smats,Sreal
      complex(8),dimension(:,:,:),allocatable :: Gmats,Greal
      !hamiltonian input:
      complex(8),dimension(:,:,:),allocatable :: Hk
      complex(8),dimension(:,:),allocatable   :: Hloc
      complex(8),dimension(:,:),allocatable   :: sigmaBHZ
      real(8),dimension(:,:),allocatable      :: Zbhz
      !variables for the model:
      integer                                 :: Nx,Nkpath
      real(8)                                 :: mh,lambda,wmixing,z2
      character(len=16)                       :: finput
      !
      complex(8),dimension(4,4)               :: Gamma1,Gamma2,Gamma5,GammaN
      !MPI Vars:
      integer                                 :: comm,rank
      logical                                 :: master,fhtop
      
      
      call init_MPI()
      comm = MPI_COMM_WORLD
      call StartMsg_MPI(comm)
      rank = get_Rank_MPI(comm)
      master = get_Master_MPI(comm)
      
      
      !Parse additional variables && read Input && read H(k)^4x4
      call parse_cmd_variable(finput,"FINPUT",default='inputED_BHZ.in')  
      call parse_input_variable(Nx,"NX",finput,default=100)
      call parse_input_variable(nkpath,"NKPATH",finput,default=500)
      call parse_input_variable(mh,"MH",finput,default=0.d0)
      call parse_input_variable(lambda,"LAMBDA",finput,default=0.d0)
      call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
      !
      call ed_read_input(trim(finput))
            

As now clear from reading the quickstart guide, the preamble define
all the required local variables of the program. Note that here we
also load MPI module and :code:`DMFT_TOOLS` library to perform some
routine tasks required in DMFT calculations such as getting local
Greens's functions or performing self-consistency conditions for
matrix functions.
We init the MPI universe using `SciFortran` MPI interface in
:f:mod:`sf_mpi`. Then we read local variables and the input file. 
