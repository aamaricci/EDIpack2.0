Attractive Hubbard model: a superconductive example
#######################################################





Similarly to the previous case, the preamble of the code contains
definition of  the required local variables. We load all the required
modules, including  `DMFT_TOOLS`_ library to perform some tasks
related to DMFT implementation. We read the input file using
:f:func:`ed_read_input`. Given the simplicity of this program we do
**not** make use of MPI setup here. In the absence of MPI
initialization, the `EDIpack2.0` code will automatically fall back to
serial execution. 
   

.. code-block:: fortran

   program ed_ahm_2d
      USE EDIPACK2
      USE SCIFOR
      USE DMFT_TOOLS
      implicit none
      !> a list of the variable used in this code
      
      ...
            
      !> Parse some variables
      call parse_cmd_variable(finput,"FINPUT",default='inputAHM.conf')
      call parse_input_variable(wmixing,"wmixing",finput,default=0.5d0,comment="Mixing bath parameter")

      !> read ED input
      call ed_read_input(trim(finput))




In presence of superconductivity we have to deal with Nambu structure
of the local functions and Hamiltonians, for instance the Green's
function has the matrix form:

.. math::

   \begin{pmatrix}
   \hat{G}_{\uparrow\uparrow} & \hat{F}_{\uparrow\downarrow}\\
   \bar{\hat{F}}_{\downarrow\uparrow}  &    \bar{\hat{G}}_{\downarrow\downarrow} \\
   \end{pmatrix}

where the :math:`\hat{A}` indicates the multi-orbital character of the
function, while the :math:`\bar{A}` symbol indicates the components related to the first row ones
by particle-hole and time-reversal symmetry. The actual implementation
of such symmetries depends on the argument of the functions, i.e. Matsubara or
Real-axis. Exploiting the relation among components, we need to keep track only of
the first row elements. To this end, we consider a further leading
dimension of size 2 for all local functions and Hamiltonian.  

.. code-block:: fortran

   !> Allocate local functions. The leading dimension indicates the Nambu 11,12 components.
   !    We get the other components by symmetry. 
   allocate(Gmats(2,Nso,Nso,Lmats))
   allocate(Smats(2,Nso,Nso,Lmats))
   allocate(Weiss(2,Nso,Nso,Lmats))
   allocate(Greal(2,Nso,Nso,Lreal))
   allocate(Sreal(2,Nso,Nso,Lreal))
   


Next we construct the dispersion and the exact DOS of the 2d square
lattice using `SciFortran` functions. Note again the leading dimension
2 of the dispersion.  

.. code-block:: fortran
		
   !> Build the 2d square lattice exact DOS and dispersion for the Nambu structure (use DMFT_TOOLS)
   !    We ask for separate dispersions, or two H(k), for the two diagonal Nambu components 11,22
   allocate(Edos(2,1,Le))  ![Nambu 11,12, Number of orbitals, Number of energy levels]
   Edos(1,1,:) = linspace(-D,D,Le,mesh=de)
   Edos(2,1,:) =-linspace(-D,D,Le,mesh=de)
   allocate(Ddos(1,Le))
   do i=1,Le
     Ddos(1,i) = dens_2dsquare(Edos(1,1,i),D)*de
   enddo
   !> Get local Hamiltonian (to be used in DMFT_TOOLS)
   allocate(H0(2,1))
   H0=0d0





We start the DMFT calculation by query the bath dimension (which
include superconductive parameters) and allocating the  bath on the
user side. Finally we initialize the ED solver. 


.. code-block:: fortran

   !> Get bath dimension and allocate user bath to this size
   Nb=ed_get_bath_dimension()
   allocate(Bath(Nb))
   allocate(Bath_prev(Nb))

   !> Initialize the ED solver (bath is guessed or read from file) 
   call ed_init_solver(bath)
		 


Then we implement the DMFT loop, using a similar structure as we
discussed in the previous sections:

  #. Solve the quantum impurity problem for a given user bath :math:`\vec{b}`.
  #. Retrieve Matsubara self-energies :math:`\Sigma` and :math:`S`,
     get the local Nambu Green's function
     :math:`\hat{G}_{loc}=\sum_k [i\omega_n + \mu - \hat{H}(k) - \hat{\Sigma}]^{-1}` where the :math:`\hat{A}` indicates the multi-orbital and the Nambu structure. 
  #. Implement DMFT self-consistency to update Nambu Weiss fields:
     :math:`\hat{{\cal G}}_0 = [\hat{G}^{-1}_{loc} + \hat{\Sigma}]^{-1}`
  #. Update the user bath :math:`\vec{b}` using :math:`\chi^2`
     optimization against the obtained updated Nambu Weiss fields.
  #. Check error and restart.

     
.. code-block:: fortran

   !DMFT loop
   iloop=0;converged=.false.
   do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1

     !> Solve the impurity problem, retrieve matsubara self-energy 
     call ed_solve(bath)

     !> Retrieve impurity self-energies (normal, anomalous)
     call ed_get_sigma(Smats(1,:,:,:),axis='m',type='n')
     call ed_get_sigma(Smats(2,:,:,:),axis='m',type='a')

     !> Get Gloc (using DMFT_TOOLS)
     call get_gloc(Edos,Ddos,H0,Gmats,Smats,axis='m')

     !> Update the Weiss field: (using DMFT_TOOLS).
     call dmft_self_consistency(&
          Gmats(1,:,:,:),Gmats(2,:,:,:),&
          Smats(1,:,:,:),Smats(2,:,:,:),&
          Weiss(1,:,:,:),Weiss(2,:,:,:) )

     !> Fit the new bath, starting from the current bath + the update Weiss field
     call ed_chi2_fitgf(Weiss(1,:,:,:),Weiss(2,:,:,:),bath,ispin=1)

     
     !> Linear mixing of the bath (this can be done in alternative of mixing Weiss)
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath

     !Check convergence (using DMFT_TOOLS)
     converged = check_convergence(Weiss(1,1,1,:),dmft_error,nsuccess,nloop)
   enddo




Finally, once converged DMFT solution has been obtained we can get
real-axis Green's functions using the retrieved normal and anomalous self-energies:

.. code-block:: fortran

   !> retrieve real-axis self-energies (n, a) and build local Green's function
   call ed_get_sigma(Sreal(1,:,:,:),axis='r',type='n')
   call ed_get_sigma(Sreal(2,:,:,:),axis='r',type='a')
   call get_gloc(Edos,Ddos,H0,Greal,Sreal,axis='r')




   

.. raw:: html

   <hr>


The program to solve the main model is available here:
:download:`Attractive Hubbard Code <02_ahm.f90>`

Here is a list of bath files:

  * Bath  :math:`U=-1`:   :download:`hamiltonian.restart <U1_2d_hamiltonian.restart>`
  * Bath  :math:`U=-2`:   :download:`hamiltonian.restart <U2_2d_hamiltonian.restart>`
  * Bath  :math:`U=-3`:   :download:`hamiltonian.restart <U3_2d_hamiltonian.restart>`
  * Bath  :math:`U=-4`:   :download:`hamiltonian.restart <U4_2d_hamiltonian.restart>`  
  * Bath  :math:`U=-5`:   :download:`hamiltonian.restart <U5_2d_hamiltonian.restart>`
  * Bath  :math:`U=-6`:   :download:`hamiltonian.restart <U6_2d_hamiltonian.restart>`
  * Bath  :math:`U=-8`:   :download:`hamiltonian.restart <U8_2d_hamiltonian.restart>`


An example of the input file used in the calculations can be find here:  :download:`InputFile <inputAHM.conf>`







.. _DMFT_TOOLS: https://github.com/aamaricci/DMFTtools
