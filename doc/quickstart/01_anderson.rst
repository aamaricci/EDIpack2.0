Single Impurity Anderson problem
#########################################

In this section we use the methods in `EDIpack2.0`  to solve a simple
example of Anderson quantum impurity problem.

Looking forward for a DMFT application, here we consider the Bethe
lattice DOS :math:`\rho(x)=\frac{1}{2D}\sqrt{D^2-x^2}` and build the
corresponding non-interacting  Green's function :math:`G_0(z) =
\int_{\mathbb{R}}\rho(\epsilon)\left[ z -\epsilon \right]^{-1}`.

We construct a discretized bath by fitting such function on the
Matsubara frequencies  :math:`G_0(i\omega_n)` using the methods in
:ref:`fit`. Finally we input this bath into the :f:func:`ed_solve` solver of `EDIpack2.0` in presence of local interaction on the impurity.  


The initialization of the code is:

.. code-block:: fortran
		
   program lancED
      USE EDIPACK2
      USE SCIFOR
      implicit none
      integer,parameter                       :: Le=5000  !
      real(8),parameter                       :: D=1d0  !
      integer                                 :: Nb,Nso  !
      real(8),dimension(Le)                   :: Dbands  !
      real(8),dimension(Le)                   :: Ebands  !
      real(8),allocatable                     :: Bath(:)  !
      complex(8),allocatable,dimension(:,:,:) :: Weiss  !

      !> READ THE input using EDIpack procedure: 
      call ed_read_input('inputED.conf')

      
where we load both the `EDIpack2.0` and `SciFortran` libraries through
their main module :f:mod:`edipac2` and :f:mod:`scifor`. We also define
some local variables and proceed with reading the input file
:code:`"inputED.conf"` using the `EDIpack2.0` function
:f:func:`ed_read_input`.

Next we construct the Bethe lattice DOS and non-interacting Green's
function, using procedures available in :f:var:`SciFortran`:


.. code-block:: fortran
		
   !> Bethe Lattice linear disperion
   Ebands = linspace(-D,D,Le,mesh=de)
   !> Bethe Lattice DOS (rescaled to mesh)
   Dbands = dens_bethe(Ebands,D)


   !> Get the Bethe lattice non-interacting Matsubara GF as a guess for the bath 
   allocate(Weiss(Nso,Nso,Lmats))
   call bethe_guess_g0(Weiss(1,1,:),D,beta,hloc=0d0)



Then we initialize the solver. This step requires the user to pass
the user bath as a rank-1 double precision array of a given size. The
correct size of the bath array is evaluated internally by the
`EDIpack2.0` code through the function
:f:func:`ed_get_bath_dimension`.

.. code-block:: fortran

   !> Init solver: 
   Nb=ed_get_bath_dimension()
   allocate(bath(Nb))
   call ed_init_solver(bath)


Upon initialization the bath is guessed from a flat distribution
centered around zero and with half-width :f:var:`ed_hw_bath`. Here we
update the bath optimizing it against the non-interacting Bethe
lattice Green's function:

.. code-block:: fortran

   !> Fit the bath against G0 guess: the outcome is a bath discretizing the Bethe DOS.
   call ed_chi2_fitgf(Weiss,bath,ispin=1,iorb=1)


We are now ready to solve the quantum impurity problem for a given set
of parameters specified in the input file (see below)


.. code-block:: fortran

   !> Solve SIAM with this given bath
   call ed_solve(bath)





.. raw:: html

   <hr>


Here is an example of input file used in actual calculations:

.. code-block::

   NORB=1                                        !Number of impurity orbitals (max 5).
   NBATH=9                                       !Number of bath sites
   ULOC=1d0				       
   ED_MODE=normal                                !Flag to set ED type
   BATH_TYPE=normal                              !flag to set bath type
   BETA=1000.000000000                           !Inverse temperature, at T=0 is used as a IR cut-off.
   ED_VERBOSE=3                                  !Verbosity level: 0=almost nothing --> 5:all. Really: all
   LMATS=4096                                    !Number of Matsubara frequencies.
   LFIT=2048
   EPS=1.000000000E-02                           !Broadening on the real-axis.
   CG_FTOL=1.000000000E-10                       !Conjugate-Gradient tolerance.
   CG_NITER=2048                                 !Max. number of Conjugate-Gradient iterations.
   ED_TWIN=T
   LANC_NGFITER=500                              !Number of Lanczos iteration in GF determination. Number of momenta.
   ED_HW_BATH=1.000000000                        !half-bandwidth for the bath initialization: flat in -ed_hw_bath:ed_hw_bath


   
