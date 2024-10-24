Single Impurity Anderson problem
#########################################

In this section we shall use the methods in `EDIpack2.0` as tools to solve a simple example of Anderson quantum impurity problem. We consider the Bethe lattice DOS :math:`\rho(x)=\frac{1}{2D}\sqrt{D^2-x^2}` and build the corresponding non-interacting  Green's function :math:`G_0(z) = \int_{\mathbb{R}}\rho(\epsilon)\left[ z -\epsilon \right]^{-1}`.
We construct a discretized bath by fitting such function on the Matsubara frequencies  :math:`G_0(i\omega_n)`. Finally we use input this bath into the :code:`ED` solver of `EDIpack2.0` in presence of local interaction on the impurity.  

.. code-block:: fortran
		
   program lancED
      USE EDIPACK2
      USE SCIFOR
      implicit none
      integer,parameter                       :: Le=5000
      real(8),parameter                       :: D=1d0
      integer                                 :: Nb,Nso
      real(8),dimension(Le)                   :: Dbands
      real(8),dimension(Le)                   :: Ebands
      real(8),allocatable                     :: Bath(:)
      complex(8),allocatable,dimension(:,:,:) :: Weiss


