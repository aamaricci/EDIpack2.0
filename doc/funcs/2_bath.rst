.. _bath:

Bath Manipulation
===========================

These functions manipulate the user-accessible bath array

.. function:: bath_inspect(bath=None,e=None,v=None,d=None,u=None)

   This function translates between the user-accessible continuous bath array and the bath components (energy level, hybridization and so on). It functions in both ways, given the array returns the components and vice-versa. It autonomously determines the type of bath and ED mode.

   :type bath: np.array(dtype=float)
   :param bath: The user-accessible bath array
   
   :type e: np.array(dtype=float)
   :param e: an array for the bath levels (:code:`ED_MODE = NORMAL, NONSU2, SUPERC`) It has dimension :code:`[ed.Nspin, ed.Norb, ed.Nbath]` for :code:`NORMAL` bath, :code:`[ed.Nspin, ed.Nbath]` for :code:`HYBRID` bath 
   
   :type v: np.array(dtype=float)
   :param v: an array for the bath hybridizations (:code:`ED_MODE = NORMAL, NONSU2, SUPERC`) It has dimension :code:`[ed.Nspin, ed.Norb, ed.Nbath]` for :code:`NORMAL` and :code:`HYBRID` bath
   
   :type d: np.array(dtype=float)
   :param d: an array for the bath anomalous enery levels(:code:`ED_MODE = SUPERC`) It has dimension :code:`[ed.Nspin, ed.Norb, ed.Nbath]` for :code:`NORMAL` bath, :code:`[ed.Nspin, ed.Nbath]` for :code:`HYBRID` bath
   
   :type u: np.array(dtype=float)
   :param u: an array for the bath spin off-diagonal hybridization (:code:`ED_MODE = NONSU2`). It has dimension :code:`[ed.Nspin, ed.Norb, ed.Nbath]` for :code:`NORMAL` and :code:`HYBRID` bath

   :raise ValueError: if both :code:`bath` and some among :code:`e,u,v,d` are provided, or the shapes are inconsistent

   :return: 
     - if :code:`bath` is provided, returns :code:`e,v`, :code:`e,d,v` or :code:`e,v,u` depending on :code:`ED_MODE`
     - if :code:`e,v`, :code:`e,d,v` or :code:`e,v,u` depending on :code:`ED_MODE` are provided, returns :code:`bath` 
   :rtype: np.array(dtype=float) 





.. function:: break_symmetry_bath(bath, field, sign, save=True)

   This function breaks the spin symmetry of the bath, useful for magnetic calculations to incite symmetry breaking. Not compatible with :code:`REPLICA` or :code:`GENERAL` bath types.

   :type bath: np.array(dtype=float)
   :param bath: The user-accessible bath array
   
   :type field: float
   :param field: the magnitude of the symmetry-breaking shift
   
   :type sign: float
   :param sign: the sign of the symmetry-breaking shift
   
   :type save: bool
   :param save: whether to save the symmetry-broken bath for reading
   
   :return: the modified bath array
   :rtype: np.array(dtype=float) 


.. function:: chi2_fitgf(*args,ispin=0,iorb=None,fmpi=True)

   This function fits the Weiss field or Hybridization function (delta) with a discretized version. The fit parameters are the bath parameters contained in the user-accessible array. Depending on the type of system we are considering (normal, superconductive, non-SU(2)) a different set of inputs has to be passed. The specifics of the numerical fitting routines are controlled in the input file.
   
    
   :type args: [np.array(dtype=complex,np.array(dtype=complex),np.array(dtype=float)] **or** [np.array(dtype=complex,np.array(dtype=float)]
   :param args: The positional arguments are the function(s) to fit and the bath array. 
   
    If the system is not superconductive (:code:`ED_MODE=NORMAL` or :code:`ED_MODE=NONSU2`) the argumens are
   
    * :code:`g`: the function to fit
    * :code:`bath`: the bath
   
    If the system is superconductive (:code:`ED_MODE=SUPERC`) the arguments are

    * :code:`g`: the normal function to fit
    * :code:`f`: the anomalous function to fit
    * :code:`bath`: the bath 
   
    The dimensions of the previous arrays can vary:
   
    The dimension of :code:`bath` can be
 
    * :code:`Nb`: if single-impurity, the output of :func:`get_bath_gimension`
    * :code:`[Nlat,Nb]`: if real-space DMFT
   
    Accordingly, the dimension of g (and f) can be:
   
    * :code:`3`: in the single-impurity case,  an array of the shape :code:`[ed.Nspin*ed.Norb, ed.Nspin*ed.Norb, Lfit]`. 
    * :code:`3`: in the real-space DMFT case, an array of the shape :code:`[Nlat*ed.Nspin*ed.Norb, Nlat*ed.Nspin*ed.Norb, Lfit]`
    * :code:`4`: in the real-space DMFT case, an array of the shape :code:`[Nlat, ed.Nspin*ed.Norb, ed.Nspin*ed.Norb, Lfit]`
    * :code:`5`: in the single-impurity case, an array of the shape :code:`[ed.Nspin, ed.Nspin, ed.Norb, ed.Norb, Lfit]`
    * :code:`6`: in the real-space DMFT case, an array of the shape :code:`[Nlat, ed.Nspin, ed.Nspin, ed.Norb, ed.Norb, Lfit]`
   
    where :code:`Lfit` is a given number of frequencies.

   
 
   :type ispin: int 
   :param ispin: spin species to be fitted. For the normal and superconducting cases, if :code:`ed.Nspin=2`, the fitting function needs to be called twice. Only the corresponding elements of :code:`bath` will be updated each time. For the non-SU(2) case, this argument is irrelevant, since all the elements of the Weiss/Delta function need to be fitted.
    
   :type iorb: int 
   :param iorb: the orbital to be fitted. If omitted, all orbitals will be fitted
   
   :type fmpi: bool 
   :param fmpi: flag to automatically do and broadcast the fit over MPI, if defined

   :raise ValueError: if the shapes of the positional arguments are incompatible
   :raise ValueError: if a number of positional arguments different from 2 or 3 are passed   
     
   :return: An array of floats that contains the bath parameters for the impurity problem. This is a required input of :func:`solve` and :func:`chi2_fitgf`. Its elements are ordered differently depending on the bath geometry. They are (de)compactified for user interaction via :func:`bath_packaging`. Specific symmetrization operations are implemented and listed in the :ref:`bath` section.
   :rtype: np.array(dtype=float) 


.. function:: get_bath_dimension()

   This function returns the correct dimension for the bath to be allocated (for each impurity) given the parameters of the system.
   
   :return: a number which is the dimension of the bath array for each impurity.
   :rtype: int  
    
   
.. function:: orb_equality_bath(bath, indx, save=True)

   This function sets every orbital component to be equal to the one of orbital :code:`indx`. Not compatible with :code:`REPLICA` or :code:`GENERAL` bath types.

   :type bath: np.array(dtype=float)
   :param bath: The user-accessible bath array
   
   :type iorb: int 
   :param iorb: the orbital index to which every other will be set as equal
      
   :type save: bool
   :param save: whether to save the symmetry-broken bath for reading
   
   :raise ValueError: if the orbital index is out of bounds
   
   :return: the modified bath array
   :rtype: np.array(dtype=float) 
   
   

.. function:: orb_symmetrize_bath(bath, orb1, orb2, save=True)

   This function enforces equality of the different-orbital components of the bath array. Not compatible with :code:`REPLICA` or :code:`GENERAL` bath types.

   :type bath: np.array(dtype=float)
   :param bath: The user-accessible bath array
   
   :type orb1: int
   :param orb1: first orbital index
   
   :type orb2: int
   :param orb2: second orbital index
      
   :type save: bool
   :param save: whether to save the symmetry-broken bath for reading
   
   :return: the modified bath array
   :rtype: np.array(dtype=float) 
 
 
 

.. function:: ph_symmetrize_bath(bath, save=True)

   This function enforces particle-hole symmetry of the bath hybridization function. Not compatible with :code:`REPLICA` or :code:`GENERAL` bath types.

   :type bath: np.array(dtype=float)
   :param bath: The user-accessible bath array
      
   :type save: bool
   :param save: whether to save the symmetry-broken bath for reading
   
   :return: the modified bath array
   :rtype: np.array(dtype=float) 
 

.. function:: save_array_as_bath(bath)

   This function takes the user-accessible array and saves it in the correct format for every bath type in the file :code:`hamiltonian.restart`

   :type bath: np.array(dtype=float)
   :param bath: The user-accessible bath array
   
   :return: Nothing
   :rtype: None

   
.. function:: set_hgeneral(hvec,lambdavec)

   This function is specific to :code:`BATH_TYPE=GENERAL`. It sets the basis of matrices and scalar parameters that, upon linear combination, make up the bath replica. The difference between the bath types :code:`REPLICA` and :code:`GENERAL` is only in the hybridization, which is independently set, so the input of this function is the same as that of :func:`set_hreplica`.
    
   :type hvec: np.array(dtype=complex)
   :param hvec: array of bath matrices. They decompose the nonzero part of the replica in a set. Each element of the set correspond to a variational parameter. That way the bath replica matrix is updated while preserving symmetries of the user's choosing. The array can have the following shapes:

    * :code:`[(Nnambu)*ed.Nspin*ed.Norb, (Nnambu)*ed.Nspin*ed.Norb, Nsym]`: 3-dimensional, where Nnambu refers to the superconducting case and Nsym is the number of matrices that make up the linear combination 
    * :code:`[(Nnambu)*ed.Nspin*, (Nnambu)*ed.Nspin, ed.Norb, ed.Norb, Nsym]`:5-dimensional, where Nnambu refers to the superconducting case and Nsym is the number of matrices that make up the linear combination 
    
   :type lambdavec: np.array(dtype=float) 
   :param iorb: the array of coefficients of the linear combination. This, along with the hybridizations V, are the fitting parameters of the bath. The array has the following shape
    * :code:`[ed.Nbath, Nsym]`: for single-impurity DMFT, 2-dimensional, where Nsym is the number of matrices that make up the linear combination 
    * :code:`[Nlat, ed.Nbath, Nsym]`: for real-space DMFT, 3-dimensional, where Nlat is the number of inequivalent impurity sites and Nsym is the number of matrices that make up the linear combination 

   :raise ValueError: if the shapes of the arrays are inconsistent
     
   :return: Nothing
   :rtype: None
   

.. function:: set_hreplica(hvec,lambdavec)

   This function is specific to :code:`BATH_TYPE=REPLICA`. It sets the basis of matrices and scalar parameters that, upon linear combination, make up the bath replica.
    
   :type hvec: np.array(dtype=complex)
   :param hvec: array of bath matrices. They decompose the nonzero part of the replica in a set. Each element of the set correspond to a variational parameter. That way the bath replica matrix is updated while preserving symmetries of the user's choosing. The array can have the following shapes:

    * :code:`[(Nnambu)*ed.Nspin*ed.Norb, (Nnambu)*ed.Nspin*ed.Norb, Nsym]`: 3-dimensional, where Nnambu refers to the superconducting case and Nsym is the number of matrices that make up the linear combination 
    * :code:`[(Nnambu)*ed.Nspin*, (Nnambu)*ed.Nspin, ed.Norb, ed.Norb, Nsym]`:5-dimensional, where Nnambu refers to the superconducting case and Nsym is the number of matrices that make up the linear combination 
    
   :type lambdavec: np.array(dtype=float) 
   :param iorb: the array of coefficients of the linear combination. This, along with the hybridizations V, are the fitting parameters of the bath. The array has the following shape
    * :code:`[ed.Nbath, Nsym]`: for single-impurity DMFT, 2-dimensional, where Nsym is the number of matrices that make up the linear combination 
    * :code:`[Nlat, ed.Nbath, Nsym]`: for real-space DMFT, 3-dimensional, where Nlat is the number of inequivalent impurity sites and Nsym is the number of matrices that make up the linear combination 

   :raise ValueError: if the shapes of the arrays are inconsistent
     
   :return: Nothing
   :rtype: None

   
.. function:: spin_symmetrize_bath(bath, save=True)

   This function enforces equality of the opposite-spin components of the bath array. Not compatible with :code:`REPLICA` or :code:`GENERAL` bath types.

   :type bath: np.array(dtype=float)
   :param bath: The user-accessible bath array
      
   :type save: bool
   :param save: whether to save the symmetry-broken bath for reading
   
   :return: the modified bath array
   :rtype: np.array(dtype=float) 
