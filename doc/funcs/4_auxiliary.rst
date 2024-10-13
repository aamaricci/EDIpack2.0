Auxiliary functions
===========================

.. function:: check_convergence(func,threshold,N1,N2)

   This function checks the variation of a given quantity (Weiss field, Delta, ...) against the one for the previous step. It is used to determined whether the DMFT loop has converged. If a maximum number of loops is exceeded, returns True with a warning.

   :type func: np.array(dtype=complex) 
   :param func: the quantity to be checked. It is one-dimensional, with its length being a number of frequencies
   
   :type threshold: float 
   :param threshold: the error threshold
   
   :type N1: int
   :param N1: minimum number of loops

   :type N2: int
   :param N2: maximum number of loops
   
   :return: 
    - the error
    - a boolean signifying convergence
   :rtype: float, bool


.. function:: get_bath_type(func,threshold,N1,N2)

   This function returns an integer number related to the value of :code:`BATH_TYPE` in the input file
    - :code:`1` for **normal** bath
    - :code:`2` for **hybrid** bath
    - :code:`3` for **replica** bath
    - :code:`4` for **general** bath
   
   :return: the integer index
   :rtype: int
   

.. function:: get_ed_mode(func,threshold,N1,N2)

   This function returns an integer number related to the value of :code:`ED_MODE` in the input file
    - :code:`1` for **normal** mode
    - :code:`2` for **superc** mode
    - :code:`3` for **nonsu2** mode
   
   :return: the integer index
   :rtype: int


.. function:: search_variable(var,ntmp,converged)

   This function checks the value of the read density :code:`ntmp` against the desired value :code:`ed.nread` (if different from zero) and adjusts :code:`var` accordingly (in a monotonous way).
   
   :type var: float
   :param var: the variable to be adjusted (usually :code:`ed.xmu`)

   :type ntmp: float
   :param ntmp: the density value at the given iteration
   
   :type converged: bool
   :param converged: whether the DMFT loop has achieved a sufficiently small error independently on the density
   
   :return: 
    - the new value of :code:`var`
    - a boolean signifying convergence
   :rtype: float, bool 


