Main 
=======================

The module :f:mod:`ED_MAIN` contains the functions that initializes,
launch and finalize the `EDIpack2.0` solver for the quantum impurity
problem. 

The initialization, :f:func:`ed_init` setups and allocates all the
internal variables and memory used in the code,  which remain
available to the user until :f:func:`ed_finalize` is called.  
The main function is :f:func:`ed_solve:` which aim to diagonalize the
impurity problem, evaluate the dynamical response functions and local
observables making them available to user through input/output
procedures. 

.. f:automodule::   ed_main

