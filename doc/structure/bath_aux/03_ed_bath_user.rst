User Bath symmetry operations
=================================

..
 .. raw:: html
    :file:  ../graphs/bath_dmft/02_ed_bath_user.html

 |


In this module we implement few functions to let the user implement
common symmetry operation on the bath array. In some cases using such
operations can help the convergence of DMFT calculations.


.. f:automodule::   ed_bath_user
   :members: break_symmetry_bath,spin_symmetrize_bath,orb_symmetrize_bath,orb_equality_bath,ph_symmetrize_bath,ph_trans_bath,enforce_normal_bath,impose_equal_lambda,save_array_as_bath
