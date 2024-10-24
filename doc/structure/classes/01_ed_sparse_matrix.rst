Sparse Matrices 
=======================

.. raw:: html
   :file:  ../graphs/classes/01_ed_sparse_matrix.html

|


This class defines a data structure to efficiently store sparse
matrices into dedicated CSR matrices. The class features support to
MPI parallel storage, so that each matrix is spread across the
threads.

Each instance of :f:var:`sparse_matrix_csr` corresponds to a rank-1
array of dynamically reallocated rows, :f:var:`sparse_row_csr`,
including one array for columns indices and one array for non-zero
values of the dense matrix. Access to each element is :math:`O(1)`.

.. note::
   A more modern, object oriented, class for sparse matrices is
   available in SciFortran_

.. _SciFortran: https://github.com/SciFortran/SciFortran/tree/master/src/SF_SPARSE


.. f:automodule::   ed_sparse_matrix
   :members: sparse_row_csr,sparse_matrix_csr, sp_init_matrix, sp_insert_element,sp_dump_matrix, sp_insert_element,sp_set_mpi_matrix




