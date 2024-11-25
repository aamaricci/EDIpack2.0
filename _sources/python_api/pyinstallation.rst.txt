EDIpy2 Install
==============

The `python` module `EDIpy2` requires:

* `SciFortran <https://github.com/aamaricci/SciFortran>`_

* `EDIpack2 <https://github.com/aamaricci/EDIpack2.0>`_

libraries to be installed beforehand, please see related documentation
to install such Fortran libraries. Once both are set up, the python module can be installed from the EDIpack2 folder via

.. code-block:: shell

   pip install . --break-system-packages
   
The latter option may not be required in all cases, but it is in recent versions of Debian and OSX. Since no edipy2 package is provided by any distro, this will not create problems. If the user is using a virtual environment, the option is not necessary.





