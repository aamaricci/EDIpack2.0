Install
#####################

EDIpack2.0 is available in the form of a static Fortran library
`libedipack2.a` and the related Fortran module `EDIPACK2`.
Installation is available using CMake. In the current release the
library enables standard Fortran API and Python API `EDIpy2`. 


Building
======================

We assume that `SciFortran` and `MPI` have been correctly installed
and are available in the system. See related documentation. Note that
the installation of `EDIpack2.0` closely follows the `SciFortran`
template.


Clone the repo:

.. code-block:: bash
		
   git clone https://github.com/aamaricci/EDIpack2.0 EDIpack2



Optionally define the fortran compiler:

.. code-block:: bash
		
   export FC=mpif90/gfortran/ifort


From the repository directory (`cd EDIpack2`) make a standard
out-of-source CMake compilation:

**GNU Make**

Using GNU `make` is the default CMake workflow, with widest version
support (CMake > 3.0). Note that parallel `make` execution is tested
and working.

.. code-block:: bash
		
   mkdir build 
   cd build  
   cmake .. 
   make -j



**Ninja**

Using `ninja` if a fortran-capable version of `ninja
<https://ninja-build.org>`_ is available in your system (and CMake can
take advantage of it), you can use it to build the library at lightning, multi-threaded, speed. 

.. code-block:: bash
		
   mkdir build    
   cd build  
   cmake -GNinja ..  
   ninja

The `CMake` compilation can be customized using the following
additional variables (default values between `< >`, optional in `[ ]`):   

.. code-block:: bash
		
   -DPREFIX=prefix directory  #<~/opt/EDIpack2/VERSION/PLATFORM/[GIT_BRANCH]>` 

   -DUSE_MPI=<yes>/no    

   -DVERBOSE=yes/<no>

  -DBUILD_TYPE = <RELEASE>/TESTING/DEBUG/AGGRESSIVE  #compilation options. TESTING:mild or no optimization, DEBUG:relevant debugging options, AGGRESSIVE: all debug options of (might not compile on  some systems) 



Install
======================

System-wide installation is completed after the build step using either:

.. code-block:: bash

   make install

or

.. code-block:: bash
		
   ninja install

  
Please follow the instructions on the screen to complete installation on your environment.  
The library can be loaded using one of the following, automatically generated, files :  

*  A generated `environment module`_ , installed to`~/.modules.d/EDIpack2/<PLAT>`
  
* A generated `bash` script at `<PREFIX>/bin/configvars.sh`, to be sourced for permanent loading.

*  A generated `pkg-config`_ file to, installed to `~/.pkg-config.d/EDIpack2.pc`  

.. _environment module: https://github.com/cea-hpc/modules
.. _pkg-config: https://github.com/freedesktop/pkg-config


Uninstall
===================

Although CMake does not officially provide uninstall procedures in the
generated Make/Ninja files. Hence SciFortran supplies a homebrew
method to remove the generated files by calling (from the relevant
build folder):

.. code-block:: bash
		
   make uninstall

or

.. code-block:: bash
		
   ninja uninstall




Install Python API
======================
The `edipy2` python module is installable from this folder via:

.. code-block:: bash
		
    pip install .

on some systems such as Debian >= 11 and Mac Os

.. note:: 
   If a virtual environment is not in use, the flag
   `--break-system-packages` has to be set. This creates no issue
   since no distro is packaging this library.
   

To remove the module, run:

.. code-block:: bash
		
   pip uninstall -y edipy2

with same caveat for the `--break-system-packages` flag.


.. tip::

   See `EDIpy2` documentation for more details on installing the `python` API. 
