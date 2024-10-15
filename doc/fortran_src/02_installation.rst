Installation
#####################

EDIpack2.0 is available in the form of a static Fortran library
`libedipack2.a` and the related Fortran module `EDIPACK2`.
Installation is available using CMake. In the current release the
library enables standard Fortran API and Python API `EDIpy2`. 


Installing EDIpack2.0
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

GNU Make
------------
Using GNU `make` is the default CMake workflow, with widest version
support (CMake > 3.0). Note that parallel `make` execution is tested
and working.

.. code-block:: bash
		
   mkdir build 
   cd build  
   cmake .. 
   make -j



Ninja
------------
Using `ninja` if a fortran-capable version of `ninja <https://ninja-build.org>`_ is available in your system (and CMake can take advantage of it), you can use it to build the library at lightning, multi-threaded, speed. 

.. code-block:: bash
		
   mkdir build    
   cd build  
   cmake -GNinja ..  
   ninja




The `CMake` compilation can be customized using the following
additional variables (default values between `< >`):   

* `-DPREFIX=prefix directory <~/opt/EDIpack2/VERSION/PLATFORM/[GIT_BRANCH]>` 

* `-DUSE_MPI=<yes>/no`  

* `-DVERBOSE=yes/<no>`  

* `-DBUILD_TYPE=<RELEASE>/TESTING/DEBUG/AGGRESSIVE`  


Please follow the instructions on the screen to complete installation on your environment.  
The library can be loaded using one of the following, automatically generated, files :  

* pkg-config file in `~/.pkg-config.d/EDIpack2.pc`  
* environment module file `~/.modules.d/EDIpack2/<PLAT>`  
* homebrew `bash` script `<PREFIX>/bin/configvars.sh`


Python API: EDIpy2 
======================
The `edipy2` python module is installable from this folder via:
o R
.. code-block:: bash
		
   pip install .

on some systems such as Debian >= 11 and Mac Os, if a virtual environment is not in use, the flag `--break-system-packages` has to be set. This creates no issue since no distro is packaging this library.
To remove the module, run:

.. code-block:: bash
		
   pip uninstall -y edipy2

with same caveat for the `--break-system-packages` flag.

See `EDIpy2` documentation for more information. 
