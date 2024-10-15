# EDIpack2.0: Massively parallel Exact Diagonalization for generic Quantum Impurity problems

[![TestSuite](https://img.shields.io/github/actions/workflow/status/aamaricci/EDIpack2.0/PushWorkflow.yml?label=TestSuite&logo=Fortran&style=flat-square)](https://github.com/aamaricci/EDIpack2.0/actions/workflows/PushWorkflow.yml) 

<!-- TO BE SETUP ASAP
[![Coverage]()]()
[![api docs](https://img.shields.io/static/v1?label=API&message=documentation&color=734f96&logo=read-the-docs&logoColor=white&style=flat-square)](https://qcmplab.github.io/DMFT_ED)
-->
A suitable extension of [EDIpack](https://github.com/aamaricci/EDIpack): a  Lanczos based method 
for the solution of generic Quantum Impurity problems,  exploiting distributed memory MPI parallelisation.
This updated version, aims to solve single-site, multi-orbital models, in either  *normal*, *superconducting* (s-wave) or *Spin-non-conserving* (e.g. with Spin-Orbit Coupling or in-plane magnetization) phases, including electron-phonons coupling. The code works at zero and low temperatures.   
 
See [j.cpc.2021.108261](https://doi.org/10.1016/j.cpc.2021.108261) for further information about the underlying algorithms. Yet, suitable modifications have been developed to address the Superconducting and non-SU(2) channels.  

### Documentation
The documentation for the EDIpack2 library and its Python API (EDIpy2) is available at [https://aamaricci.github.io/EDIpack2.0/](https://aamaricci.github.io/EDIpack2.0/)  

The documentation is under construction:  
- [x] Setup the general structure  
- [x] Merge with the EDIpy2 documentation (from L.Crippa)  
- [ ] Draft the main index of the documentation  


### Dependencies

The code is based on:  

* SciFortran [https://github.com/aamaricci/SciFortran](https://github.com/aamaricci/SciFortran)  

* MPI 

  


### Installation

Installation is available using CMake. In the current version API are only provided in Fortran and Python.  
The software gives acces to the static library `libedipack2.a` and the related modules `EDIPACK2`

Clone the repo:

`git clone https://github.com/aamaricci/EDIpack2.0 EDIpack2`

And from the repository directory (`cd EDIpack2`) make a standard out-of-source CMake compilation:

`mkdir build`  
`cd build`  
`cmake ..`      
`make`     
`make install`   

Please follow the instructions on the screen to complete installation on your environment.  
The library can be loaded using one of the following, automatically generated, files :  

* pkg-config file in `~/.pkg-config.d/EDIpack2.pc`  
* environment module file `~/.modules.d/EDIpack2/<PLAT>`  
* homebrew `bash` script `<PREFIX>/bin/configvars.sh`


The `CMake` compilation can be controlled using the following additional variables, default values between `< >`:   

* `-DPREFIX=prefix directory <~/opt/EDIpack2/VERSION/PLAT/[GIT_BRANCH]>` 

* `-DUSE_MPI=<yes>/no`  

* `-DVERBOSE=yes/<no> `  

* `-DBUILD_TYPE=<RELEASE>/TESTING/DEBUG`  


### Python module

The `edipy2` python module is installable from this folder via

`pip install .`

on some systems such as Debian >= 11 and Mac Os, if a virtual environment is not in use, the flag `--break-system-packages` has to be set. This creates no issue since no distro is packaging this library.
To remove the module, run

`pip uninstall -y edipy2`

with same caveat for the `--break-system-packages` flag.


--

### CONTACT

If you encounter bugs or difficulties, please [file an issue](https://github.com/aamaricci/SciFortran/issues/new/choose). For any other communication, please reach out any of the contributors or developers:         
[Adriano Amaricci](https://github.com/aamaricci)\
[Lorenzo Crippa](https://github.com/lcrippa)\
[Samuele Giuli](https://github.com/SamueleGiuli)\
[Gabriele Bellomia](https://github.com/beddalumia)\
[Giacomo Mazza](https://github.com/GiacMazza)\
[Francesco Petocchi](mailto:francesco.petocchi@gmail.com)\
[Alberto Scazzola](mailto:alberto.scazzola@polito.it)\
[Massimo Capone](mailto:capone@sissa.it)

--

***LICENSE***  
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL) as published by
the Free Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU LGPL for more details.

You should have received a copy of the GNU LGPL along with this program.  If not, see <http://www.gnu.org/licenses/>.
