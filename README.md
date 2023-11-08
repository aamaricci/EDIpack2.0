# EDIpack2.0: Massively parallel Exact Diagonalization for generic Quantum Impurity problems

[![TestSuite](https://img.shields.io/github/actions/workflow/status/aamaricci/DMFT_ED/PushWorkflow.yml?label=TestSuite&logo=Fortran&style=flat-square)](https://github.com/aamaricci/DMFT_ED/actions/workflows/PushWorkflow.yml) 

<!-- TO BE SETUP ASAP
[![Coverage]()]()
[![api docs](https://img.shields.io/static/v1?label=API&message=documentation&color=734f96&logo=read-the-docs&logoColor=white&style=flat-square)](https://qcmplab.github.io/DMFT_ED)
-->
A suitable extension of ![EDIpack](https://github.com/aamaricci/EDIpack): a  Lanczos based method 
for the solution of generic Quantum Impurity problems,  exploiting distributed memory MPI parallelisation.
This updated version, aims to solve single-site, multi-orbital models, in either  *normal*, *superconducting* (s-wave) or *Spin-non-conserving* (e.g. with Spin-Orbit Coupling or in-plane magnetization) phases, including electron-phonons coupling. The code works at zero and low temperatures.   
 
See [j.cpc.2021.108261](https://doi.org/10.1016/j.cpc.2021.108261) for further information about the underlying algorithms. Yet, suitable modifications have been developed to address the Superconducting and non-SU(2) channels.  


### Dependencies

The code is based on:  

* SciFortran [https://github.com/aamaricci/SciFortran](https://github.com/aamaricci/SciFortran)  

* MPI 

  


### Installation

Installation is available using CMake. In the current v0.0.1 API are only provided in Fortran. In a future release Python and C/C++ API will be included.  
The software gives acces to the static library `libdmft_ed.a` and the related modules `DMFT_ED`

Clone the repo:

`git clone https://github.com/aamaricci/EDIpack2.0`

And from the repository directory (`cd EDIpack2.0`) make a standard out-of-source CMake compilation:

`mkdir build`
`cd build`
`cmake ..`     
`make`     
`make install`   
`make post-install`    

Please follow the instructions on the screen to complete installation on your environment.  
The library can be loaded using one of the following, automatically generated, files :  

* pkg-config file in `~/.pkg-config.d/dmft_ed.pc`  
* environment module file `~/.modules.d/dmft_ed/<PLAT>`  
* homebrew `bash` script `<PREFIX>/bin/configvars.sh`


The `CMake` compilation can be controlled using the following additional variables, default values between `< >`:   

* `-DPREFIX=prefix directory <~/opt/dmft_ed/VERSION/PLAT/[GIT_BRANCH]>` 

* `-DUSE_MPI=<yes>/no`  

* `-DVERBOSE=yes/<no> `  

* `-DBUILD_TYPE=<RELEASE>/TESTING/DEBUG`  


For any information contact the author as:  
adriano DOT amaricci @ gmail DOT com

--

***LICENSE***  
Copyright 2020- (C) Adriano Amaricci, Lorenzo Crippa, Alberto Scazzola, Gabriele Bellomia, Samuele Giuli, Giacomo Mazza, Francesco Petocchi

The software is provided with no license, as such it is protected by copyright.
The software is provided as it is and can be read and copied, in agreement with 
the Terms of Service of GITHUB. Use of the code is constrained to author agreement.   


This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.


