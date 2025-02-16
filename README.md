# EDIpack2.0: Massively parallel Exact Diagonalization for generic Quantum Impurity problems

[![TestSuite](https://img.shields.io/github/actions/workflow/status/edipack/EDIpack2.0/PushWorkflow.yml?label=TestSuite&logo=Fortran&style=flat-square)](https://github.com/edipack/EDIpack2.0/actions/workflows/PushWorkflow.yml) 
[![api docs](https://img.shields.io/static/v1?label=API&message=documentation&color=734f96&logo=read-the-docs&logoColor=white&style=flat-square)](https://edipack.github.io/EDIpack2.0/)

<!-- TO BE SETUP ASAP
[![Coverage]()]()
[![api docs](https://img.shields.io/static/v1?label=API&message=documentation&color=734f96&logo=read-the-docs&logoColor=white&style=flat-square)](https://qcmplab.github.io/DMFT_ED)
-->


A suitable extension of [EDIpack](https://github.com/edipack/EDIpack): a  Lanczos based method 
for the solution of generic Quantum Impurity problems,  exploiting distributed memory MPI parallelisation.
This updated version, aims to solve single-site, multi-orbital models, in either  *normal*, *superconducting* (s-wave) or *Spin-non-conserving* (e.g. with Spin-Orbit Coupling or in-plane magnetization) phases, including electron-phonons coupling. The code works at zero and low temperatures.   
 
See [j.cpc.2021.108261](https://doi.org/10.1016/j.cpc.2021.108261) for further information about the underlying algorithms. Yet, suitable modifications have been developed to address the Superconducting and non-SU(2) channels.  

### Install & Use

*EDIpack2.0* builds and get installed using CMake. Loading into the operative systemisis provided by different methods, including module environment.    
Further informations and a guided proceudre are available in the documentation.


### Documentation
All the informations about the structure of the library and its use, together with the Python API *EDIpy2*, are documented at [edipack.github.io/EDIpack2.0/](https://edipack.github.io/EDIpack2.0/)  

NOTE: The documentation is currently under construction. 



### Authors
[Adriano Amaricci](https://github.com/aamaricci)  
[Lorenzo Crippa](https://github.com/lcrippa)  
[Samuele Giuli](https://github.com/SamueleGiuli)  
[Gabriele Bellomia](https://github.com/beddalumia)  
[Giacomo Mazza](https://github.com/GiacMazza)  
[Francesco Petocchi](mailto:francesco.petocchi@gmail.com)  
[Alberto Scazzola](mailto:alberto.scazzola@polito.it)  
[Massimo Capone](mailto:capone@sissa.it)


### Issues
If you encounter bugs or difficulties, please [file an issue](https://github.com/edipack/EDIpack2.0/issues/new/choose). For any other communication, please reach out any of the developers.          
