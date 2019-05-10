# HP_HEM
High Performance Hybrid Electromagnetic Model, a computer implementation.

Cite the current release: [![DOI](https://zenodo.org/badge/151085118.svg)](https://zenodo.org/badge/latestdoi/151085118)

This library is intended to be easy to use. Interfaces to use it from other programming languages are provided for (work in progress, Julia is the recommended one):

  Julia https://github.com/pedrohnv/hp_hem_julia  
  Mathematica https://github.com/pedrohnv/hp_hem_mathematica  
  Matlab https://github.com/pedrohnv/hp_hem_matlab  
  Octave https://github.com/pedrohnv/hp_hem_octave  

The `examples` folder contains various C files which reproduce results published in the technical literature. Use them as starting point to build your own cases in pure C if you want maximum performance. Examples of use from other programming languages are in their respective repository.

These codes make use of [cubature](https://github.com/stevengj/cubature). It also uses [Intel's MKL](https://software.intel.com/en-us/mkl).

Don't forget to define `#MKLROOT` (use the bash script `mklvars.sh` that is shipped together with intel's MKL). In my case I use the command `source mklvars.sh intel64`.
