# HP-HEM
High Performance implementation of the Hybrid Electromagnetic Model and its variants

The Hybrid Electromagnetic Model can be used to simulate any situation where the
thin wire hypothesis is valid. See [1] for the mathematical formulation.

[1] VISACRO, S.; SOARES, A. HEM: A model for simulation of lightning-related engineering problems. IEEE Transactions on power delivery, v. 20, n. 2, p. 1206-1208, 2005.

Cite the current release: [![DOI](https://zenodo.org/badge/151085118.svg)](https://zenodo.org/badge/latestdoi/151085118)

The permanent DOI is https://doi.org/10.5281/zenodo.2644010, which always resolves to the last release.

This library is intended to be easy to use. Interfaces to use it from other programming languages are provided for (work in progress, Julia is the recommended one):

  Julia https://github.com/pedrohnv/hp_hem_julia
  Matlab https://github.com/pedrohnv/hp_hem_matlab

The `examples` folder contains various C files which reproduce results published in the technical literature. Use them as starting point to build your own cases in pure C if you want maximum performance. Examples of use from other programming languages are in their respective repository.

Dependencies:
 - [Cubature](https://github.com/stevengj/cubature)
 - [OpenBLAS](https://www.openblas.net/)
 - [LAPACK](http://www.netlib.org/lapack/)
 - [FFTW3](http://www.fftw.org/)
