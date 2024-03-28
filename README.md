[See this Code Ocean capsule for a running example.](https://codeocean.com/capsule/0309688/tree)

# TAGS
Formely HP-HEM (High Performance implementation of the Hybrid Electromagnetic Model and its variants), but we changed the name to TAGS (Transient Analysis of Grounding Systems) to avoid any confusion with Hewlett-Packard Company.

The documentation of the code can be read at https://pedrohnv.github.io/transient-analysis-grounding-systems

The Hybrid Electromagnetic Model can be used to simulate any situation where the
thin wire hypothesis is valid. See [1] for the mathematical formulation.

[1] VISACRO, S.; SOARES, A. HEM: A model for simulation of lightning-related engineering problems. IEEE Transactions on power delivery, v. 20, n. 2, p. 1206-1208, 2005.

Cite the current release: [![DOI](https://zenodo.org/badge/151085118.svg)](https://zenodo.org/badge/latestdoi/151085118)

The permanent DOI is https://doi.org/10.5281/zenodo.2644010, which always resolves to the last release.

There is a pure Julia version for those who don't want to use C: https://github.com/pedrohnv/transient-analysis-grounding-systems-julia

There is also an interface to use a dynamic library from Julia (outdated, but it is a starting point): https://github.com/pedrohnv/hp_hem_julia

The `examples` folder contains various C files which reproduce results published in the technical literature. Use them as starting point to build your own cases in pure C if you want maximum performance. Examples of use from other programming languages are in their respective repository.

Dependencies:
 - A C99 compliant compiler (such as gcc, try MinGW on Windows - Visual Studio won't work)
 - [Cubature](https://github.com/stevengj/cubature)
 - [OpenBLAS](https://www.openblas.net/)
 - [LAPACK](http://www.netlib.org/lapack/)
 - [FFTW3](http://www.fftw.org/)
