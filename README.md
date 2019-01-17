# HP_HEM
High performance Hybrid Electromagnetic Model implementation.

This library is intended to be easy to use. Interfaces to use it from other programming language are provided (work in progress) for Julia, Mathematica, Matlab, Octave, Python.

The `examples` folder contains various C files which reproduce results published in the technical literature. Use them as starting point to build your own cases in pure C if you want maximum performance. Examples of uses from other programming languages are in their respective interface subfolder: `interfaces/language`.

These codes make use of [cubature](https://github.com/stevengj/cubature) and [SLATEC](http://netlib.org/slatec/). It also uses [Intel's MKL](https://software.intel.com/en-us/mkl).

It is assumed that the user has an optmized [BLAS](https://www.netlib.org/blas/) in their machine (install [ATLAS](http://math-atlas.sourceforge.net/) if you are unsure).

Don't forget to define `#MKLROOT` (use the bash script `mklvars.sh` that is shipped together with intel's MKL). In my case I use the command `source mklvars.sh intel64`.
