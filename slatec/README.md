SLATEC Common Mathematical Library
==================================

**Quick links:** [guide][gde], [releases][rel].

[SLATEC][slt] is a comprehensive software library containing over 1400 general
purpose mathematical and statistical routines written in Fortran 77.

For technical reasons, the library is broken into three parts:

  - `libslatec`: contains most of the functions, except for the routines for
	solving boundary-value problems.

  - `libslatec-dbvp`: contains double-precision routines for solving
    boundary-value problems.  They are kept separate from the rest of the
    library as they require the user to define the input functions or will
    fail during link time.

  - `libslatec-sbvp`: similar to `libslatec-dbvp` but for single-precision.

# Installation

Run the `Makefile` using:

    make FC=gfortran all

where `gfortran` may be replaced by the name of your Fortran 77 compiler.

[gde]: https://raw.githubusercontent.com/Rufflewind/slatec/master/guide
[rel]: https://github.com/Rufflewind/slatec/releases
[slt]: http://netlib.org/slatec
