# HP_HEM
High performance Hybrid Electromagnetic Model implementation.

These codes make use of [cubature](https://github.com/stevengj/cubature). Put it in the source directory like:
```
/HP_HEM
  /cubature
  other files
```
It is assumed that the user has an optmized [BLAS](https://www.netlib.org/blas/) in their machine (install [ATLAS](http://math-atlas.sourceforge.net/) if your are unsure).
This also makes use of [slatec](http://www.netlib.org/slatec/) and [LAPACK](http://www.netlib.org/lapack/), both assumed to be in a path your compiler can find (add paths in the Makefile, if you need)..
