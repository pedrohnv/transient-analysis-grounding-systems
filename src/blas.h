/* High Performance implementation of the Hybrid Electromagnetic Model
Released under the General Public License 3 (GPLv3).
All parameters' units are in the SI base if omitted.

This file declares which routines from BLAS are used with the appropriate name
mangling.
*/
#ifndef BLAS_H_
#define BLAS_H_

/** ZGEMM performs one of the matrix-matrix operations
C := alpha*op( A )*op( B ) + beta*C,
where  op( X ) is one of
op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H,
alpha and beta are scalars, and A, B and C are matrices, with op( A )
an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
@see [netlib zgemm](http://www.netlib.org/lapack/explore-html/dc/d17/group__complex16__blas__level3_ga4ef748ade85e685b8b2241a7c56dd21c.html#ga4ef748ade85e685b8b2241a7c56dd21c) for full documentation
*/
extern void
zgemm_ (char* transa, char* transb, int* m, int* n, int* k,
        _Complex double* alpha, _Complex double* a, int* lda,
        _Complex double* b, int* ldb, _Complex double* beta,
        _Complex double* c, int* ldc);

/** ZSYMM  performs one of the matrix-matrix operations
C := alpha*A*B + beta*C,
or
C := alpha*B*A + beta*C,
where  alpha and beta are scalars, A is a symmetric matrix and  B and
C are m by n matrices.
@see [netlib zsymm](http://www.netlib.org/lapack/explore-html/dc/d17/group__complex16__blas__level3_ga263a46a500f5c7f04bee1b75ea7f64f6.html#ga263a46a500f5c7f04bee1b75ea7f64f6) for full documentation
*/
extern void
zsymm_ (char* side, char* uplo, int* m, int* n, _Complex double* alpha,
        _Complex double* a, int* lda, _Complex double* b, int* ldb,
        _Complex double* beta, _Complex double* c, int* ldc);

/** ZGEMV  performs one of the matrix-vector operations
y := alpha*A*x + beta*y, or
y := alpha*A**T*x + beta*y, or
y := alpha*A**H*x + beta*y,
where alpha and beta are scalars, x and y are vectors and A is an
m by n matrix.
*/
extern void
zgemv_ (char *trans, int *m, int *n, _Complex double *alpha,
        _Complex double *a, int *lda, _Complex double *x, int *incx,
        _Complex double *beta, _Complex double *y, int *incy);


#endif /* BLAS_H_ */
