/* High Performance implementation of the Hybrid Electromagnetic Model
Released under the General Public License 3 (GPLv3).
All parameters' units are in the SI base if omitted.

This file declares which routines from LAPACK are used.
*/
#ifndef LAPACK_H_
#define LAPACK_H_

/** ZSYSV computes the solution to a complex system of linear equations
A * X = B,
where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
matrices.
The diagonal pivoting method is used to factor A as
A = U * D * U**T,  if UPLO = 'U', or
A = L * D * L**T,  if UPLO = 'L',
where U (or L) is a product of permutation and unit upper (lower)
triangular matrices, and D is symmetric and block diagonal with
1-by-1 and 2-by-2 diagonal blocks.  The factored form of A is then
used to solve the system of equations A * X = B.
@see http://www.netlib.org/lapack/explore-html/
*/
extern void
zsysv_ (char* uplo, int* n, int* nrhs, _Complex double* a, int* lda, int* ipiv,
        _Complex double* b, int* ldb, _Complex double* work, int* lwork, int* info);

extern void
csysv_ (char* uplo, int* n, int* nrhs, _Complex float* a, int* lda, int* ipiv,
        _Complex float* b, int* ldb, _Complex float* work, int* lwork, int* info);

/** ZGESV computes the solution to a complex system of linear equations
A * X = B,
where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
The LU decomposition with partial pivoting and row interchanges is
used to factor A as
A = P * L * U,
where P is a permutation matrix, L is unit lower triangular, and U is
upper triangular.  The factored form of A is then used to solve the
system of equations A * X = B.
@see http://www.netlib.org/lapack/explore-html/
*/
extern void
zgesv_ (int* n, int* nrhs, _Complex double* a, int* lda, int* ipiv,
        _Complex double* b, int* ldb, int* info);

extern void
cgesv_ (int* n, int* nrhs, _Complex float* a, int* lda, int* ipiv,
        _Complex float* b, int* ldb, int* info);

/** ZGETRF computes an LU factorization of a general M-by-N matrix A
using partial pivoting with row interchanges.
The factorization has the form
A = P * L * U
where P is a permutation matrix, L is lower triangular with unit
diagonal elements (lower trapezoidal if m > n), and U is upper
triangular (upper trapezoidal if m < n).
This is the right-looking Level 3 BLAS version of the algorithm.
@see http://www.netlib.org/lapack/explore-html/
*/
extern void
zgetrf_ (int* m, int* n, _Complex double* a, int* lda, int* ipiv, int* info);

extern void
cgetrf_ (int* m, int* n, _Complex float* a, int* lda, int* ipiv, int* info);

/** ZGETRI computes the inverse of a matrix using the LU factorization
computed by ZGETRF.
This method inverts U and then computes inv(A) by solving the system
inv(A)*L = inv(U) for inv(A).
@see http://www.netlib.org/lapack/explore-html/
*/
extern void
zgetri_ (int* n, _Complex double* a, int* lda, int* ipiv,
         _Complex double* work, int* lwork, int* info);

extern void
cgetri_ (int* n, _Complex float* a, int* lda, int* ipiv,
         _Complex float* work, int* lwork, int* info);

/** ZSYTRF computes the factorization of a complex symmetric matrix A
using the Bunch-Kaufman diagonal pivoting method.  The form of the
factorization is
A = U*D*U**T  or  A = L*D*L**T
where U (or L) is a product of permutation and unit upper (lower)
triangular matrices, and D is symmetric and block diagonal with
with 1-by-1 and 2-by-2 diagonal blocks.
This is the blocked version of the algorithm, calling Level 3 BLAS.
@see http://www.netlib.org/lapack/explore-html/
*/
extern void
zsytrf_ (char* uplo, int* n, _Complex double* a, int* lda, int* ipiv,
         _Complex double* work, int* lwork, int* info);

extern void
csytrf_ (char* uplo, int* n, _Complex float* a, int* lda, int* ipiv,
         _Complex float* work, int* lwork, int* info);

/** ZSYTRI computes the inverse of a complex symmetric indefinite matrix
A using the factorization A = U*D*U**T or A = L*D*L**T computed by
ZSYTRF.
@see http://www.netlib.org/lapack/explore-html/
*/
extern void
zsytri_ (char* uplo, int* n, _Complex double* a, int* lda, int* ipiv,
         _Complex double* work, int* info);

extern void
csytri_ (char* uplo, int* n, _Complex float* a, int* lda, int* ipiv,
         _Complex float* work, int* info);

#endif /* LAPACK_H_ */
