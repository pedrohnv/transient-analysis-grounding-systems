/** High Performance Hybrid Electromagnetic Model calculations in C.

All parameters' units are in the SI base units if omitted.

Routines to build and solve the electrode system that depend on linear
algebra libraries, e.g., LAPACK and BLAS.

Every matrix is stored as a flat array in column major format.
*/
#ifndef LINALG_H_
#define LINALG_H_

#include <complex.h>
#include <string.h>
#include "electrode.h"

// BLAS routines ==========
/**
ZGEMM  performs one of the matrix-matrix operations
C := alpha*op( A )*op( B ) + beta*C,
where  op( X ) is one of
op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H,
alpha and beta are scalars, and A, B and C are matrices, with op( A )
an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
@see http://www.netlib.org/lapack/explore-html/
*/
extern void
zgemm_ (char *transa, char *transb, int *m, int *n, int *k,
        _Complex double *alpha, _Complex double *a, int *lda,
        _Complex double *b, int *ldb, _Complex double *beta,
        _Complex double *c, int *ldc);

/**
ZSYMM  performs one of the matrix-matrix operations
C := alpha*A*B + beta*C,
or
C := alpha*B*A + beta*C,
where  alpha and beta are scalars, A is a symmetric matrix and  B and
C are m by n matrices.
@see http://www.netlib.org/lapack/explore-html/
*/
extern void
zsymm_ (char *side, char *uplo, int *m, int *n, _Complex double *alpha,
        _Complex double *a, int *lda, _Complex double *b, int *ldb,
        _Complex double *beta, _Complex double *c, int *ldc);

// LAPACK routines =========
/**
ZSYSV computes the solution to a complex system of linear equations
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
zsysv_ (char *uplo, int *n, int *nrhs, _Complex double *a, int *lda, int *ipiv,
        _Complex double *b, int *ldb, _Complex double *work, int *lwork, int *info);

/**
ZGESV computes the solution to a complex system of linear equations
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
zgesv_ (int *n, int *nrhs, _Complex double *a, int *lda, int* ipiv,
        _Complex double *b, int *ldb, int* info);

/**
ZGETRF computes an LU factorization of a general M-by-N matrix A
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
zgetrf_ (int *m, int *n, _Complex double *a, int *lda, int *ipiv, int *info);

/**
ZGETRI computes the inverse of a matrix using the LU factorization
computed by ZGETRF.
This method inverts U and then computes inv(A) by solving the system
inv(A)*L = inv(U) for inv(A).
@see http://www.netlib.org/lapack/explore-html/
*/
extern void
zgetri_ (int *n, _Complex double *a, int *lda, int *ipiv,
         _Complex double *work, int *lwork, int *info);

/**
ZSYTRF computes the factorization of a complex symmetric matrix A
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
zsytrf_ (char *uplo, int *n, _Complex double *a, int *lda, int *ipiv,
         _Complex double *work, int *lwork, int *info);

/**
ZSYTRI computes the inverse of a complex symmetric indefinite matrix
A using the factorization A = U*D*U**T or A = L*D*L**T computed by
ZSYTRF.
@see http://www.netlib.org/lapack/explore-html/
*/
extern void
zsytri_ (char *uplo, int *n, _Complex double *a, int *lda, int *ipiv,
         _Complex double *work, int *info);

// Immittance =================
/** fill_incidence_imm
Fills the immittance matrix `WE = [[Ye, C, D], [A, ZL/2, -ZL/2], [B, ZT, ZT]]`
with the incidence matrices `A`, `B`, `C` and `D`. This function is separated from
fill_impedance because `A, B, C, D` only depends on geometry, while `ZT` and `ZL`
will vary with frequency.
@param we immittance matrix as flat array of size
`(2*num_electrodes + num_nodes)^2`
@param electrodes array of electrodes
@param num_electrodes number of electrodes
@param nodes array of nodes
@param num_nodes number of nodes
@return 0 on success
@see fill_impedance_imm
@see solve_immittance
*/
int
fill_incidence_imm (_Complex double *we, const Electrode *electrodes,
                    size_t num_electrodes, const double *nodes,
                    size_t num_nodes);

int
fill_incidence_imm2 (_Complex double *we, const Electrode *electrodes,
                     size_t num_electrodes, const double *nodes,
                     size_t num_nodes);

/** fill_impedance_imm
Fills the immittance matrix `WE` = `[[Ye, C, D], [A, ZL/2, -ZL/2], [B, ZT, ZT]]`
with the impedance matrices `ZT` and `ZL`, and the nodal admittance Yn.
@param we immittance matrix as flat array of size
`(2*num_electrodes + num_nodes)^2`
@param num_electrodes number of electrodes
@param num_nodes number of nodes
@param zt transversal impedance matrix as a flat array of size
`num_electrodes^2`
@param zl longitudinal impedance matrix as a flat array of size
`num_electrodes^2`
@param ye "external" nodal admittance matrix as a flat array of size `num_nodes^2`
@return 0 on success
@see fill_incidence_we
@see fill_incidence_yn
@see solve_immittance
*/
int
fill_impedance_imm (_Complex double *we, size_t num_electrodes, size_t num_nodes,
                    const _Complex double *zl, const _Complex double *zt,
                    const _Complex double *ye);

/** solve_immittance
Solves the system of equations of the Global Immitance formulation.
@param we global immittance matrix `[[Ye, C, D], [A, ZL/2, -ZL/2], [B, ZT, ZT]]`.
The LU decomposition is done in-place.
@param ie RHS vector `[ic, 0, 0]^T` where `ic` are injected currents in each node
and 0 are null vectors of size `num_electrodes` each. The solution replaces
this array in-place becoming `[u, i1, i2]`. Where `u` is the nodal potentials
and the transversal and longitudinal currents, `IT` and `IL` respectively,
can then be calculated as
        IT = i1 + i2
        IL = (i1 - i2)/2
@param num_electrodes number of electrodes
@param num_nodes number of nodes
@return 0 on sucess
@see solve_admittance
*/
int
solve_immittance (_Complex double *we, _Complex double *ie, size_t num_electrodes,
                  size_t num_nodes);

// Admittance =================
/** fill_incidence_adm
Fills the nodal incidence matrices `α` and `β` used in the calculation of the nodal
admittance matrix `Yn` = `α^T.YT.α + β^T.YL.β`
with the impedance matrices `ZT` and `ZL`, and the nodal admittance Yn.
@param a transversal incidence `α` of size `num_electrodes * num_nodes`
@param b longitudinal incidence `β` of size `num_electrodes * num_nodes`
@param electrodes array of electrodes
@param num_electrodes number of electrodes
@param nodes array of nodes
@param num_nodes number of nodes
@return 0 on success

*/
int
fill_incidence_adm (double *a, double *b, const Electrode *electrodes,
                    size_t num_electrodes, const double *nodes,
                    size_t num_nodes);

/** fill_impedance_adm
Calculates the nodal admittance matrix `Yn` = `α^T.YT.α + β^T.YL.β`
@param yn nodal admmitance matrix as flat array of size `num_nodes^2`
@param a transversal incidence `α`
@param b longitudinal incidence `β`
@param zl longitudinal impedance matrix as a flat array of size
`num_electrodes^2`
@param zt transversal impedance matrix as a flat array of size
`num_electrodes^2`
@param num_electrodes number of electrodes
@param num_nodes number of nodes
@param ye "external" nodal admittance matrix as a flat array of size `num_nodes^2`
@return 0 on success
@see fill_incidence_we
@see fill_incidence_yn
@see solve_admittance
*/
int
fill_impedance_adm (_Complex double *yn, const double *a, const double *b,
                    _Complex double *zl, _Complex double *zt, size_t num_electrodes,
                    size_t num_nodes, const _Complex double *ye);

/** solve_admittance
Solves the system of equations of the Nodal Admitance formulation.
@param yn nodal admittance matrix `[[Ye, C, D], [A, ZL/2, -ZL/2], [B, ZT, ZT]]`. The
LU decomposition is done in-place on we
@param ic RHS vector of injected currents in each node. The solution replaces
this array in-place.
@param num_electrodes number of electrodes
@param num_nodes number of nodes
@return 0 on sucess
@see solve_immittance
*/
int
solve_admittance (_Complex double *yn, _Complex double *ic, size_t num_nodes);

/** zh_immittance
Calculates the harmonic impedance of a copper electrode system buried in a
single layer medium. No segmentation of the electrodes is done.
This calculation is done considering the current injected at every node,
one at a time. Immitance formulation is used with double integral.
@param ns number of frequencies
@param s array of complex frequencies `c + I*w`
@param sigma medium conductivity in S/m
@param epsr medium relative permittivity
@param mur medium relative permeability
@param electrodes array of electrodes
@param images array of images
@param num_electrodes number of electrodes
@param nodes array of nodes
@param num_nodes number of nodes
@param max_eval specifies a maximum number of function evaluations (0 for no
limit)
@param reqAbsError the absolute error requested (0 to ignore)
@param reqRelError the relative error requested (0 to ignore)
@param zh harmonic impedance matrix `num_nodes x ns` as a flat array
@return 0 on success
*/
int
zh_immittance (size_t ns, const _Complex double *s, double sigma, double epsr,
               double mur, const Electrode *electrodes, const Electrode *images,
               size_t num_electrodes, const double *nodes, size_t num_nodes,
               size_t max_eval, double req_abs_error, double req_rel_error,
               _Complex double *zh);

/** sim_immittance
Calculates the node potentials and conductors' transversal and longitudinal
currents of a copper electrode system buried in a single layer medium.
No segmentation of the electrodes is done.
This calculation is done considering the current injected at a single node `n`.
Immitance formulation is used with double integral.
@param ns number of frequencies
@param s array of complex frequencies `c + I*w`
@param sigma medium conductivity in S/m
@param epsr medium relative permittivity
@param mur medium relative permeability
@param electrodes array of electrodes
@param images array of images
@param num_electrodes number of electrodes
@param nodes array of nodes
@param num_nodes number of nodes
@param max_eval specifies a maximum number of function evaluations (0 for no
limit)
@param reqAbsError the absolute error requested (0 to ignore)
@param reqRelError the relative error requested (0 to ignore)
@param inj_node node number in which the current is injected
@param inj_current array of size `ns` of injected currents for each `s`
@param inj_adm array of size `ns` of the source admittance for each `s`
@param u node potential matrix `num_nodes x ns` as a flat array
@param il longitudinal currents matrix `num_electrodes x ns` as a flat array
@param it transversal currents matrix `num_electrodes x ns` as a flat array
@return 0 on success
*/
int
sim_immittance (size_t ns, const _Complex double *s, double sigma, double epsr,
                double mur, const Electrode *electrodes, const Electrode *images,
                size_t num_electrodes, const double *nodes, size_t num_nodes,
                size_t max_eval, double req_abs_error, double req_rel_error,
                size_t inj_node, const _Complex double *inj_current,
                const _Complex double *inj_adm, _Complex double *u,
                _Complex double *il, _Complex double *it);

#endif /* LINALG_H_ */
