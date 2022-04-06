/* High Performance implementation of the Hybrid Electromagnetic Model
Released under the General Public License 3 (GPLv3).
All parameters' units are in the SI base if omitted.

Routines to build and solve the electrode system that depend on linear
algebra libraries (e.g. LAPACK and BLAS).

Every matrix is stored as a flat array in column major format.
*/
#ifndef LINALG_H_
#define LINALG_H_

#include <complex.h>
#include <string.h>
#include "electrode.h"
#include "blas.h"
#include "lapack.h"
//#include "mkl.h"
//#include "mkl_lapacke.h"

/* ================= Immittance Formulation ================= */
/** Fills the Global Immittance matrix
\f$ W_G = \begin{bmatrix} Y_E & A^T & B^T \\ A & -Z_L & 0 \\ B & 0 & -Z_T \end{bmatrix} \f$
with the incidence matrices \f$ A,B \f$. This function is separated from
fill_impedance because \f$ A,B \f$ depend only on geometry, while \f$ Z_L,Z_T \f$
will vary with frequency. It is expected that \f$ W_G \f$ was zero-ed before
calling this function (or allocated through calloc). As \f$ W_G \f$ is symmetric,
only its lower half is stored (set).
@param wg Global Immittance matrix \f$ W_G \f$ as flat array of size \f$ (2m + n)^2 \f$
@param electrodes array of electrodes
@param num_electrodes number of electrodes \f$ m \f$
@param nodes array of nodes
@param num_nodes number of nodes \f$ n \f$
@return 0 on success
@see fill_impedance_imm
@see solve_immittance
*/
int
fill_incidence_imm (_Complex float* wg, const Electrode* electrodes,
                    size_t num_electrodes, const float* nodes,
                    size_t num_nodes);

/** Fills the Global Immittance matrix
\f$ W_G = \begin{bmatrix} Y_E & A^T & B^T \\ A & -Z_L & 0 \\ B & 0 & -Z_T \end{bmatrix} \f$
with the impedance matrices \f$ Z_L,Z_T \f$. As it is symmetric, only its lower
half is stored (set).
@param wg Global Immittance matrix as flat array of size \f$ (2m + n)^2 \f$
@param zl longitudinal impedance matrix \f$ Z_L \f$ as a flat array of size \f$ m^2 \f$
@param zt transversal impedance matrix \f$ Z_T \f$ as a flat array of size \f$ m^2 \f$
@param num_electrodes number of electrodes \f$ m \f$
@param num_nodes number of nodes \f$ n \f$
@return 0 on success
@see fill_incidence_imm
@see solve_immittance
*/
int
fill_impedance_imm (_Complex float* wg, const _Complex float* zl,
                    const _Complex float* zt, size_t num_electrodes,
                    size_t num_nodes);

/** Solves the system of equations of the Global Immitance formulation.
@param wg global immittance matrix \f$ W_G \f$, assumed to be symmetric with
lower half stored. The LU decomposition is done in-place.
@param ie RHS vector \f$[I_E, 0, 0]^T\f$ where \f$ I_E \f$ are injected currents
in each node and 0 are null vectors of size \f$ m \f$ each. The solution replaces
this array in-place becoming \f$ [U, I_L, I_T] \f$. Where \f$ U \f$ is the nodal
potentials and \f$ I_L \f$ and \f$ I_T \f$ the longitudinal and transversal
currents, respectively.
It is recommended to use the raw LAPACK routines instead of this function.
That way, the same system can be used to solve multiple right hand sides
(injection vectors), the optimal workspace can be queried just once when solving
for multiple frequencies. See the source code of this function to have an idea.
@param num_electrodes number of electrodes \f$ m \f$
@param num_nodes number of nodes \f$ n \f$
@return 0 on sucess
@see fill_incidence_imm
@see fill_impedance_imm
*/
int
solve_immittance (_Complex float* wg, _Complex float* ie, size_t num_electrodes,
                  size_t num_nodes);

/* ================= Admittance Formulation ================= */
/** Fills the nodal incidence matrices \f$ A \f$ and \f$ B \f$ used in the
calculation of the nodal admittance matrix \f$ Y_N = A^T Z_L^{-1} A + B^T Z_T^{-1} B \f$.\n
Note that \f$ B = |A|/2 \f$ and if it is desired to save memory, A and B can be
set to be the same pointer and, in that case, it is guaranteed that the end result
will be \f$ A \f$. In other words: when only \f$ A \f$ is desired, \f$ B \f$ can
point to \f$ A, (B := A)\f$, saving memory space. Then call the function as:
```fill_incidence_adm(a, a, electrodes, m, nodes, n);``` so that \f$ a:=A \f$.
@param a longitudinal incidence matrix \f$ A \f$ as flat array of size
\f$ m \times n \f$ in COLUMN MAJOR ordering
@param b transversal incidence matrix \f$ B \f$ as flat array of size
\f$ m \times n \f$ in COLUMN MAJOR ordering
@param electrodes array of electrodes
@param num_electrodes number of electrodes \f$ m \f$
@param nodes array of nodes
@param num_nodes number of nodes \f$ n \f$
@return 0 on success
@see fill_impedance_adm
@see fill_impedance_adm2
*/
int
fill_incidence_adm (_Complex float* a, _Complex float* b,
                    const Electrode* electrodes, size_t num_electrodes,
                    const float* nodes, size_t num_nodes);

/** Calculates the nodal admittance matrix \f$ Y_N = A^T Z_L^{-1} A + B^T Z_T^{-1} B \f$.
As it is symmetric, only its lower half is stored (set).
@param yn nodal admmitance matrix \f$ Y_N \f$ as flat array of size \f$ n^2 \f$
@param a longitudinal incidence matrix \f$ A \f$ as flat array of size
\f$ m \times n \f$ in COLUMN MAJOR ordering
@param b transversal incidence matrix \f$ B \f$ as flat array of size
\f$ m \times n \f$ in COLUMN MAJOR ordering
@param zl longitudinal impedance matrix \f$ Z_L \f$ as a flat array of size \f$ m^2 \f$
@param zt transversal impedance matrix \f$ Z_T \f$ as a flat array of size \f$ m^2 \f$
@param num_electrodes number of electrodes \f$ m \f$
@param num_nodes number of nodes \f$ n \f$
@return 0 on success
@see fill_impedance_adm2
@see solve_admittance
*/
int
fill_impedance_adm (_Complex float* yn, _Complex float* zl, _Complex float* zt,
                    _Complex float* a, _Complex float* b, size_t num_electrodes,
                    size_t num_nodes);

/** Calculates the modified longitudinal and transversal admittances \f$ Z_L^{-1} A \f$
and \f$ Z_T^{-1} B \f$. \f$ Z_L \f$ and \f$ Z_T \f$ are replaced by their inverses
in the process. Useful when they are needed to calculate the currents
\f$ I_L = Z_L^{-1} A U \f$ and \f$ I_T = Z_T^{-1} B U \f$.
@param yla modified longitudinal admittance \f$ Z_L^{-1} A \f$ of size \f$ m \times n \f$
@param ytb modified transversal admittance \f$ Z_T^{-1} B \f$ of size \f$ m \times n \f$
@param zl longitudinal impedance matrix \f$ Z_L \f$ as a flat array of size \f$ m^2 \f$
@param zt transversal impedance matrix \f$ Z_T \f$ as a flat array of size \f$ m^2 \f$
@param a longitudinal incidence matrix \f$ A \f$ as flat array of size
\f$ m \times n \f$ in COLUMN MAJOR ordering
@param b transversal incidence matrix \f$ B \f$ as flat array of size
\f$ m \times n \f$ in COLUMN MAJOR ordering
@param num_electrodes number of electrodes \f$ m \f$
@param num_nodes number of nodes \f$ n \f$
@return 0 on success
@see fill_impedance_adm2
*/
int
calculate_yla_ytb (_Complex float* yla, _Complex float* ytb,
                   _Complex float* zl, _Complex float* zt,
                   _Complex float* a, _Complex float* b,
                   size_t num_electrodes, size_t num_nodes);

/** Calculates the nodal admittance matrix \f$ Y_N = A^T Z_L^{-1} A + B^T Z_T^{-1} B \f$,
were \f$ Z_L^{-1} A \f$ and \f$ Z_T^{-1} B \f$ are already calculated. Those
matrices are preserved. Useful when they must be calculated to find the currents
\f$ I_L = Z_L^{-1} A U \f$ and \f$ I_T = Z_T^{-1} B U \f$.
@param yn nodal admmitance matrix \f$ Y_N \f$ as flat array of size \f$ n^2 \f$
@param a longitudinal incidence matrix \f$ A \f$ as flat array of size
\f$ m \times n \f$ in COLUMN MAJOR ordering
@param b transversal incidence matrix \f$ B \f$ as flat array of size
\f$ m \times n \f$ in COLUMN MAJOR ordering
@param yla modified longitudinal admittance \f$ Z_L^{-1} A \f$
@param ytb modified transversal admittance \f$ Z_T^{-1} B \f$
@param num_electrodes number of electrodes \f$ m \f$
@param num_nodes number of nodes \f$ n \f$
@return 0 on success
@see fill_impedance_adm
@see calculate_yla_ylb
@see solve_admittance
*/
int
fill_impedance_adm2 (_Complex float* yn, _Complex float* yla,
                     _Complex float* ytb, _Complex float* a,
                     _Complex float* b, size_t num_electrodes,
                     size_t num_nodes);

/** Solves the system of equations of the Nodal Admitance formulation
\f$ Y_N U = I_E \f$.
It is recommended to use the raw LAPACK routines instead of this function.
That way, the same system can be used to solve multiple right hand sides
(injection vectors), the optimal workspace can be queried just once when solving
for multiple frequencies. See the source code of this function to have an idea.
@param yn nodal admittance matrix \f$ Y_N = A^T Z_L^{-1} A + B^T Z_T^{-1} B \f$,
assumed to be symmetric with lower half stored. The LU decomposition is done in-place.
@param ie RHS vector of injected currents in each node. The solution replaces
this array in-place.
@param num_electrodes number of electrodes \f$ m \f$
@param num_nodes number of nodes \f$ n \f$
@return 0 on sucess
@see fill_incidence_adm
@see fill_impedance_adm
@see fill_impedance_adm2
*/
int
solve_admittance (_Complex float* yn, _Complex float* ic, size_t num_nodes);

#endif /* LINALG_H_ */
