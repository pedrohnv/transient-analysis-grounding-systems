/** High Performance Hybrid Electromagnetic Model calculations in C.

All parameters' units are in the SI base units if omitted.

Routines to build and solve the electrode system that depend on linear
algebra libraries, e.g., Intel MKL or LAPACK.
*/
#ifndef LINALG_H_
#define LINALG_H_

#include "electrode.h"
#include <complex.h>
#include <string.h>
#include "mkl_types.h"
#define MKL_Complex16 _Complex double //overwrite type

// Imitance matrix WE building
/** fill_incidence
Fills the imitance matrix `we = [[ZL/2, -ZL/2, A], [ZT, ZT, B], [C, D, Yn]]`
with the incidence matrices `A`, `B`, `C` and `D`. This function is separated from
fill_impedance because `A, B, C, D` only depends on geometry, while `ZT` and `ZL`
will vary with frequency.
@param we imitance matrix as flat array of size
`(2*num_electrodes + num_nodes)^2`
@param electrodes array of electrodes
@param num_electrodes number of electrodes
@param nodes array of nodes
@param num_nodes number of nodes
@return 0 on success
@see fill_impedance
*/
int
fill_incidence (_Complex double *we, const Electrode *electrodes,
                size_t num_electrodes, const double nodes[][3], size_t num_nodes);

/** fill_impedance
Fills the imitance matrix `we` = `[[ZL/2, -ZL/2, A], [ZT, ZT, B], [C, D, Yn]]`
with the impedance matrices `ZT` and `ZL`, and the nodal admittance Yn.
@param we imitance matrix as flat array of size
`(2*num_electrodes + num_nodes)^2`
@param electrodes array of electrodes
@param num_electrodes number of electrodes
@param num_nodes number of nodes
@param zt transversal impedance matrix as a flat array of size
`num_electrodes^2`
@param zl longitudinal impedance matrix as a flat array of size
`num_electrodes^2`
@param yn nodal admittance matrix as a flat array of size `num_nodes^2`
@return 0 on success
@see fill_incidence
*/
int
fill_impedance (_Complex double *we, const Electrode *electrodes,
                size_t num_electrodes, size_t num_nodes, const _Complex double *zl,
                const _Complex double *zt, const _Complex double *yn);

/** incidence_alt
Alternative approach using an equivalent Ynodal
*/
int
incidence_alt (double *a, double *b, const Electrode *electrodes,
               size_t num_electrodes, const double nodes[][3], size_t num_nodes);

/** solve_electrodes
Solves the system of equations that defines the electrode system.
Uses Intel MKL LAPACKE_zgesv.
@param we imitance matrix `[[ZL/2, -ZL/2, A], [ZT, ZT, B], [C, D, Yn]]`. The
LU decomposition is done in-place on we
@param ie RHS vector `[0, 0, ic]^T` where `ic` are injected currents in each node
and 0 are null vectors of size `num_electrodes` each. The solution replaces
this array in-place
@param num_electrodes number of electrodes
@param num_nodes number of nodes
@return 0 on sucess
@see https://software.intel.com/en-us/mkl-developer-reference-c-gesv
*/
int
solve_electrodes (_Complex double *we, _Complex double *ie, size_t num_electrodes,
                  size_t num_nodes);

/** ynodal_eq
Finds the equivalent nodal admittance matrix of the electrode system.
TODO put const qualifier to a, b and zl?
*/
int
ynodal_eq (_Complex double *yn, const double *a, const double *b,
           _Complex double *zl, _Complex double *zt, size_t num_electrodes,
           size_t num_nodes);

/** harmonic_impedance1
Calculates the harmonic impedance of a copper electrode system in a
two layer medium. Electrodes are considered to be in medium 1.
No segmentation of the electrodes is done.
Injection node is considered the first.
@param ns number of frequencies
@param s array of complex frequencies `c + I*w`
@param kappa1 medium 1 complex conductivity `(sigma + I*w*eps)` in S/m
@param kappa2 medium 2 complex conductivity `(sigma + I*w*eps)` in S/m
@param gamma1 medium 1 propagation constant
@param electrodes array of electrodes
@param num_electrodes number of electrodes
@param nodes array of nodes
@param num_nodes number of nodes
@param max_eval specifies a maximum number of function evaluations (0 for no
limit)
@param reqAbsError the absolute error requested (0 to ignore)
@param reqRelError the relative error requested (0 to ignore)
@param error_norm (enumeration defined in cubature.h) error checking scheme
@param rsource source internal resistence to consider.
@param zh harmonic impedance array of size `ns`
@return 0 on success
*/
int
harmonic_impedance1 (size_t ns, const _Complex double *s,
                     const _Complex double *kappa1, const _Complex double *kappa2,
                     const _Complex double *gamma1, const Electrode *electrodes,
                     const Electrode *images, size_t num_electrodes,
                     const double nodes[][3], size_t num_nodes, size_t max_eval,
                     double req_abs_error, double req_rel_error, int error_norm,
                     double rsource, _Complex double *zh);

int
harmonic_impedance1_alt (size_t ns, const _Complex double *s,
                         const _Complex double *kappa1,
                         const _Complex double *kappa2,
                         const _Complex double *gamma1,
                         const Electrode *electrodes, const Electrode *images,
                         size_t num_electrodes, const double nodes[][3],
                         size_t num_nodes, size_t max_eval, double req_abs_error,
                         double req_rel_error, int error_norm, double rsource,
                         _Complex double *zh);

#endif /* LINALG_H_ */
