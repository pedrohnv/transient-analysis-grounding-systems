/**
High performance Hybrid Electromagnetic Model calculations in C.

All parameters' units are in SI if omitted.

@author Pedro Henrique Nascimento Vieira

TODO insert condition to check if sender == receiver?
    during integration and calculate distance from center radius
*/

#ifndef ELECTRODE_H_
#define ELECTRODE_H_

#include <complex.h>
#include <stdlib.h>
#include <mkl_types.h>
#define MKL_Complex16 _Complex double //overwrite type

//default integration options ===============================================
/**
Type of integration and simplification thereof to be done.
@param NONE integration is not done
@param INTG_DOUBLE \f$ \int_0^{L_s} \int_0^{L_r} \frac{e^{-\gamma r}}{r} dl_r dl_s \f$
@param INTG_EXP_LOGNF \f$ \int_0^{L_r} e^{-\gamma \bar r} Log(N_f) dl_r \f$
@param INTG_LOGNF \f$ e^{-\gamma \bar r} \int_0^{L_r} Log(N_f) dl_r \f$
*/
enum Integration_type
{
    NONE,
    INTG_DOUBLE,
    INTG_EXP_LOGNF,
    INTG_LOGNF
};

/**
Structure that defines an electrode.
@param start_point array `(x,y,z)` that defines the starting point of the
electrode
@param end_point array `(x,y,z)` defining the ending point
@param middle_point array `(x,y,z)` of the middle point:
`(start_point + end_point)/2`
@param length electrode length `Norm(start_point - end_point)`
@param radius electrode radius
@param zi total internal impedance of the electrode (Ohm)
*/
typedef struct
{
    double start_point[3];
    double end_point[3];
    double middle_point[3];
    double length;
    double radius;
    _Complex double zi;
} Electrode;

/**
Structure to make the integration ot the "potential" between two electrodes.
@param sender electrode that generates the excitation
@param receiver electrode that is excitated
@param gamma complex propagation constant of the medium
@param simplified integration being done is INTG_EXP_LOGNF (false) or
INTG_LOGNF (true); no use otherwise
@see Integration_type
*/
typedef struct {
    Electrode* sender;
    Electrode* receiver;
    _Complex double gamma;
    int simplified;
} Integration_data;

/**
Populates an Electrode structure.
@param electrode pointer to an allocated memory
@param start_point array `(x,y,z)` that defines the starting point of the
electrode
@param end_point array `(x,y,z)` defining the ending point
@param radius electrode radius
@param zi total internal impedance of the electrode
@return 0 on success
*/
int populate_electrode(
    Electrode* electrode, double start_point[3], double end_point[3],
    double radius, _Complex double zi);

/**
Segments an electrode (conductor) populating a passed array of electrodes and
an array of nodes.
@param electrodes pointer to an array of Electrode to be filled
@param nodes array[`(num_segments + 1)`][`3`] to be filled by the new nodes
that are created
@param num_segments number of segments to do
@param start_point electrode's start point
@param end_point electrode's end point
@param radius electrode's radius
@param unit_zi internal impedance per unit length (Ohm/m)
@return 0 on success
*/
int segment_electrode(
    Electrode* electrodes, double nodes[][3], int num_segments,
    double* start_point, double* end_point, double radius,
    _Complex double unit_zi);

/**
Segments a group o Electrodes.
*/
int segment_group();

/**
Calculates the integrand \f$ \frac{e^{-\gamma r}}{r} \f$ to be integrated
using Cubature.
@param ndim must be = 2
@param t array of electrodes length's percentage (0 to 1)
@param auxdata Integration_data
@param fdim must be = 2
@param fval pointer to where the result is stored
@return 0 on success
@see Integration_type
@see Integration_data
@see integral
@see https://github.com/stevengj/cubature
*/
int integrand_double(
    unsigned ndim, const double *t, void *auxdata, unsigned fdim,double *fval);

/**
Calculates the simplified integrand \f$ e^{-\gamma \bar r} Log(N_f) \f$ to be used
in Cubature integration.
@param ndim must be = 1
@param t receiver's length percentage (0 to 1)
@param auxdata Integration_data
@param fdim must be = 2
@param fval pointer to where the result is stored
@return 0 on success
@see Integration_type
@see Integration_data
@see integral
@see https://github.com/stevengj/cubature
*/
int exp_logNf(
    unsigned ndim, const double *t, void *auxdata, unsigned fdim, double *fval);

/**
Calculates the integral along the sender and receiver using Cubature.
@param sender Electrode
@param receiver Electrode
@param gamma medium propagation constant
@param max_eval specifies a maximum number of function evaluations (0 for no
limit)
@param reqAbsError the absolute error requested (0 to ignore)
@param reqRelError the relative error requested (0 to ignore)
@param error_norm (enumeration defined in cubature.h) error checking scheme
@param integration_type type of integration to be done. If NONE, the value
`result[0] + result[1]*I` will be unaltered and no integration is done.
@param result array[2] where to store the integral result (real and imaginary
parts)
@param error array[2] where to store the integral error (real and imaginary
parts)
@return 0 on success
@see Integration_type
@see https://github.com/stevengj/cubature
*/
int integral(
    Electrode* sender, Electrode* receiver, _Complex double gamma,
    size_t max_eval, double req_abs_error, double req_rel_error,
    int error_norm, int integration_type, double result[2], double error[2]);

/**
Calculates the internal impedance per unit length of cylindrical conductors.
@param w angular frequency in rad/s
@param rho conductor resistivity
@param radius conductor radius
@param mu relative magnetic permeability of the conductor
@return zin (Ohm/m)
*/
_Complex double internal_impedance(
    double w, double rho, double radius, double mu);

// Longitudinal impedance
/**
Calculates the self longitudinal impedance of a given electrode.
@param electrode
@param w angular frequency in rad/s
@param mu magnetic permeability of the medium
@return zlp
*/
_Complex double longitudinal_self(Electrode* electrode, double w, double mu);

/**
Calculates the mutual longitudinal impedance of given electrodes.
@param sender Electrode
@param receiver Electrode
@param w angular frequency in rad/s
@param mu magnetic permeability of the medium
@param gamma medium propagation constant
@param max_eval specifies a maximum number of function evaluations (0 for no
limit)
@param reqAbsError the absolute error requested (0 to ignore)
@param reqRelError the relative error requested (0 to ignore)
@param error_norm (enumeration defined in cubature.h) error checking scheme
@param integration_type type of integration to be done. If NONE, the value
`result[0] + result[1]*I` will be used and no integration is done.
@param result array[2] where to store the integral result (real and imaginary
parts)
@param error array[2] where to store the integral error (real and imaginary
parts)
@return zlm
@see integral
@see Integration_type
@see https://github.com/stevengj/cubature
*/
_Complex double longitudinal_mutual(
    Electrode* sender, Electrode* receiver, double w, double mu,
    _Complex double gamma, size_t max_eval, double req_abs_error,
    double req_rel_error, int error_norm, int integration_type,
    double result[2], double error[2]);

// Transveral impedance
/**
Calculates the self transversal impedance of a given electrode.
@param electrode
@param kappa medium complex conductivity `(sigma + j*w*eps)` in S/m
@return ztp
*/
_Complex double transversal_self(Electrode* electrode, _Complex double kappa);

/**
Calculates the mutual transversal impedance of given electrodes.
@param sender Electrode
@param receiver Electrode
@param kappa medium complex conductivity `(sigma + j*w*eps)` in S/m
@param gamma medium propagation constant
@param max_eval specifies a maximum number of function evaluations (0 for no
limit)
@param reqAbsError the absolute error requested (0 to ignore)
@param reqRelError the relative error requested (0 to ignore)
@param error_norm (enumeration defined in cubature.h) error checking scheme
@param integration_type type of integration to be done. If NONE, the value
`result[0] + result[1]*I` will be used and no integration is done.
@param result array[2] where to store the integral result (real and imaginary
parts)
@param error array[2] where to store the integral error (real and imaginary
parts)
@return ztm
@see integral
@see Integration_type
@see https://github.com/stevengj/cubature
*/
_Complex double transversal_mutual(
    Electrode* sender, Electrode* receiver, _Complex double kappa,
    _Complex double gamma, size_t max_eval, double req_abs_error,
    double req_rel_error, int error_norm, int integration_type,
    double result[2], double error[2]);

// Electrode system

/**
Calculates the impedance matrices.
@param electrodes array of electrodes
@param num_electrodes number of electrodes
@param zt transversal impedance matrix as a flat array of size
`num_electrodes^2`
@param zl longitudinal impedance matrix as a flat array of size
`num_electrodes^2`
@param w angular frequency in rad/s
@param mu magnetic permeability of the medium
@param kappa medium complex conductivity `(sigma + j*w*eps)` in S/m
@param max_eval specifies a maximum number of function evaluations (0 for no
limit)
@param reqAbsError the absolute error requested (0 to ignore)
@param reqRelError the relative error requested (0 to ignore)
@param error_norm (enumeration defined in cubature.h) error checking scheme
@param integration_type type of integration to be done. If NONE, the value
`result[0] + result[1]*I` will be used and no integration is done.
@param result array[2] where to store the integral result (real and imaginary
parts)
@param error array[2] where to store the integral error (real and imaginary
parts)
@return 0 on sucess
@see integral
@see Integration_type
@see https://github.com/stevengj/cubature
*/
int calculate_impedances(
    Electrode* electrodes, int num_electrodes, _Complex double* zl,
    _Complex double* zt, _Complex double gamma, double w, double mu,
    _Complex double kappa, size_t max_eval, double req_abs_error,
    double req_rel_error, int error_norm, int integration_type);

/**
Add the images effect to the impedance matrices.
@param electrodes array of electrodes
@param images array of electrodes
@param num_electrodes number of electrodes
@param zt transversal impedance matrix as a flat array of size
`num_electrodes^2`
@param zl longitudinal impedance matrix as a flat array of size
`num_electrodes^2`
@param w angular frequency in rad/s
@param mu magnetic permeability of the medium
@param kappa medium complex conductivity `(sigma + j*w*eps)` in S/m
@param ref_l longitudinal current reflection coefficient
@param ref_t transversal current reflection coefficient
@param max_eval specifies a maximum number of function evaluations (0 for no
limit)
@param reqAbsError the absolute error requested (0 to ignore)
@param reqRelError the relative error requested (0 to ignore)
@param error_norm (enumeration defined in cubature.h) error checking scheme
@param integration_type type of integration to be done. If NONE, the value
`result[0] + result[1]*I` will be used and no integration is done.
@param result array[2] where to store the integral result (real and imaginary
parts)
@param error array[2] where to store the integral error (real and imaginary
parts)
@return 0 on sucess
@see integral
@see Integration_type
@see https://github.com/stevengj/cubature
*/
int impedances_images(
    Electrode* electrodes, Electrode* images, int num_electrodes,
    _Complex double* zl, _Complex double* zt, _Complex double gamma, double w,
    double mu, _Complex double kappa, _Complex double ref_l,
    _Complex double ref_t, size_t max_eval, double req_abs_error,
    double req_rel_error, int error_norm, int integration_type);

// Imitance matrix WE building
/**
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
int fill_incidence(
    _Complex double* we, Electrode* electrodes, int num_electrodes,
    double nodes[][3], int num_nodes);

/**
Fills the imitance matrix `we` = `[[ZL/2, -ZL/2, A], [ZT, ZT, B], [C, D, Yn]]`
with the impedance matrices `ZT` and `ZL`, and the nodal admitance Yn.
@param we imitance matrix as flat array of size
`(2*num_electrodes + num_nodes)^2`
@param electrodes array of electrodes
@param num_electrodes number of electrodes
@param num_nodes number of nodes
@param zt transversal impedance matrix as a flat array of size
`num_electrodes^2`
@param zl longitudinal impedance matrix as a flat array of size
`num_electrodes^2`
@param yn nodal admitance matrix as a flat array of size `num_nodes^2`
@return 0 on success
@see fill_incidence
*/
int fill_impedance(
    _Complex double* we, Electrode* electrodes, int num_electrodes,
    int num_nodes, _Complex double* zt, _Complex double* zl,
    _Complex double* yn);

/**
Solves the system of equations that defines the electrode system.
Uses Intel MKL LAPACKE_zgesv.
@param we imitance matrix `[[ZL/2, -ZL/2, A], [ZT, ZT, B], [C, D, Yn]]`. The
LU decomposition is done in-place on we
@param ie RHS vector `[0, 0, ic]^T` where ic are injected currents in each node
and 0 are null vectors of size `num_electrodes` each. The solution replaces
this array in-place
@param num_electrodes number of electrodes
@param num_nodes number of nodes
@return 0 on sucess
@see https://software.intel.com/en-us/mkl-developer-reference-c-gesv
*/
int solve_electrodes(
    _Complex double* we, _Complex double* ie, int num_electrodes, int num_nodes);


#endif /* ELECTRODE_H_ */