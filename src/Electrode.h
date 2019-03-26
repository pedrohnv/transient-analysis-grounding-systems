/** High Performance Hybrid Electromagnetic Model calculations in C.

All parameters' units are in the SI base units if omitted.

Routines to manipulate the Electrode struct and do the base calculations, i.e.,
numerical integration and impedances calculation. It also includes any routines
to build the linear system to be solved that do not depend on linear algebra
libraries.

TODO insert condition to check if sender == receiver?
    during integration and calculate distance from center radius
*/
#ifndef ELECTRODE_H_
#define ELECTRODE_H_

#include <complex.h>
#include <stdlib.h>

#define MKL_Complex16 _Complex double //overwrite type

//default integration options ===============================================
/** Integration_type
Type of integration and simplification thereof to be done.
@param NONE integration has closed form solution \f e^{-\gamma \bar r}}{\bar r} L_s L_r \f
@param INTG_DOUBLE \f$ \int_0^{L_s} \int_0^{L_r} \frac{e^{-\gamma r}}{r} dl_r dl_s \f$
@param INTG_EXP_LOGNF \f$ \int_0^{L_r} e^{-\gamma \bar r} Log(N_f) dl_r \f$
@param INTG_LOGNF \f$ e^{-\gamma \bar r} \int_0^{L_r} Log(N_f) dl_r \f$
*/
enum
Integration_type {
    INTG_NONE = 1,
    INTG_DOUBLE,
    INTG_EXP_LOGNF,
    INTG_LOGNF
};

/** Electrode
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
typedef struct {
    double start_point[3];
    double end_point[3];
    double middle_point[3];
    double length;
    double radius;
    _Complex double zi;
} Electrode;

/** Integration_data
Structure to make the integration ot the "potential" between two electrodes.
@param sender electrode that generates the excitation
@param receiver electrode that is excitated
@param gamma complex propagation constant of the medium
@param simplified integration being done is INTG_EXP_LOGNF (false) or
INTG_LOGNF (true); no use otherwise
@see Integration_type
*/
typedef struct {
    const Electrode *sender;
    const Electrode *receiver;
    _Complex double gamma;
    int simplified;
} Integration_data;

/** populate_electrode
Populates an Electrode structure.
@param electrode pointer to an allocated memory
@param start_point array `(x,y,z)` that defines the starting point of the
electrode
@param end_point array `(x,y,z)` defining the ending point
@param radius electrode radius
@param zi total internal impedance of the electrode
@return 0 on success
*/
int
populate_electrode (Electrode *electrode, const double start_point[3],
                    const double end_point[3], double radius, _Complex double zi);

/** electrodes_file
Populates an Electrode array from a file. Each Electrode must be defined in
a single line with parameters: "x0 y0 z0 x1 y1 z1 radius Re(zi) Im(zi)".
The @param `num_electrodes` is the number of lines in the file to be read.
Any line number greater than `num_electrodes` will be ignored.
@param file_name path to the file
@param electrodes pointer to an array of Electrode to be filled
@param num_electrodes number of electrodes to read (number of lines in file)
@return error
    error == 0 on success
    error < 0: number of missing arguments from last line read (as a negative integer)
    error > 0: bad input
*/
int
electrodes_file (const char file_name[], Electrode *electrodes,
                 size_t num_electrodes);

/** nodes_file
Fill a nodes array from a file. Each node must be defined in a single line
with parameters: "x y z".
The @param `num_nodes` is the number of lines in the file to be read.
Any line number greater than `num_nodes` will be ignored.
@param file_name path to the file
@param nodes array of double[3] to be filled
@param num_nodes number of nodes to read (number of lines in file)
@return number of missing arguments from last line read (as a negative integer)
*/
int
nodes_file (const char file_name[], double nodes[][3], size_t num_nodes);

/** segment_electrode
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
int
segment_electrode (Electrode *electrodes, double nodes[][3], size_t num_segments,
                   const double *start_point, const double *end_point,
                   double radius, _Complex double unit_zi);

/** integrand_double
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
int
integrand_double (unsigned ndim, const double *t, void *auxdata, unsigned fdim,
                  double *fval);

/** exp_logNf
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
int
exp_logNf (unsigned ndim, const double *t, void *auxdata, unsigned fdim,
           double *fval);

/** integral
Calculates the integral along the sender and receiver using Cubature.
@param sender Electrode
@param receiver Electrode
@param gamma medium propagation constant
@param max_eval specifies a maximum number of function evaluations (0 for no
limit)
@param reqAbsError the absolute error requested (0 to ignore)
@param reqRelError the relative error requested (0 to ignore)
@param error_norm (enumeration defined in cubature.h) error checking scheme
@param integration_type type of integration to be done.
@param result array[2] where to store the integral result (real and imaginary
parts)
@param error array[2] where to store the integral error (real and imaginary
parts)
@return 0 on success
@see Integration_type
@see https://github.com/stevengj/cubature
*/
int
integral (const Electrode *sender, const Electrode *receiver, _Complex double gamma,
          size_t max_eval, double req_abs_error, double req_rel_error,
          int error_norm, int integration_type, double result[2], double error[2]);

/** internal_impedance
Calculates the internal impedance per unit length of cylindrical conductors.
@param s complex angular frequency `c + I*w` (rad/s)
@param rho conductor resistivity
@param radius conductor radius
@param mur relative magnetic permeability of the conductor
@return zin (Ohm/m)
*/
_Complex double
internal_impedance (_Complex double s, double rho, double radius, double mur);

// Longitudinal impedance
/** longitudinal_self
Calculates the self longitudinal impedance of a given electrode.
@param electrode
@param s complex angular frequency `c + I*w` (rad/s)
@param mur relative magnetic permeability of the medium
@return zlp
*/
_Complex double
longitudinal_self (const Electrode *electrode, _Complex double s, double mur);

/** longitudinal_mutual
Calculates the mutual longitudinal impedance of given electrodes.
@param sender Electrode
@param receiver Electrode
@param s complex angular frequency `c + I*w` (rad/s)
@param mur relative magnetic permeability of the medium
@param gamma medium propagation constant
@param max_eval specifies a maximum number of function evaluations (0 for no
limit)
@param reqAbsError the absolute error requested (0 to ignore)
@param reqRelError the relative error requested (0 to ignore)
@param error_norm (enumeration defined in cubature.h) error checking scheme
@param integration_type type of integration to be done.
@param result array[2] where to store the integral result (real and imaginary
parts)
@param error array[2] where to store the integral error (real and imaginary
parts)
@return zlm
@see integral
@see Integration_type
@see https://github.com/stevengj/cubature
*/
_Complex double
longitudinal_mutual (const Electrode *sender, const Electrode *receiver,
                     _Complex double s, double mur, _Complex double gamma,
                     size_t max_eval, double req_abs_error, double req_rel_error,
                     int error_norm, int integration_type, double result[2],
                     double error[2]);

// Transveral impedance
/** transversal_self
Calculates the self transversal impedance of a given electrode.
@param electrode
@param kappa medium complex conductivity `(sigma + I*w*eps)` in S/m
@return ztp
*/
_Complex double
transversal_self (const Electrode *electrode, _Complex double kappa);

/** transversal_mutual
Calculates the mutual transversal impedance of given electrodes.
@param sender Electrode
@param receiver Electrode
@param kappa medium complex conductivity `(sigma + I*w*eps)` in S/m
@param gamma medium propagation constant
@param max_eval specifies a maximum number of function evaluations (0 for no
limit)
@param reqAbsError the absolute error requested (0 to ignore)
@param reqRelError the relative error requested (0 to ignore)
@param error_norm (enumeration defined in cubature.h) error checking scheme
@param integration_type type of integration to be done.
@param result array[2] where to store the integral result (real and imaginary
parts)
@param error array[2] where to store the integral error (real and imaginary
parts)
@return ztm
@see integral
@see Integration_type
@see https://github.com/stevengj/cubature
*/
_Complex double
transversal_mutual (const Electrode *sender, const Electrode *receiver,
                    _Complex double kappa, _Complex double gamma, size_t max_eval,
                    double req_abs_error, double req_rel_error, int error_norm,
                    int integration_type, double result[2], double error[2]);

// Electrode system

/** calculate_impedances
Calculates the impedance matrices.
@param electrodes array of electrodes
@param num_electrodes number of electrodes
@param zl longitudinal impedance matrix as a flat array of size
`num_electrodes^2`
@param zt transversal impedance matrix as a flat array of size
`num_electrodes^2`
@param gamma medium propagation constant
@param s complex angular frequency `c + I*w` (rad/s)
@param mur relative magnetic permeability of the medium
@param kappa medium complex conductivity `(sigma + I*w*eps)` in S/m
@param max_eval specifies a maximum number of function evaluations (0 for no
limit)
@param reqAbsError the absolute error requested (0 to ignore)
@param reqRelError the relative error requested (0 to ignore)
@param error_norm (enumeration defined in cubature.h) error checking scheme
@param integration_type type of integration to be done.
@return 0 on sucess
@see integral
@see Integration_type
@see https://github.com/stevengj/cubature
*/
int
calculate_impedances (const Electrode *electrodes, size_t num_electrodes,
                      _Complex double *zl, _Complex double *zt,
                      _Complex double gamma, _Complex double s, double mur,
                      _Complex double kappa, size_t max_eval, double req_abs_error,
                      double req_rel_error, int error_norm, int integration_type);

/** impedances_images
Add the images effect to the impedance matrices.
@param electrodes array of electrodes
@param images array of electrodes
@param num_electrodes number of electrodes
@param zl longitudinal impedance matrix as a flat array of size
`num_electrodes^2`s
@param zt transversal impedance matrix as a flat array of size
`num_electrodes^2`
@param gamma medium propagation constant
@param s complex angular frequency `c + I*w` (rad/s)
@param mur relative magnetic permeability of the medium
@param kappa medium complex conductivity `(sigma + I*w*eps)` in S/m
@param ref_l longitudinal current reflection coefficient
@param ref_t transversal current reflection coefficient
@param max_eval specifies a maximum number of function evaluations (0 for no
limit)
@param reqAbsError the absolute error requested (0 to ignore)
@param reqRelError the relative error requested (0 to ignore)
@param error_norm (enumeration defined in cubature.h) error checking scheme
@param integration_type type of integration to be done.
@return 0 on sucess
@see integral
@see Integration_type
@see https://github.com/stevengj/cubature
*/
int
impedances_images (const Electrode *electrodes, const Electrode *images,
                   size_t num_electrodes, _Complex double *zl, _Complex double *zt,
                   _Complex double gamma, _Complex double s, double mur,
                   _Complex double kappa, _Complex double ref_l,
                   _Complex double ref_t, size_t max_eval, double req_abs_error,
                   double req_rel_error, int error_norm, int integration_type);

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

#endif /* ELECTRODE_H_ */
