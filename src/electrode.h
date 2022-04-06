/* High Performance implementation of the Hybrid Electromagnetic Model
Released under the General Public License 3 (GPLv3).
All parameters' units are in the SI base if omitted.

Routines to manipulate the Electrode struct and do the base calculations, i.e.,
numerical integration and impedances calculation. It also includes any routines
to build the linear system to be solved that do not depend on linear algebra
libraries.
*/
#ifndef ELECTRODE_H_
#define ELECTRODE_H_

#include <complex.h>
#include <stdlib.h>
#include <stdbool.h>

/** Type of integration and simplification thereof to be done. */
enum
Integration_type {
    /** No integration, returns the geometric distance between middle points */
    INTG_NONE = 1,
    /** Traditional integral along both electrodes
    \f$ \int_0^{L_i} \int_0^{L_k} \frac{e^{-\gamma r}}{r} d\ell_k \, d\ell_i \f$ */
    INTG_DOUBLE = 2,
    /** Traditional integral along sender electrodes
    \f$ \int_0^{L_k} \frac{e^{-\gamma r}}{r} d\ell_k \f$ */
    INTG_SINGLE = 3,
    /** Modified HEM integral
    \f$ \int_0^{L_k} \log_e \left( \frac{R_1 + R_2 + L_i}{R_1 + R_2 - L_i} \right) d\ell_k \f$ */
    INTG_MHEM = 4
};

/** Structure that defines an electrode (conductor segment). */
typedef struct {
    /** Array \f$(x,y,z)_0\f$ that defines the starting point of the electrode */
    float start_point[3];
    /** Array \f$(x,y,z)_1\f$ that defines the ending point of the electrode */
    float end_point[3];
    /** Array \f$\frac{(x,y,z)_1 + (x,y,z)_0}{2}\f$ of the middle point of the electrode */
    float middle_point[3];
    /** The electrode length \f$ \| (x,y,z)_1 - (x,y,z)_0 \|_2 \f$ */
    float length;
    /** The electrode radius */
    float radius;
} Electrode;

/** Structure to make the integration of the "potential" between two electrodes.
@see Integration_type
*/
typedef struct {
    /** Current carrying Electrode that generates the excitation. */
    const Electrode* sender;
    /** Potential receiving Electrode that receives the excitation. */
    const Electrode* receiver;
    /** Complex propagation constant of the medium
    \f$ \gamma = \sqrt{j\omega\mu(\sigma + j\omega\varepsilon)} \f$ */
    _Complex double gamma;
} Integration_data;

/** Structure to pass all needed arguments by magnetic_potential to mag_pot_integral.
@see magnetic_potential
@see https://github.com/stevengj/cubature
*/
typedef struct {
    /** Line integral start point */
    const float* point1;
    /** Line integral end point */
    const float* point2;
    /** array of electrodes */
    const Electrode* electrodes;
    /** number of electrodes \f$ m \f$ */
    size_t num_electrodes;
    /** longitudinal currents array \f$ I_L \f$ */
    const _Complex float* il;
    /** transversal currents array \f$ I_T \f$ */
    const _Complex float* it;
    /** medium propagation constant
        \f$ \gamma = \sqrt{j\omega\mu(\sigma + j\omega\varepsilon)} \f$ */
    _Complex double gamma;
    /** complex frequency
        \f$ s = c + j\omega \f$ */
    _Complex double s;
    /** relative magnetic permeability of the medium \f$ \mu_r \f$ */
    double mur;
    /** medium complex conductivity
        \f$ \sigma + j\omega\varepsilon \f$ in [S/m] */
    _Complex double kappa;
    /** specifies a maximum number of function evaluations (0 for no limit) */
    size_t max_eval;
    /** the absolute error requested (0 to ignore) */
    double req_abs_error;
    /** req_rel_error the relative error requested (0 to ignore) */
    double req_rel_error;
} Field_integrand_data;


/** Populates an Electrode structure.
@param pointer to an allocated memory
@param start_point Array \f$(x,y,z)_0\f$ that defines the starting point of the electrode
@param start_point Array \f$(x,y,z)_0\f$ that defines the ending point of the electrode
@param radius electrode radius The Electrode radius
@return 0 on success
*/
int
populate_electrode (Electrode *electrode, const float start_point[3],
                    const float end_point[3], float radius);

/** Compares two electrodes for equality: have the same radius and coincide the
starts and end points (independent of direction).
@param sender Electrode
@param receiver Electrode
@return true or false
*/
bool
equal_electrodes (const Electrode *sender, const Electrode *receiver);

/** Populates an Electrode array from a CSV file. Each Electrode must be defined in
a single line with parameters: "x0, y0, z0, x1, y1, z1, radius".
@param file_name path to the file
@param electrodes pointer to an array of Electrode to be filled
@param num_electrodes number of electrodes to read (number of lines in file).
    Any line number greater than num_electrodes will be ignored.
@return error == 0 on success \n
        error < 0: number of missing arguments from last line read (as a negative integer) \n
        error > 0: bad input
*/
int
electrodes_file (const char file_name[], Electrode *electrodes, size_t num_electrodes);

/** Fill a nodes array from a CSV file. Each node must be defined in a single line
with parameters: "x, y, z".
@param file_name path to the file
@param nodes array of float[3] to be filled
@param num_nodes number of nodes to read (number of lines in file).
    Any line number greater than num_nodes will be ignored.
@return error == 0 on success \n
        error < 0: number of missing arguments from last line read (as a negative integer) \n
        error > 0: bad input
*/
int
nodes_file (const char file_name[], float *nodes, size_t num_nodes);

/** Segments an electrode populating a passed array of electrodes and an array of nodes.
@param electrodes pointer to an array of Electrode to be filled
@param nodes flat array \f$3 (n + 1)\f$ to be filled by the new nodes
that are created
@param num_segments number \f$ n \f$ of segments to create
@param start_point electrode's start point
@param end_point electrode's end point
@param radius electrode's radius
@return 0 on success
*/
int
segment_electrode (Electrode *electrodes, float *nodes, size_t num_segments,
                   const float *start_point, const float *end_point, float radius);

/** DON'T USE THIS FUNCTIONS, IT'S BUGGED.

Given an array of electrodes, populates the nodes array with all the unique
nodes.
@param nodes flat array of size \f$n \ge 3(m + 1)\f$ to be filled with
   the unique nodes.
@param electrodes array of electrodes
@param num_electrodes number $m$ of electrodes
@returns number of unique nodes \f$n_u\f$; can be used to realloc the nodes array
*/
size_t
nodes_from_elecs (float *nodes, Electrode *electrodes, size_t num_electrodes);

/** Calculates the integrand \f$ \frac{e^{-\gamma r}}{r} \f$ to be integrated
using Cubature.
@param ndim must be = 2
@param t array of size 2 of the electrodes' length as a percentage (0 to 1)
@param auxdata Integration_data with electrodes and propagation constant
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

/** Calculates the integrand \f$ \frac{e^{-\gamma r}}{r} \f$ to be integrated
using Cubature; but the point on the receiver electrode is fixed at its middle.
@param ndim must be = 1
@param t array of size 1 of the sender electrode's length as a percentage (0 to 1)
@param auxdata Integration_data with electrodes and propagation constant
@param fdim must be = 2
@param fval pointer to where the result is stored
@return 0 on success
@see Integration_type
@see Integration_data
@see integral
@see https://github.com/stevengj/cubature
*/
int
integrand_single (unsigned ndim, const double *t, void *auxdata, unsigned fdim,
                  double *fval);

/** Calculates the mHEM integrand
\f$ \int_0^{L_k} \log_e \left( \frac{R_1 + R_2 + L_i}{R_1 + R_2 - L_i} \right) d\ell_k \f$
to be used in Cubature integration.
@param ndim must be = 1
@param t array of size 1 with receiver's length percentage (0 to 1)
@param auxdata Integration_data with electrodes and propagation constant
@param fdim must be = 2
@param fval pointer to where the result is stored
@return 0 on success
@see Integration_type
@see Integration_data
@see integral
@see https://github.com/stevengj/cubature
*/
int
logNf (unsigned ndim, const double *t, void *auxdata, unsigned fdim, double *fval);

/** Calculates the mHEM closed form solution of the integral when the sender and receiver
segments are the same (self impedance).
\f$ 2 L_k \left[ \log_e \left( \frac{\sqrt{1 + (b/L_k)^2} + 1}{b/L_k} \right)
- \sqrt{1 + \left(\frac{b}{L_k}\right)^2} + \frac{b}{L_k}  \right] \f$
*/
_Complex double
self_integral (const Electrode *sender);

/** Calculates the integral along the sender and receiver using Cubature.
Note that the integration stops when either max_eval or req_abs_error
or req_rel_error is reached.
@param sender Current carrying Electrode that generates the excitation.
@param receiver Potential receiving Electrode that receives the excitation.
@param gamma medium propagation constant
  \f$ \gamma = \sqrt{j\omega\mu(\sigma + j\omega\varepsilon)} \f$
@param max_eval specifies a maximum number of function evaluations (0 for no limit)
@param req_abs_error the absolute error requested (0 to ignore)
@param req_rel_error the relative error requested (0 to ignore)
@param integration_type type of integration to be done.
@param result array[2] where to store the integral result (real and imaginary parts)
@param error array[2] where to store the integral error (real and imaginary parts)
@return 0 on success \n
      -10 if integration_type is unrecognized.
@see Integration_type
@see https://github.com/stevengj/cubature
*/
int
integral (const Electrode *sender, const Electrode *receiver, _Complex double gamma,
          size_t max_eval, double req_abs_error, double req_rel_error,
          int integration_type, double result[2], double error[2]);

/** Calculates the impedance matrices \f$ Z_L, Z_T \f$. As they are symmetric,
only their lower half is stored (set).  If pointer zl = zt, then the resulting
matrix will be filled with zt.\n
If integration_type == INTG_MHEM or integration_type == INTG_NONE, then parameters
gamma, s, mur and kappa are ignored such that
\f$ \frac{j\omega\mu}{4\pi} = \frac{1}{4\pi\,(\sigma + j\omega\varepsilon)} = 1 \f$.
@param zl longitudinal impedance matrix \f$ Z_L \f$ as a flat array of size \f$ m^2 \f$
@param zt transversal impedance matrix \f$ Z_T \f$ as a flat array of size \f$ m^2 \f$
@param electrodes array of electrodes
@param num_electrodes number of electrodes \f$ m \f$
@param gamma medium propagation constant
    \f$ \gamma = \sqrt{j\omega\mu(\sigma + j\omega\varepsilon)} \f$
@param s complex angular frequency \f$ s = c + j\omega \f$ in [rad/s]
@param mur relative magnetic permeability of the medium \f$ \mu_r \f$
@param kappa medium complex conductivity \f$ \sigma + j\omega\varepsilon \f$ in [S/m]
@param max_eval specifies a maximum number of function evaluations (0 for no limit)
@param req_abs_error the absolute error requested (0 to ignore)
@param req_rel_error the relative error requested (0 to ignore)
@param integration_type type of integration to be done.
@return 0 on success
@see integral
@see Integration_type
@see https://github.com/stevengj/cubature
*/
int
calculate_impedances (_Complex float *zl, _Complex float *zt,
                      const Electrode *electrodes, size_t num_electrodes,
                      _Complex double gamma, _Complex double s, double mur,
                      _Complex double kappa, size_t max_eval, double req_abs_error,
                      double req_rel_error, int integration_type);

/** Add the images' effect to the impedance matrices \f$ Z_L \f$ and \f$ Z_T \f$.
As they are symmetric, only their lower half is stored (set). If pointer
zl = zt, then the resulting matrix will be filled with zt.\n
If integration_type == INTG_MHEM or integration_type == INTG_NONE, then parameters
gamma, s, mur, kappa, ref_l and ref_t are ignored such that
\f$ \Gamma_L \frac{j\omega\mu}{4\pi} = \Gamma_T \frac{1}{4\pi\,(\sigma + j\omega\varepsilon)} = 1 \f$.
@param zl longitudinal impedance matrix \f$ Z_L \f$ as a flat array of size \f$ m^2 \f$
@param zt transversal impedance matrix \f$ Z_T \f$ as a flat array of size \f$ m^2 \f$
@param electrodes array of electrodes
@param images array of electrodes
@param num_electrodes number of electrodes \f$ m \f$
@param gamma medium propagation constant
    \f$ \gamma = \sqrt{j\omega\mu(\sigma + j\omega\varepsilon)} \f$
@param s complex angular frequency \f$ s = c + j\omega \f$ in [rad/s]
@param mur relative magnetic permeability of the medium \f$ \mu_r \f$
@param kappa medium complex conductivity \f$ \sigma + j\omega\varepsilon \f$ in [S/m]
@param ref_l longitudinal current reflection coefficient \f$ \Gamma_L \f$
@param ref_t transversal current reflection coefficient \f$ \Gamma_T \f$
@param max_eval specifies a maximum number of function evaluations (0 for no limit)
@param req_abs_error the absolute error requested (0 to ignore)
@param req_rel_error the relative error requested (0 to ignore)
@param integration_type type of integration to be done.
@return 0 on success
@see integral
@see Integration_type
@see https://github.com/stevengj/cubature
*/
int
impedances_images (_Complex float *zl, _Complex float *zt,
                   const Electrode *electrodes, const Electrode *images,
                   size_t num_electrodes, _Complex double gamma,
                   _Complex double s, float mur, _Complex double kappa,
                   _Complex double ref_l, _Complex double ref_t,
                   size_t max_eval, double req_abs_error, double req_rel_error,
                   int integration_type);

/** Calculates the scalar electric potential \f$ u \f$ to remote earth at a point.
@param point array \f$(x, y, z)\f$
@param electrodes array of electrodes
@param num_electrodes number of electrodes \f$ m \f$
@param it transversal currents array \f$ I_T \f$
@param gamma medium propagation constant
   \f$ \gamma = \sqrt{j\omega\mu(\sigma + j\omega\varepsilon)} \f$
@param kappa medium complex conductivity \f$ \sigma + j\omega\varepsilon \f$ in [S/m]
@param max_eval specifies a maximum number of function evaluations (0 for no limit)
@param req_abs_error the absolute error requested (0 to ignore)
@param req_rel_error the relative error requested (0 to ignore)
@return potential \f$ u \f$ to remote earth
@see https://github.com/stevengj/cubature
*/
_Complex float
electric_potential (const float *point, const Electrode *electrodes,
                    size_t num_electrodes, const _Complex float *it,
                    _Complex double gamma, _Complex double kappa,
                    size_t max_eval, double req_abs_error, double req_rel_error);

/** Calculates the magnetic vector potential \f$ \vec A \f$.
@param point array \f$(x, y, z)\f$
@param electrodes array of electrodes
@param num_electrodes number of electrodes \f$ m \f$
@param il longitudinal currents array \f$ I_L \f$
@param gamma medium propagation constant
    \f$ \gamma = \sqrt{j\omega\mu(\sigma + j\omega\varepsilon)} \f$
@param mur relative magnetic permeability of the medium \f$ \mu_r \f$
@param max_eval specifies a maximum number of function evaluations (0 for no limit)
@param req_abs_error the absolute error requested (0 to ignore)
@param req_rel_error the relative error requested (0 to ignore)
@param va pointer to array of size 3 in which the magnetic vector potential
components \f$ (A_x,A_y,A_z) \f$ will be stored on return.
@return 0 on success
@see https://github.com/stevengj/cubature
*/
int
magnetic_potential (const float *point, const Electrode *electrodes,
                    size_t num_electrodes, const _Complex float *il,
                    _Complex double gamma, double mur, size_t max_eval,
                    double req_abs_error, double req_rel_error,
                    _Complex float *va);

/** Integrand to calculate the electric field caused by a differential
transversal current \f$ dI_T \f$.
@param ndim should be 1
@param t integration variable (electrode percentage 0 to 1)
@param auxdata pointer to Field_integrand_data structure
@param fdim should be 6 (Real and Imag.):
    \f$ (\Re(E_x), \Im(E_x), \Re(E_y), \Im(E_y), \Re(E_z), \Im(E_z)) \f$
@param fval where the results are stored
    \f$ (\Re(E_x), \Im(E_x), \Re(E_y), \Im(E_y), \Re(E_z), \Im(E_z)) \f$
return 0 on success
@see https://github.com/stevengj/cubature
@see electric_field
*/
int
elec_field_integrand (unsigned ndim, const double *t, void *auxdata,
                      unsigned fdim, double *fval);

/** Calculates the electric field at a point. \f$ \vec E = -\nabla u - s \vec A \f$
@param point array \f$ (x, y, z) \f$
@param electrodes array of electrodes
@param num_electrodes number of electrodes \f$ m \f$
@param il longitudinal currents array \f$ I_L \f$
@param it transversal currents array \f$ I_T \f$
@param gamma medium propagation constant
    \f$ \gamma = \sqrt{j\omega\mu(\sigma + j\omega\varepsilon)} \f$
@param s complex angular frequency \f$ s = c + j\omega \f$ in [rad/s]
@param mur relative magnetic permeability of the medium \f$ \mu_r \f$
@param kappa medium complex conductivity \f$ \sigma + j\omega\varepsilon \f$ in [S/m]
@param max_eval specifies a maximum number of function evaluations (0 for no limit)
@param req_abs_error the absolute error requested (0 to ignore)
@param req_rel_error the relative error requested (0 to ignore)
@param ve pointer to array of size 3 in which the electric field
components \f$ (E_x, E_y, E_z) \f$ will be add to on return.
@return 0 on success
@see elec_field_integrand
*/
int
electric_field (const float *point, const Electrode *electrodes,
                size_t num_electrodes, const _Complex float *il,
                const _Complex float *it, _Complex double gamma,
                _Complex double s, double mur, _Complex double kappa,
                size_t max_eval, double req_abs_error, double req_rel_error,
                _Complex float *ve);

/** Interface to pass magnetic_potential as an integrand to Cubature when
calculating the voltage along a path.
@param ndim should be 1
@param t integration variable (electrode percentage 0 to 1)
@param auxdata pointer to Field_integrand_data structure
@param fdim should be 2: \f$ \Re(\vec A \cdot d\vec\ell) + \Im(\vec A \cdot d\vec\ell) \f$
@param fval where the results are stored (array of size 2)
return 0 on success
@see https://github.com/stevengj/cubature
*/
int
v_mag_pot_integrand (unsigned ndim, const double *t, void *auxdata, unsigned fdim,
                     double *fval);

/** Calculates the voltage \f$ \Delta U_{12} \f$ along a straight line from
\f$ \vec p_1 \f$ to \f$ \vec p_2 \f$.\n
\f$ \Delta U_{12} = \int_{\vec p_1}^{\vec p_2}
\left( -\nabla u - s \vec A \right) d\vec\ell \f$
@param point1 Line integral start point \f$ \vec p_1 = (x,y,z)_1 \f$
@param point2 Line integral start point \f$ \vec p_2 = (x,y,z)_2 \f$
@param electrodes array of electrodes
@param num_electrodes number of electrodes \f$ m \f$
@param il longitudinal currents array \f$ I_L \f$
@param gamma medium propagation constant
    \f$ \gamma = \sqrt{j\omega\mu(\sigma + j\omega\varepsilon)} \f$
@param s complex angular frequency \f$ s = c + j\omega \f$ in [rad/s]
@param mur relative magnetic permeability of the medium \f$ \mu_r \f$
@param kappa medium complex conductivity \f$ \sigma + j\omega\varepsilon \f$ in [S/m]
@param max_eval specifies a maximum number of function evaluations (0 for no limit)
@param req_abs_error the absolute error requested (0 to ignore)
@param req_rel_error the relative error requested (0 to ignore)
@return the voltage along the line
@see https://github.com/stevengj/cubature
*/
_Complex float
voltage (const float *point1, const float *point2,
         const Electrode *electrodes, size_t num_electrodes,
         const _Complex float *il, const _Complex float *it,
         _Complex double gamma, _Complex double s, double mur,
         _Complex double kappa, size_t max_eval, double req_abs_error,
         double req_rel_error);

#endif /* ELECTRODE_H_ */
