/** Funtions to interface the C code to MATLAB.
Uses Matlab's "C Matrix API" to run on Matlab R2017b and
earlier.
*/
#ifndef INTERFACE_MATLAB_H_
#define INTERFACE_MATLAB_H_

#include "mex.h"
#include "electrode.h"
//#include <complex.h>

/** cast_electrode
Cast the matlab Electrode structure to the C struct.
@param matlab_elect pointer to Matlab structure that holds the electrode structure
@param electrode pointer that will be used as the C struct
*/
int
cast_electrode (const mxArray *matlab_elect, mwIndex index, Electrode *electrode);

/*
//auxiliary.h
#define PI 3.1415926535897932384626433832795029L
#define TWO_PI 6.283185307179586
#define FOUR_PI 12.56637061435917
#define MU0 1.256637061435917e-6 //permeability vac.
#define EPS0 8.854187817620e-12 //permittivity vac.
#define RHO_CU 1.689e-8 //copper resistivity
//#define RHO_CU 1.9e-6 //copper resistivity FIXME what is the true value?

int equal_points(const double point_1[3], const double point_2[3]);

double vector_norm(const double start_point[3], const double end_point[3]);

//Electrode.h
enum Integration_type {
    INTG_NONE = 0,
    INTG_DOUBLE,
    INTG_EXP_LOGNF,
    INTG_LOGNF
};

typedef struct {
    double start_point[3];
    double end_point[3];
    double middle_point[3];
    double length;
    double radius;
    _Complex double zi;
} Electrode;

typedef struct {
    const Electrode* sender;
    const Electrode* receiver;
    _Complex double gamma;
    int simplified;
} Integration_data;

int integrand_double(unsigned ndim, const double *t, void *auxdata,
                     unsigned fdim, double *fval);

int exp_logNf(unsigned ndim, const double *t, void *auxdata, unsigned fdim,
              double *fval);

int integral(const Electrode* sender, const Electrode* receiver,
             _Complex double gamma, size_t max_eval, double req_abs_error,
             double req_rel_error, int error_norm, int integration_type,
             double result[2], double error[2]);

int calculate_impedances(const Electrode* electrodes, size_t num_electrodes,
                         _Complex double *zl, _Complex double *zt,
                         _Complex double gamma, _Complex double s, double mur,
                         _Complex double kappa, size_t max_eval,
                         double req_abs_error, double req_rel_error,
                         int error_norm, int integration_type);

int impedances_images(const Electrode* electrodes, const Electrode* images,
                      size_t num_electrodes, _Complex double *zl,
                      _Complex double *zt, _Complex double gamma,
                      _Complex double s, double mur, _Complex double kappa,
                      _Complex double ref_l, _Complex double ref_t,
                      size_t max_eval, double req_abs_error,
                      double req_rel_error, int error_norm, int integration_type);
*/

#endif /* INTERFACE_MATLAB_H_ */
