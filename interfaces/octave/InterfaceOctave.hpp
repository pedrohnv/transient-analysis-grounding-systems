/**
Funtions to interface the C code to GNU Octave using SWIG.
*/
extern "C" {
    #include "electrode.h"
}
#include <complex>

// double array
/** new_doublep
Allocates a `double` array of size `n` and return a pointer to it.
*/
double *
new_doublep (size_t n);

/** delete_doublep
Frees the memory pointed by `arr`.
@return 0 on success
*/
int
delete_doublep (double *arr);

/** doublep_assign
Stores value `val` at position `pos` of array `arr`.
@return 0 on success
*/
int
doublep_assign (double *arr, size_t pos, double val);

/** complexp_value
Return the value stored at position `pos` of array `z`.
*/
double
doublep_value (double *arr, size_t pos);

// complex array
/** new_complexp
Allocates a `complex<double>` array of size `n` and return a pointer to it.
*/
std::complex<double> *
new_complexp (size_t n);

/** delete_complexp
Frees the memory pointed by `z`.
@return 0 on success
*/
int
delete_complexp (std::complex<double> *z);

/** complexp_assign
Stores value `val` at position `pos` of array `z`.
@return 0 on success
*/
int
complexp_assign (std::complex<double> *z, size_t pos, std::complex<double> val);

/** complexp_value
Return the value stored at position `pos` of array `z`.
*/
std::complex<double>
complexp_value (std::complex<double> *z, size_t pos);

// Electrode array
/** new_electrodep
Allocate a `Electrode` array of size `n` and return a pointer to it.
*/
Electrode *
new_electrodep (size_t n);

/** delete_electrodep
Frees the memory pointed by `electrodes`.
@return 0 on success
*/
int
delete_electrodep (Electrode *electrodes);

/** electrodep_assign
Populates struct at position `pos` of array `electrodes`.
@param electrodes array of electrodes
@param pos position in array
@param electrode pointer to an allocated memory
@param start_point array `(x,y,z)` that defines the starting point of the
electrode
@param end_point array `(x,y,z)` defining the ending point
@param radius electrode radius
@param zi total internal impedance of the electrode
@return 0 on success
*/
int
electrodep_assign (Electrode *electrodes, size_t pos, double start_point[3],
                   double end_point[3], double radius, std::complex<double> zi);

/** electrodep_value
Return the struct stored at position `pos` of array `electrodes`.
*/
Electrode
electrodep_value (Electrode *electrodes, size_t pos);

int
Ocalculate_impedances (Electrode *electrodes, size_t num_electrodes,
                       std::complex<double> *zl, std::complex<double> *zt,
                       std::complex<double> gamma, std::complex<double> s,
                       double mur, std::complex<double> kappa, size_t max_eval,
                       double req_abs_error, double req_rel_error, int error_norm,
                       int integration_type);

int
Oimpedances_images (Electrode *electrodes, Electrode *images,
                    size_t num_electrodes, std::complex<double> *zl,
                    std::complex<double> *zt, std::complex<double> gamma,
                    std::complex<double> s, double mur,
                    std::complex<double> kappa, std::complex<double> ref_l,
                    std::complex<double> ref_t, size_t max_eval,
                    double req_abs_error, double req_rel_error,
                   int error_norm, int integration_type);

int
Oharmonic_impedance1 (size_t ns, std::complex<double> *s,
                      std::complex<double> *kappa1, std::complex<double> *kappa2,
                      std::complex<double> *gamma1, Electrode *electrodes,
                      Electrode *images, size_t num_electrodes, double nodes[][3],
                      size_t num_nodes, size_t max_eval, double req_abs_error,
                      double req_rel_error, int error_norm, double rsource,
                      std::complex<double> *zh);

std::complex<double>
sumall (std::complex<double> *z1, size_t n);
