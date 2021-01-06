/** High Performance Hybrid Electromagnetic Model calculations in C.

All parameters' units are in the SI base if omitted.

Constants, auxiliary functions and routines.
*/
#ifndef AUXILIARY_H_
#define AUXILIARY_H_

#include <complex.h>
#include <stdio.h>
#include <stdbool.h>
#include <fftw3.h>

//constants
/** math constant \f$ \pi \f$ */
#define PI 3.1415926535897932384626433832795029L
/** math constant \f$ 2\pi \f$ */
#define TWO_PI (2.0 * PI)
/** math constant \f$ 4\pi \f$ */
#define FOUR_PI (4.0 * PI)
/** Vacuum permeability \f$ \mu_0 \f$ */
#define MU0 1.256637061435917e-6
/** Vacuum permittivity \f$ \varepsilon_0 \f$ */
#define EPS0 8.854187817620e-12
/** Copper resistivity \f$ \rho_{cu} \f$ */
#define RHO_CU 1.689e-8


/** Fill array \f$ u \f$ with \f$ n \f$ linearly spaced numbers between \f$ a \f$
and \f$ b \f$ (included).
@param a first number
@param b last number
@param n size of array
@param u array to be filled
@return pointer to filled array
@see Source code taken from https://github.com/ntessore/algo
*/
double*
linspace (double a, double b, size_t n, double u[]);

/** Fill array \f$ u \f$ with \f$ n \f$ logarithmically spaced numbers between
\f$ 10^a \f$ and \f$ 10^b \f$ (included).
@param a first power
@param b last power
@param n size of array
@param u array to be filled
@return pointer to filled array
@see Source code taken from https://github.com/ntessore/algo
*/
double*
logspace (double a, double b, size_t n, double u[]);

/** Computes the harmonic electromagnetic wave length \f$ \lambda \f$ in a lossy medium.
@param f frequency (Hz)
@param sigma medium (real) conductivity \f$ \sigma \f$ in \f$ \Omega \cdot m \f$
@param ep medium electric permittivity \f$ \varepsilon = \varepsilon_r \, \varepsilon_0\f$
@param mur medium relative magnetic permeability \f$ \mu_r \f$
@return \f$ \lambda \f$ in \f$ m \f$
*/
double
wave_length (double f, double sigma, double ep, double mur);

/** Calculates the soil parameters \f$ \sigma(s) \f$ and \f$ \varepsilon_r(s) \f$
based on the Alipio-Visacro model [1].

\f$ \sigma = \sigma_0 + \sigma_0 \times h(\sigma_0)\left( \frac{s}{\text{1 MHz}} \right)^g \f$

\f$
\varepsilon_r = \frac{\varepsilon_\infty'}{\varepsilon_0} +
\frac{\tan(\pi g / 2) \times 10^{-3}}{2\pi\varepsilon_0(\text{1 MHz})^g}
\sigma_0 \times h(\sigma_0) s^{g - 1}
\f$

Recommended values of \f$ h(\sigma_0) \f$, \f$ g \f$ and
\f$ \varepsilon_\infty' / \varepsilon_0  \f$ are given in Fig. 8 of [1]:

<table>
<tr><th>Results                            <th>\f$ h(\sigma_0) \f$                  <th>\f$ g \f$  <th>\f$ \varepsilon_\infty' / \varepsilon_0 \f$
<tr><td rowspan="1">mean                   <td> \f$ 1.26 \times (1000 \sigma_0)^{−0.73} \f$<td> 0.54      <td> 12
<tr><td rowspan="1">relatively conservative<td> \f$ 0.95 \times (1000 \sigma_0)^{−0.73} \f$<td> 0.58      <td> 8
<tr><td rowspan="1">conservative           <td> \f$ 0.70 \times (1000 \sigma_0)^{−0.73} \f$<td> 0.62      <td> 4
</table>

[1] R. Alipio and S. Visacro, "Modeling the Frequency Dependence of Electrical
Parameters of Soil," in IEEE Transactions on Electromagnetic Compatibility,
vol. 56, no. 5, pp. 1163-1171, Oct. 2014, doi: 10.1109/TEMC.2014.2313977.
@param sigma pointer to where the conductivity \$f \sigma(s) \$f is written in S/m
@param epsr pointer to where the relative permitivitty \f$ \varepsilon_r(s) \f$ is written
@param sigma0 value of the soil conductivity in low frequency \f$ \sigma_0 \f$ in S/m
@param s complex frequency \f$ s = c + j\omega \f$ of interest in rad/s
@param h parameter \f$ h(\sigma_0) \f$
@param g parameter \f$ g \f$
@param eps_ratio parameter \f$ \varepsilon_\infty' / \varepsilon_0 \f$
@return 0 on success
*/
int
alipio_soil (_Complex double* sigma, _Complex double* epsr, double sigma0,
	           _Complex double s, double h, double g, double eps_ratio);

/* Calculates the soil parameters \f$ \sigma(s) \f$ and \f$ \varepsilon_r(s) \f$
based on the Smith-Longmire model as presented in [1].

[1] D. Cavka, N. Mora, F. Rachidi, A comparison of frequency-dependent soil
models: application to the analysis of grounding systems, IEEE Trans.
Electromagn. Compat. 56 (February (1)) (2014) 177–187.

@param sigma pointer to where the conductivity \$f \sigma(s) \$f is written in S/m
@param epsr pointer to where the relative permitivitty \f$ \varepsilon_r(s) \f$ is written
@param sigma0 value of the soil conductivity in low frequency \f$ \sigma_0 \f$ in S/m
@param s complex frequency \f$ s = c + j\omega \f$ of interest in rad/s
@param erinf parameter \f$ \varepsilon_\infty'\f$
@return 0 on success
*/
int
smith_longmire_soil (_Complex double* sigma, _Complex double* epsr, double sigma0,
                     _Complex double s, double erinf);

/** Check if two points are the same using DBL_EPSILON as tolerance.
@param point1 first point in \f$ \mathbf{R}^3 \f$
@param point2 second point in \f$ \mathbf{R}^3 \f$
@return true or false
@see equal_points_tol
*/
bool
equal_points (const double* point_1, const double* point_2);

/** Check if two points are the same within a tolerance.
@param point1 first point in \f$ \mathbf{R}^3 \f$
@param point2 second point in \f$ \mathbf{R}^3 \f$
@param tol the tolerance \f$ \epsilon \f$ applied to each coordinate.
@return false if the difference between any coordinates is greater than tol
(e.g. \f$ |x_1 - x_2| > \epsilon \f$).
@see equal_points
*/
bool
equal_points_tol (const double* point_1, const double* point_2, double tol);

/** Computes the \f$ \|\cdot\|_2 \f$ (euclidian length) of the vector which starts at
start_point and ends at end_point.
@param start_point array \f$(x,y,z)_1\f$ that defines the starting point
@param end_point array \f$(x,y,z)_0\f$ that defines the ending point
@return \f$ \| (x,y,z)_1 - (x,y,z)_0 \|_2 \f$
*/
double
vector_length (const double start_point[3], const double end_point[3]);

/** LAPACK routine to compute the Euclidean norm of a vector \f$ \|x\|_2 \f$.
@param n Specifies the number of elements in vector \f$x\f$.
@param x Array of size at least \f$(1 + (n-1) \times i)\f$
@param incx Specifies the increment \f$ i \f$ for the elements of \f$x\f$
@return \f$ \|x\|_2 \f$
*/
extern double
snrm2_ (int* n, double* x, int* incx);

/** Prints a complex (double) matrix to a file.
@param m number of rows
@param n number of columns
@param a pointer to array where matrix is stored
@param lda leading dimension of \f$ a \f$
@param fp pointer to file stream
@return 0 on success
*/
int
complex_matrix_file (size_t m, size_t n, const _Complex double* a, int lda, FILE* fp);

/** Prints a real (double) matrix to a file.
@param m number of rows
@param n number of columns
@param a pointer to array where matrix is stored
@param lda leading dimension of \f$ a \f$
@param fp pointer to file stream
@return 0 on success
*/
int
double_matrix_file (size_t m, size_t n, const double* a, int lda, FILE* fp);

/** FORTRAN subroutine (from SLATEC) to calculate I-Bessel function, i.e.,
modified Bessel function of the first kind, with complex argument.
    \f$ cy = I_{fnu} (z) \f$
@param zr argument's real part \f$\Re(z)\f$
@param zi argument's imaginary part \f$\Im(z)\f$
@param fnu order of initial I function
@param kode scaling \n
        1: no scaling \n
        2: \f$e^{-|x|}\f$
@param n number of members of the sequence
@param cyr result's real part \f$\Re(cy)\f$
@param cyi result's imaginary part \f$\Im(cy)\f$
@param nz number of components set to zero due to underflow
@param ierr error flag
@see http://netlib.org/amos/zbesi.f
*/
extern int
zbesi_ (double* zr, double* zi, double* fnu, int* kode, int* n, double* cyr,
        double* cyi, int* nz, int* ierr);

/** Prints a COLUMN MAJOR matrix to stdio.
@param desc description to print before the matrix
@param m number of rows
@param n number of columns
@param a column major matrix to be printed
@param lda leading dimension of a
*/
void
print_matrix (char *desc, int m, int n, const _Complex double* a, int lda);

/** Prints a ROW MAJOR matrix to stdio.
@param desc description to print before the matrix
@param m number of rows
@param n number of columns
@param a row major matrix to be printed
@param lda leading dimension of a
*/
void
print_matrix_row (char *desc, int m, int n, const _Complex double* a, int lda);

/** Copies the source matrix into target matrix considering COLUMN MAJOR storage.
\f$ B := A \f$
@param source pointer to be copied
@param target pointer where the copy will be stored
@param lds leading dimension of the source array
@param ldt leading dimension of the target array
@param nline number of lines in source array to copy
@param ncol number of columns in source array to copy
@return 0 on success
*/
int
matrix_copy (const _Complex double* source, _Complex double* target,
			       size_t lds, size_t ldt, size_t nline, size_t ncol);

/** Copies the transpose of source matrix into target matrix considering
COLUMN MAJOR storage.
\f$ B := A^{T} \f$
@param source pointer to be copied
@param target pointer where the copy will be stored
@param lds leading dimension of the source array
@param ldt leading dimension of the target array
@param nline number of lines in source array to copy
@param ncol number of columns in source array to copy
@return 0 on success
*/
int
transpose_copy (const _Complex double* source, _Complex double* target,
                size_t lds, size_t ldt, size_t nline, size_t ncol);

/** pc_copy
Copies the column permutation of `source` matrix into `target` matrix
considering COLUMN MAJOR storage.
@param source pointer to be copied
@param target pointer where the copy will be stored
@param lds leading dimension of the source array
@param ldt leading dimension of the target array
@param nline number of lines in source array to copy
@param ncol number of columns in source array to copy
@return 0 on success
*/
int
pc_copy (const _Complex double* source, _Complex double* target,
         size_t lds, size_t ldt, size_t nline, size_t ncol);

/** Copies the line permutation of source matrix into target matrix
considering COLUMN MAJOR storage.
\f$ B := A^{Pl} \f$
@param source pointer to be copied
@param target pointer where the copy will be stored
@param lds leading dimension of the source array
@param ldt leading dimension of the target array
@param nline number of lines in source array to copy
@param ncol number of columns in source array to copy
@return 0 on success
*/
int
pl_copy (const _Complex double* source, _Complex double* target,
         size_t lds, size_t ldt, size_t nline, size_t ncol);

/** Copies the column and line permutation of source matrix into target matrix
considering COLUMN MAJOR storage.
\f$ B := A^{Pcl} \f$
@param source pointer to be copied
@param target pointer where the copy will be stored
@param lds leading dimension of the source array
@param ldt leading dimension of the target array
@param nline number of lines in source array to copy
@param ncol number of columns in source array to copy
@return 0 on success
*/
int
pcl_copy (const _Complex double* source, _Complex double* target,
          size_t lds, size_t ldt, size_t nline, size_t ncol);


/** Numerical Laplace Transform \f$ g(s) = \mathcal{L}\{f(t)\} \f$ of a real
input. Uses FFTW, therefore is not thread-safe! Based on:
GÓMEZ, Pablo; URIBE, Felipe A. The numerical Laplace transform: An accurate
technique for analyzing electromagnetic transients on power system devices.
International Journal of Electrical Power & Energy Systems, v. 31, n. 2-3,
p. 116-123, 2009.

It is used the value \f$ c=\log(n^2)/T \f$. That value affects the aliasing.
If you want a different value, change it inside the source code.

@param f sampled function \f$ f(t) \f$ of size \f$ n \f$ from
\f$ t = 0 \f$ to \f$ t = t_\max \f$ with uniform spacing \f$ \Delta t \f$
@param t time stamps \f$ t \f$ of the sampling of size \f$ n \f$, assumed uniformly spaced
@param g transformed function \f$ g(s) \f$ of size \f$ \left \lfloor{n/2}\right \rfloor + 1 \f$
@param s ``frequencies'' \f$ s = c + j\omega \f$ of size \f$ \left \lfloor{n/2}\right \rfloor + 1 \f$
of the transformed function
@param tmax maximum time \f$ t_\max \$f of the sampled \f$ f(t) \$f
@param nt size \f$ n \f$ of the sampled function \f$ f(t) \$f
@return 0 on success
@see inv_laplace_trans
@see http://www.fftw.org/fftw3_doc/
*/
int
laplace_trans (double* f, _Complex double* g, _Complex double* s, double tmax, size_t nt);

/** Filter to use in Inverse Numerical Laplace Transform. See Table 1 of:
GÓMEZ, Pablo; URIBE, Felipe A. The numerical Laplace transform: An accurate
technique for analyzing electromagnetic transients on power system devices.
International Journal of Electrical Power & Energy Systems, v. 31, n. 2-3,
p. 116-123, 2009.
*/
enum
INLT_Filter {
	  /** No filter */
		FILTER_NONE = 0,
    /** Blackman filter
    \f$ \sigma(\omega) = 0.42 + 0.5 \cos(\pi \omega/\Omega) + 0.08 \cos(2\pi \omega/\Omega) \f$ */
	  FILTER_BLACKMAN = 1,
	  /** Hanning filter
	  \f$ \sigma(\omega) = (1 +  \cos(\pi \omega/\Omega))/2 \f$ */
    FILTER_HANNING = 2,
    /** Lanczos filter
		\f$ \sigma(\omega) = \sin(\pi \omega/\Omega)/(\pi \omega/\Omega) \f$ */
    FILTER_LANCZOS = 3,
		/** Riez filter
		\f$ \sigma(\omega) = 1 - |\omega/\Omega|^2 \f$ */
    FILTER_RIESZ = 4
};

/** Numerical Inverse Laplace Transform \f$ f(t) = \mathcal{L}^{-1}\{g(s)\} \f$
with a real output. Uses FFTW, therefore is not thread-safe! Based on:
GÓMEZ, Pablo; URIBE, Felipe A. The numerical Laplace transform: An accurate
technique for analyzing electromagnetic transients on power system devices.
International Journal of Electrical Power & Energy Systems, v. 31, n. 2-3,
p. 116-123, 2009.

@param f inverse transfomed function \f$ f(t) \f$ of size \f$ n \f$ from
\f$ t = 0 \f$ to \f$ t = t_\max \f$ with uniform spacing \f$ \Delta t \f$
@param g sampled function \f$ g(s) \f$ of size \f$ \left \lfloor{n/2}\right \rfloor + 1 \f$
@param s ``frequencies'' \f$ s_k = c + j\omega_k \f$ of size \f$ \left \lfloor{n/2}\right \rfloor + 1 \f$,
the value of \f$ c \f$ is assumed constant
@param tmax maximum time \f$ t_\max \$f of the inverse transformed function \f$ f(t) \$f
@param nt size \f$ n \f$ of the inverse transformed function \f$ f(t) \$f
@param filter enum INLT_Filter that defines the \f$ \sigma(\omega) \f$
filter (data window) to apply to \f$ g(s) := \sigma(\omega) g(s)\f$
@return 0 on success
@see laplace_trans
@see INLT_Filter
@see http://www.fftw.org/fftw3_doc/
@see planned_inv_laplace
*/
int
inv_laplace_trans (double* f, _Complex double* g, _Complex double* s, double tmax,
	                 size_t nt, int filter);

/** Heidler function to create lightning current waveforms [1]. For parameters'
values, see e.g. [2]. Calculates
\f$
i(t) = \frac{I_0}{\xi} \frac{(t / \tau_1)^n}{1 + (t / \tau_1)^n} e^{-t / \tau_2}
\f$
where
\f$ \xi = e^{-(\tau_1 / \tau_2) (n\,\tau_2 / \tau_1)^{1 / n}} \f$

[1] HEIDLER, Fridolin; CVETIĆ, J. A class of analytical functions to study the
lightning effects associated with the current front. European transactions on
electrical power, v. 12, n. 2, p. 141-150, 2002. doi: 10.1002/etep.4450120209

[2] A. De Conti and S. Visacro, "Analytical Representation of Single- and
Double-Peaked Lightning Current Waveforms," in IEEE Transactions on
Electromagnetic Compatibility, vol. 49, no. 2, pp. 448-451, May 2007,
doi: 10.1109/TEMC.2007.897153.

@param t time $t$ in seconds
@param imax current peak \f$ I_0 \f$ in Amperes
@param tau1 rise time \f$ \tau_1 \f$ in seconds
@param tau2 decay time \f$ \tau_2 \f$ in seconds
@param n steepness expoent
@return current \f$ i(t) \f$ in Amperes
*/
double
heidler (double t, double imax, double tau1, double tau2, int n);

#endif /* AUXILIARY_H_ */
