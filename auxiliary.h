/**
High performance Hybrid Electromagnetic Model calculations in C.

All parameters' units are in SI if omitted.

@author Pedro Henrique Nascimento Vieira

Constants, auxiliary functions and routines.
*/

#ifndef AUXILIARY_H_
#define AUXILIARY_H_

#include <complex.h>
#include <stdio.h>

//constants
#define PI 3.1415926535897932384626433832795029L
#define TWO_PI 6.283185307179586
#define FOUR_PI 12.56637061435917
#define MU0 1.256637061435917e-6 //permeability vac.
#define EPS0 8.854187817620e-12 //permittivity vac.

/** linspace
Fill array `u` with `n` linearly spaced numbers between `a` and `b` (included).
@see Source code taken from https://github.com/ntessore/algo
@param a first number
@param b last number
@param n size of array
@param u array to be filled
@return pointer to filled array
@see https://github.com/ntessore/algo
*/
double* linspace(double a, double b, int n, double u[]);

/** logspace
Fill array `u` with `n` logarithmically spaced numbers between `10^a` and
`10^b` (included).
@see Source code taken from https://github.com/ntessore/algo
@param a first power
@param b last power
@param n size of array
@param u array to be filled
@return pointer to filled array
@see https://github.com/ntessore/algo
*/
double* logspace(double a, double b, int n, double u[]);

/** wave_length
Computes the harmonic electromagnetic wave length in a lossy medium.
@see Source code taken from https://github.com/ntessore/algo
@param f frequency (Hz)
@param sigma medium (real) conductivity (Ohm*m)
@param eps medium electric permittivity
@param mu medium magnetic permeability
@return lambda wave length (m)
*/
double wave_length(double f, double sigma, double eps, double mu);

/** equal_points
Check if two points are the same.
@param point_1 first point in \f$ R^3 \f$
@param point_2 second point in \f$ R^3 \f$
@return identity 0 if false, 1 if true
*/
int equal_points(double point_1[3], double point_2[3]);

/** vector_norm
Computes the Norm_2 (euclidian length) of the vector which starts at
start_point ends at end_point.
@param start_point array `(x,y,z)` that defines the starting point
@param end_point array `(x,y,z)` that defines the ending point
@return length Norm_2
*/
double vector_norm(double start_point[3], double end_point[3]);

/** copy_array
Copy `source` array into `target` array of complex numbers.
@param source pointer to the array to be copied
@param target pointer to array which the copy is to be stored
@para size size of the array
@return 0 on success
*/
int copy_array(_Complex double* source, _Complex double* target, int size);

/** print_zmatrix_file
Prints a complex (double) matrix to a file.
@param m number of rows
@param n number of columns
@param a pointer to array where matrix is stored
@param lda leading dimension of `a`
@param fp pointer to file stream
@return 0 on success
*/
int print_zmatrix_file(int m, int n, _Complex double* a, int lda, FILE* fp);

/** print_dmatrix_file
Prints a real (double) matrix to a file.
@param m number of rows
@param n number of columns
@param a pointer to array where matrix is stored
@param lda leading dimension of `a`
@param fp pointer to file stream
@return 0 on success
*/
int print_dmatrix_file(int m, int n, double* a, int lda, FILE* fp);

/** zbesi_
FORTRAN subroutine to calculate I-Bessel function, i.e., modified Bessel
function of the first kind, with complex argument.
    \f$ cy = I_{fnu} (z) \f$
@param zr argument's real part Re(z)
@param zi argument's imaginary part Im(z)
@param fnu order of initial I function
@param kode scaling, 1: no scaling, 2: \f$e^(-|x|)\f$
@param n number of members of the sequence
@param cyr result's real part Re(cy)
@param cyi result's imaginary part Im(cy)
@param nz number of components set to zero due to underflow
@param ierr error flag
@see http://netlib.org/amos/zbesi.f
*/
extern int zbesi_(double* zr, double* zi, double* fnu, int* kode, int* n,
    double* cyr, double* cyi, int* nz, int* ierr);

#endif /* AUXILIARY_H_ */
