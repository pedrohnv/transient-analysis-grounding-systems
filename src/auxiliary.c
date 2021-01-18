/* High Performance implementation of the Hybrid Electromagnetic Model
Released under the General Public License 3 (GPLv3).
All parameters' units are in the SI base if omitted.

Constants, auxiliary functions and routines.
*/
#include "auxiliary.h"
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include <stdio.h>
#include <stdbool.h>
#include <fftw3.h>

double*
linspace (double a, double b, size_t n, double u[])
{
    if (n < 2 || u == 0) {  // make s number of points and array are valid
        return (void*)0;
    }
    double c = (b - a) / (n - 1);  // step size
    for (size_t i = 0; i < n - 1; ++i) {
        u[i] = a + i * c;
    }
    u[n - 1] = b;  // fix last entry to b
    return u;
}

double*
logspace (double a, double b, size_t n, double u[])
{
    if (n < 2 || u == 0) {  // make sure number of points and array are valid
        return (void*)0;
    }
    double c = (b - a) / (n - 1);  // step size
    for (size_t i = 0; i < n - 1; ++i) {
        u[i] = pow(10.0, a + i * c);
    }
    u[n - 1] = pow(10.0, b);  // fix last entry to 10^b
    return u;
}

double
wave_length (double f, double sigma, double ep, double mur)
{
    double x = sqrt(1.0 + pow(sigma / (TWO_PI * f * ep), 2.0));
    double y = sqrt(ep * mur * MU0 * (1.0 + x));
    return sqrt(2.0) / (f * y);
}

int
alipio_soil (_Complex double* sigma, _Complex double* epsr, double sigma0,
             _Complex double s, double h, double g, double eps_ratio)
{

    _Complex double f = s / TWO_PI;
    *sigma = (sigma0 + sigma0 * h * cpow(f/1e6, g));
    double t = tan(PI * g / 2) / (TWO_PI * EPS0 * pow(1e6, g));
    *epsr = eps_ratio + t * sigma0 * h * cpow(f, g - 1.0);
    return 0;
}

int
smith_longmire_soil (_Complex double* sigma, _Complex double* epsr, double sigma0,
                     _Complex double s, double erinf)
{
    int N = 13;
    double a[] = {3.4e6, 2.74e5, 2.58e4, 3.38e3, 5.26e2, 1.33e2, 2.72e1, 1.25e1,
                  4.8e0, 2.17e0, 9.8e-1, 3.92e-1, 1.73e-1};
    _Complex double Fdc = cpow(125.0 * sigma0, 0.8312);
    _Complex double sum_epsr = 0.0;
    _Complex double sum_sigma = 0.0;
    _Complex double F, fratio2, den;
    for (int i = 0; i < N; i++) {
        F = Fdc * cpow(10.0, i);
        fratio2 = cpow(s / (I * TWO_PI * F), 2.0);
        den = (1.0 + fratio2);
        sum_epsr += a[i] / den;
        sum_sigma += a[i] * F * (fratio2 / den);
    }
    *epsr = erinf + sum_epsr;
    *sigma = sigma0 + TWO_PI * EPS0 * sum_sigma;
    return 0;
}

bool
equal_points (const double* point_1, const double* point_2)
{
    return equal_points_tol (point_1, point_2, DBL_EPSILON);
}

bool
equal_points_tol (const double* point_1, const double* point_2, double tol)
{
    for (int i = 0; i < 3; i++) {
        if (fabs(point_1[i] - point_2[i]) > tol) return false;
    }
    return true;
}

double
vector_length (const double start_point[3], const double end_point[3])
{
    /*int n = 3;
    double x[n];
    for (int i = 0; i < 3; i++) {
        x[i] = start_point[i] - end_point[i];
    }
    int incx = 1;
    return snrm2_(&n, x, &incx);*/
    double r = 0.0;
    for (int i = 0; i < 3; i++) {
        r += pow(start_point[i] - end_point[i], 2.0);
    }
    return sqrt(r);
}

int
complex_matrix_file (size_t m, size_t n, const _Complex double* a, int lda, FILE* fp)
{
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            fprintf(fp, "(%.6g %+.6gj)", creal(a[i * lda + j]), cimag(a[i * lda + j]) );
            if (j < n - 1) fprintf(fp, ",");
        }
        fprintf(fp, "\n");
    }
    return 0;
}

int
double_matrix_file (size_t m, size_t n, const double* a, int lda, FILE* fp)
{
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            fprintf(fp, "%.6g", a[i * lda + j]);
            if (j < n - 1) fprintf(fp, ",");
        }
        fprintf(fp, "\n");
    }
    return 0;
}

/* Auxiliary routine: printing a matrix COLUMN MAJOR*/
void
print_matrix (char *desc, int m, int n, const _Complex double* a, int lda)
{
    printf( "\n transpose of %s\n", desc );
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++)
            printf( "%10.5g %+10.5gim, ", creal(a[j + i * lda]), cimag(a[j + i * lda]) );
            //printf( "%6.2f, ", creal(a[j + i * lda]) );
        printf( "\n" );
    }
}

/* Auxiliary routine: printing a matrix ROW MAJOR*/
void
print_matrix_row (char *desc, int m, int n, const _Complex double* a, int lda)
{
    printf( "\n %s\n", desc );
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++)
            printf( "%10.5g %+10.5gim, ", creal(a[i * lda + j]), cimag(a[i * lda + j]) );
            //printf( "%6.2f, ", creal(a[i * lda + j]) );
        printf( "\n" );
    }
}

int
matrix_copy (const _Complex double* source, _Complex double* target,
             size_t lds, size_t ldt, size_t nline, size_t ncol)
{
    for (size_t i = 0; i < ncol; i++) {
        for (size_t k = 0; k < nline; k++) {
            target[i*ldt + k] = source[i * lds + k];
        }
    }
    return 0;
}

int
transpose_copy (const _Complex double* source, _Complex double* target,
                size_t lds, size_t ldt, size_t nline, size_t ncol)
{
    for (size_t i = 0; i < nline; i++) {
        for (size_t k = 0; k < ncol; k++) {
            target[i*ldt + k] = source[k * lds + i];
        }
    }
    return 0;
}

int
pc_copy (const _Complex double* source, _Complex double* target,
         size_t lds, size_t ldt, size_t nline, size_t ncol)
{
    for (size_t i = 0; i < ncol; i++) {
        for (size_t k = 0; k < nline; k++) {
            target[i*ldt + k] = source[(ncol - i - 1) * lds + k];
        }
    }
    return 0;
}

int
pl_copy (const _Complex double* source, _Complex double* target,
         size_t lds, size_t ldt, size_t nline, size_t ncol)
{
    for (size_t i = 0; i < ncol; i++) {
        for (size_t k = 0; k < nline; k++) {
            target[i*ldt + k] = source[i * lds + (nline - k - 1)];
        }
    }
    return 0;
}

int
pcl_copy (const _Complex double* source, _Complex double* target,
          size_t lds, size_t ldt, size_t nline, size_t ncol)
{
    for (size_t i = 0; i < ncol; i++) {
        for (size_t k = 0; k < nline; k++) {
            target[i*ldt + k] = source[(ncol - i - 1) * lds + (nline - k - 1)];
        }
    }
    return 0;
}

int
laplace_trans (double* f, _Complex double* g, _Complex double* s, double tmax, size_t nt)
{
    int ns = nt / 2 + 1;
    double* in = fftw_malloc(nt * sizeof(double));
    fftw_complex* out = fftw_malloc(ns * sizeof(fftw_complex));
    // You must create the plan before initializing the input
    fftw_plan p = fftw_plan_dft_r2c_1d(nt, in, out, FFTW_ESTIMATE);
    double c = log(pow(nt, 2.0)) / tmax;
    double dt = tmax / (nt - 1);
    double dw = TWO_PI / tmax;
    for (size_t k = 0; k < ns; k++) {
        s[k] = c + I * dw * k;
    }
    for (size_t k  = 0; k < nt; k++) {
        in[k] = f[k] * dt * exp(-c * k * dt);
    }
    fftw_execute(p);
    // fftw_complex is _Complex double if fftw3.h is included after complex.h
    for (int k = 0; k < ns; k++) g[k] = out[k];
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    return 0;
}

int
inv_laplace_trans (double* f, _Complex double* g, _Complex double* s, double tmax,
                   size_t nt, int filter)
{
    int ns = nt / 2 + 1;
    //double* out = malloc(nt * sizeof(double));
    //_Complex double* in = malloc(ns * sizeof(_Complex double));
    double* out = fftw_malloc(nt * sizeof(double));
    fftw_complex* in = fftw_malloc(ns * sizeof(fftw_complex));
    // You must create the plan before initializing the input
    fftw_plan plan = fftw_plan_dft_c2r_1d(nt, in, out, FFTW_ESTIMATE);
    double dt = tmax / (nt - 1);
    double omega = cimag(s[ns - 1]);
    double w, sigma, alpha;
    switch (filter) {
        case FILTER_BLACKMAN:
            for (size_t k = 0; k < ns; k++) {
                w = cimag(s[k]);
                alpha = PI * w / omega;
                sigma = 0.42 + 0.5 * cos(alpha) + 0.08 * cos(2 * alpha);
                in[k] = sigma * g[k];
            }
            break;

        case FILTER_HANNING:
            for (size_t k = 0; k < ns; k++) {
                w = cimag(s[k]);
                sigma = (1 + cos(PI * w / omega)) / 2;
                in[k] = sigma * g[k];
            }
            break;

        case FILTER_LANCZOS:
            for (size_t k = 0; k < ns; k++) {
                w = cimag(s[k]);
                alpha = PI * w / omega;
                sigma = sin(alpha) / alpha;
                in[k] = sigma * g[k];
            }
            break;

        case FILTER_RIESZ:
            for (size_t k = 0; k < ns; k++) {
                w = cimag(s[k]);
                sigma = 1.0 - pow(fabs(w / omega), 2.0);
                in[k] = sigma * g[k];
            }
            break;

        default:
            for (size_t k = 0; k < ns; k++) in[k] = g[k];
            break;
    }
    fftw_execute(plan);
    double c = creal(s[0]);
    for (size_t k = 0; k < nt; k++) {
        f[k] = out[k] * exp(c * k * dt) / (nt * dt);
    }
    fftw_destroy_plan(plan);
    //free(in);
    //free(out);
    fftw_free(in);
    fftw_free(out);
    return 0;
}

double heidler (double t, double imax, double tau1, double tau2, int n)
{
    double xi = exp( -(tau1 / tau2) * pow((n * tau2 / tau1), (1.0 / n)) );
    double tt1n = pow(t / tau1, n);
    return imax / xi * tt1n / (1 + tt1n) * exp(-t / tau2);
}
