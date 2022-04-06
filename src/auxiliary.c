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

float*
linspace (float a, float b, size_t n, float u[])
{
    if (n < 2 || u == 0) {  // make s number of points and array are valid
        return (void*)0;
    }
    float c = (b - a) / (n - 1);  // step size
    for (size_t i = 0; i < n - 1; ++i) {
        u[i] = a + i * c;
    }
    u[n - 1] = b;  // fix last entry to b
    return u;
}

float*
logspace (float a, float b, size_t n, float u[])
{
    if (n < 2 || u == 0) {  // make sure number of points and array are valid
        return (void*)0;
    }
    float c = (b - a) / (n - 1);  // step size
    for (size_t i = 0; i < n - 1; ++i) {
        u[i] = powf(10.0, a + i * c);
    }
    u[n - 1] = powf(10.0, b);  // fix last entry to 10^b
    return u;
}

float
wave_length (float f, float sigma, float ep, float mur)
{
    float x = sqrtf(1.0 + powf(sigma / (TWO_PI * f * ep), 2.0));
    float y = sqrtf(ep * mur * MU0 * (1.0 + x));
    return sqrtf(2.0) / (f * y);
}

int
alipio_soil (_Complex float* sigma, _Complex float* epsr, float sigma0,
             _Complex float s, float h, float g, float eps_ratio)
{

    _Complex float f = s / (TWO_PI * I);
    *sigma = (sigma0 + sigma0 * h * cpowf(f/1e6, g));
    float t = tanf(PI * g / 2) / (TWO_PI * EPS0 * powf(1e6, g));
    *epsr = eps_ratio + t * sigma0 * h * cpowf(f, g - 1.0);
    return 0;
}

int
smith_longmire_soil (_Complex float* sigma, _Complex float* epsr, float sigma0,
                     _Complex float s, float erinf)
{
    int N = 13;
    float a[] = {3.4e6, 2.74e5, 2.58e4, 3.38e3, 5.26e2, 1.33e2, 2.72e1, 1.25e1,
                  4.8e0, 2.17e0, 9.8e-1, 3.92e-1, 1.73e-1};
    _Complex float Fdc = cpowf(125.0 * sigma0, 0.8312);
    _Complex float sum_epsr = 0.0;
    _Complex float sum_sigma = 0.0;
    _Complex float F, fratio2, den;
    for (int i = 0; i < N; i++) {
        F = Fdc * cpowf(10.0, i);
        fratio2 = cpowf(s / (I * TWO_PI * F), 2.0);
        den = (1.0 + fratio2);
        sum_epsr += a[i] / den;
        sum_sigma += a[i] * F * (fratio2 / den);
    }
    *epsr = erinf + sum_epsr;
    *sigma = sigma0 + TWO_PI * EPS0 * sum_sigma;
    return 0;
}

bool
equal_points (const float* point_1, const float* point_2)
{
    return equal_points_tol (point_1, point_2, FLT_EPSILON);
}

bool
equal_points_tol (const float* point_1, const float* point_2, float tol)
{
    for (int i = 0; i < 3; i++) {
        if (fabsf(point_1[i] - point_2[i]) > tol) return false;
    }
    return true;
}

float
vector_length (const float start_point[3], const float end_point[3])
{
    /*int n = 3;
    float x[n];
    for (int i = 0; i < 3; i++) {
        x[i] = start_point[i] - end_point[i];
    }
    int incx = 1;
    return snrm2_(&n, x, &incx);*/
    float r = 0.0;
    for (int i = 0; i < 3; i++) {
        r += powf(start_point[i] - end_point[i], 2.0);
    }
    return sqrtf(r);
}

int
complex_matrix_file (size_t m, size_t n, const _Complex float* a, int lda, FILE* fp)
{
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            fprintf(fp, "(%.6g %+.6gj)", crealf(a[i * lda + j]), cimagf(a[i * lda + j]) );
            if (j < n - 1) fprintf(fp, ",");
        }
        fprintf(fp, "\n");
    }
    return 0;
}

int
float_matrix_file (size_t m, size_t n, const float* a, int lda, FILE* fp)
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
print_matrix (char *desc, int m, int n, const _Complex float* a, int lda)
{
    printf( "\n transpose of %s\n", desc );
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++)
            printf( "%10.5g %+10.5gim, ", crealf(a[j + i * lda]), cimagf(a[j + i * lda]) );
            //printf( "%6.2f, ", crealf(a[j + i * lda]) );
        printf( "\n" );
    }
}

/* Auxiliary routine: printing a matrix ROW MAJOR*/
void
print_matrix_row (char *desc, int m, int n, const _Complex float* a, int lda)
{
    printf( "\n %s\n", desc );
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++)
            printf( "%10.5g %+10.5gim, ", crealf(a[i * lda + j]), cimagf(a[i * lda + j]) );
            //printf( "%6.2f, ", crealf(a[i * lda + j]) );
        printf( "\n" );
    }
}

int
matrix_copy (const _Complex float* source, _Complex float* target,
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
transpose_copy (const _Complex float* source, _Complex float* target,
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
pc_copy (const _Complex float* source, _Complex float* target,
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
pl_copy (const _Complex float* source, _Complex float* target,
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
pcl_copy (const _Complex float* source, _Complex float* target,
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

double heidler (double t, double imax, double tau1, double tau2, double n)
{
    double xi = exp( -(tau1 / tau2) * pow((n * tau2 / tau1), (1.0 / n)) );
    double tt1n = pow(t / tau1, n);
    return imax / xi * tt1n / (1 + tt1n) * exp(-t / tau2);
}

int lines_in_file (const char filename[])
{
    FILE *fp;
    int count = 0;  // Line counter (result)
    char c;  // To store a character read from file
    // Open the file
    fp = fopen(filename, "r");
    // Check if file exists
    if (fp == NULL)
    {
        printf("Could not open file %s", filename);
        return 0;
    }
    // Extract characters from file and store in character c
    for (c = getc(fp); c != EOF; c = getc(fp))
        if (c == '\n')  // Increment count if this character is newline
            count = count + 1;
    // Close the file
    fclose(fp);
    //printf("The file %s has %d lines\n ", filename, count);
    return count;
}
