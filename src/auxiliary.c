#include "auxiliary.h"
#include <math.h>
#include <float.h>
#include <complex.h>
#include <stdio.h>

double *
linspace (double a, double b, size_t n, double u[])
{
    if (n < 2 || u == 0) { /* make sure number of points and array are valid */
        return (void*)0;
    }
    double c = (b - a)/(n - 1); /* step size */
    for (size_t i = 0; i < n - 1; ++i) {
        u[i] = a + i*c;
    }
    u[n - 1] = b; /* fix last entry to b */
    return u;
}

double *
logspace (double a, double b, size_t n, double u[])
{
    if (n < 2 || u == 0) { /* make sure number of points and array are valid */
        return (void*)0;
    }
    double c = (b - a)/(n - 1); /* step size */
    for (size_t i = 0; i < n - 1; ++i) {
        u[i] = pow(10.0, a + i*c);
    }
    u[n - 1] = pow(10.0, b); /* fix last entry to 10^b */
    return u;
}

double
wave_length (double f, double sigma, double eps, double mur)
{
    double x = sqrt(1.0 + pow(sigma/(TWO_PI*f*eps), 2.0));
    double y = sqrt(eps*mur*MU0*(1.0 + x));
    return sqrt(2.0)/(f*y);
}

int
equal_points (const double point_1[3], const double point_2[3])
{
    for (int i = 0; i < 3; i++) {
        if (fabs(point_1[i] - point_2[i]) > FLT_EPSILON) return 0;
    }
    return 1;
}

double
vector_norm (const double start_point[3], const double end_point[3])
{
    double length = 0.0;
	for (int i = 0; i < 3; i++) {
        length += pow(start_point[i] - end_point[i], 2.0);
	}
	length = sqrt(length);
	return length;
}

int
complex_matrix_file(size_t m, size_t n, const _Complex double *a, int lda, FILE *fp)
{
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            fprintf(fp, "(%.6f%+.6fj)", creal(a[i*lda+j]), cimag(a[i*lda+j]) );
            if (j < n - 1) fprintf(fp, ",");
        }
        fprintf(fp, "\n");
    }
    return 0;
}

int
double_matrix_file(size_t m, size_t n, const double *a, int lda, FILE *fp)
{
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            fprintf(fp, "%.6f", a[i*lda+j]);
            if (j < n - 1) fprintf(fp, ",");
        }
        fprintf(fp, "\n");
    }
    return 0;
}
