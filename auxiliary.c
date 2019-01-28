#include <auxiliary.h>
#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <float.h>

double* linspace(double a, double b, size_t n, double u[])
{
    double c;
    size_t i;
    if (n < 2 || u == 0) /* make sure number of points and array are valid */
    {
        return (void*)0;
    }
    c = (b - a)/(n - 1); /* step size */
    for(i = 0; i < n - 1; ++i) /* fill vector */
    {
        u[i] = a + i*c;
    }
    u[n - 1] = b; /* fix last entry to b */
    return u; /* done */
}

double* logspace(double a, double b, size_t n, double u[])
{
    double c;
    size_t i;
    if (n < 2 || u == 0) /* make sure number of points and array are valid */
    {
        return (void*)0;
    }
    c = (b - a)/(n - 1); /* step size */
    for(i = 0; i < n - 1; ++i) /* fill vector */
    {
        u[i] = pow(10., a + i*c);
    }
    u[n - 1] = pow(10., b); /* fix last entry to 10^b */
    return u; /* done */
}

double wave_length(double f, double sigma, double eps, double mur)
{
    double x = sqrt(1 + pow(sigma/(TWO_PI*f*eps), 2.0));
    double y = sqrt(eps*mur*MU0*(1 + x));
    return sqrt(2.)/(f*y);
}

int equal_points(double point_1[3], double point_2[3])
{
	for (size_t i = 0; i < 3; i++)
    {
        if (fabs(point_1[i] - point_2[i]) > FLT_EPSILON)
        {
            return 0;
        }
    }
    return 1;
}

double vector_norm(double start_point[3], double end_point[3])
{
    if (equal_points(start_point, end_point))
    {
        return 0.0;
    }
    else
    {
        double length = 0.0;
        for (int i = 0; i < 3; i++)
        {
            length += pow(start_point[i] - end_point[i], 2.0);
        }
        length = sqrt(length);
        return length;
    }
}

int copy_array(_Complex double* source, _Complex double* target,
	           size_t size)
{
    for (size_t i = 0; i < size; i++)
    {
        target[i] = source[i];
    }
    return 0;
}

int print_zmatrix_file(size_t m, size_t n,
	                   _Complex double* a, int lda, FILE* fp)
{
    size_t i, j;
    for( i = 0; i < m; i++ )
    {
        for( j = 0; j < n; j++ )
        {
            fprintf(fp, "(%.6f%+.6fj)", creal(a[i*lda+j]), cimag(a[i*lda+j]) );
            if (j < n - 1) fprintf(fp, ",");
        }
        fprintf(fp, "\n");
    }
    return 0;
}

int print_dmatrix_file(size_t m, size_t n,
	                   double* a, int lda, FILE* fp)
{
    size_t i, j;
    for( i = 0; i < m; i++ )
    {
        for( j = 0; j < n; j++ )
        {
            fprintf(fp, "%.6f", a[i*lda+j]);
            if (j < n - 1) fprintf(fp, ",");
        }
        fprintf(fp, "\n");
    }
    return 0;
}
