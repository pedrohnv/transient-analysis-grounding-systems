#include <Electrode.h>
#include <auxiliary.h>
#include <cubature.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <mkl.h>
#include <mkl_lapacke.h>
//#include <omp.h>

int populate_electrode(Electrode* electrode, const double start_point[3],
                       const double end_point[3], double radius,
                       _Complex double internal_impedance)
{
    if (equal_points(start_point, end_point))
    {
        printf("Error: start_point and end_point are equal.\n");
        return 1;
    }
    if (radius <= 0)
    {
        printf("Error: radius < 0.\n");
        return 2;
    }
    electrode->radius = radius;
    electrode->zi = internal_impedance;
    electrode->length = vector_norm(start_point, end_point);
    for (size_t i = 0; i < 3; i++)
    {
        electrode->start_point[i] = start_point[i];
        electrode->end_point[i] = end_point[i];
        electrode->middle_point[i] = (start_point[i] + end_point[i])/2.0;
    }
    return 0;
}

int electrodes_file(const char file_name[], Electrode* electrodes,
                    size_t num_electrodes)
{
    FILE* stream = fopen(file_name, "r");
    if (stream == NULL) {
        printf("Cannot open file %s\n", file_name);
        exit(1);
    }
    double start_point[3], end_point[3];
    double radius, rezi, imzi;
    int success = 9;
    int pop_error = 0;
    for (size_t i = 0; i < num_electrodes; i++)
    {
        success = fscanf(stream, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
                         &start_point[0], &start_point[1], &start_point[2],
                         &end_point[0], &end_point[1], &end_point[2],
                         &radius, &rezi, &imzi);
        if (success != 9)
        {
            printf("error reading line %i of file %s\n", (int) (i+1), file_name);
            break;
        }
        success = populate_electrode(&(electrodes[i]), start_point, end_point,
                                     radius, rezi + I*imzi);
        if (success != 0)
        {
            printf("Bad input: could not create electrode %i from file %s\n",
                   (int) (i+1), file_name);
            pop_error = success;
            break;
        }
    }
    fclose(stream);
    if (pop_error == 0)
    {
        return (success - 9);
    } else {
        return (pop_error);
    }
}

int nodes_file(const char file_name[], double nodes[][3], size_t num_nodes)
{
    FILE* stream = fopen(file_name, "r");
    if (stream == NULL) {
        printf("Cannot open file %s\n", file_name);
        exit(1);
    }
    int success = 3;
    for (size_t i = 0; i < num_nodes; i++)
    {
        success = fscanf(stream, "%lf %lf %lf", &nodes[i][0], &nodes[i][1],
                         &nodes[i][2]);
        if (success != 3)
        {
            printf("error reading line %i of file %s\n", (int) (i+1), file_name);
            break;
        }
    }
    fclose(stream);
    return (success - 3);
}

int segment_electrode(Electrode* electrodes, double nodes[][3],
                      size_t num_segments, const double* start_point,
                      const double* end_point, double radius,
                      _Complex double unit_zi)
{
    if (num_segments < 1)
    {
        printf("Error: number of segments should be greater than 0.\n");
        return 1;
    }
    size_t num_nodes = num_segments + 1;
    double x[num_nodes];
    double y[num_nodes];
    double z[num_nodes];
    linspace(start_point[0], end_point[0], num_nodes, x);
    linspace(start_point[1], end_point[1], num_nodes, y);
    linspace(start_point[2], end_point[2], num_nodes, z);
    size_t i;
    for (i = 0; i < num_nodes; i++)
    {
        nodes[i][0] = x[i];
        nodes[i][1] = y[i];
        nodes[i][2] = z[i];
    }
    double total_length = vector_norm(start_point, end_point);
    _Complex double zi = unit_zi*total_length/num_segments;
    for (i = 0; i < num_segments; i++)
    {
        populate_electrode(&(electrodes[i]), nodes[i], nodes[i + 1], radius, zi);
    }
    return 0;
}

int integrand_double(unsigned ndim, const double *t, void *auxdata,
                     unsigned fdim, double *fval)
{
    Integration_data *p = (Integration_data*) auxdata;
    const Electrode* sender = p->sender;
    const Electrode* receiver = p->receiver;
    _Complex double gamma = p->gamma;
    double point_r, point_s, r = 0.0;
    for (size_t i = 0; i < 3; i++)
    {
        point_r = t[0]*(receiver->end_point[i] - receiver->start_point[i]);
        point_r += receiver->start_point[i];

        point_s = t[1]*(sender->end_point[i] - sender->start_point[i]);
        point_s += sender->start_point[i];

        r += pow(point_r - point_s, 2.0);
    }
    r = sqrt(r);
    _Complex double exp_gr = cexp(-gamma*r);
    fval[0] = creal(exp_gr)/r;
    fval[1] = cimag(exp_gr)/r;
    return 0;
}

int exp_logNf(unsigned ndim, const double *t, void *auxdata, unsigned fdim,
              double *fval)
{
    Integration_data *p = (Integration_data*) auxdata;
    const Electrode* sender = p->sender;
    const Electrode* receiver = p->receiver;
    double r1, r2, eta, point_r;
    r1 = 0.0;
    r2 = 0.0;
    eta = 0.0;
    for (size_t i = 0; i < 3; i++)
    {
        point_r = t[0]*(receiver->end_point[i] - receiver->start_point[i]);
        point_r += receiver->start_point[i];
        r1 += pow(point_r - sender->start_point[i], 2.0);
        r2 += pow(point_r - sender->end_point[i], 2.0);
        eta += pow(point_r - sender->middle_point[i], 2.0);
    }
    r1 = sqrt(r1);
    r2 = sqrt(r2);
    eta = sqrt(eta);
    _Complex double exp_gr;
    if (p->simplified)
    {
        exp_gr = 1.0;
    }
    else
    {
        exp_gr = cexp(- p->gamma*eta);
    }
    double Nf = (r1 + r2 + sender->length)/(r1 + r2 - sender->length);
    exp_gr = exp_gr*log(Nf);
    fval[0] = creal(exp_gr);
    fval[1] = cimag(exp_gr);
    return 0;
}

int integral(const Electrode* sender, const Electrode* receiver,
             _Complex double gamma, size_t max_eval, double req_abs_error,
             double req_rel_error, int error_norm, int integration_type,
             double result[2], double error[2])
{
    Integration_data* auxdata = malloc(sizeof(Integration_data));
    auxdata->sender = sender;
    auxdata->receiver = receiver;
    auxdata->gamma = gamma;
    double tmin[] = {0.0, 0.0};
    double tmax[] = {1.0, 1.0};
    int failure = 1;
    switch (integration_type) {
        case NONE:
            double lslr = sender->length * receiver->length;
            double rbar = vector_norm(sender->middle_point, receiver->middle_point);
            _Complex double exp_gr = cexp(-gamma*rbar)*lslr;
            result[0] = creal(exp_gr);
            result[1] = cimag(exp_gr);
            failure = 0;
        break;

        case INTG_DOUBLE:
            failure = hcubature(2, integrand_double, auxdata, 2, tmin, tmax,
                                max_eval, req_abs_error, req_rel_error,
                                error_norm, result, error);
            result[0] = result[0] * receiver->length * sender->length;
            result[1] = result[1] * receiver->length * sender->length;
        break;

        case INTG_EXP_LOGNF:
            auxdata->simplified = 0;
            failure = hcubature(2, exp_logNf, auxdata, 1, tmin, tmax, max_eval,
                                req_abs_error, req_rel_error, error_norm,
                                result, error);
            result[0] = result[0] * receiver->length;
            result[1] = result[1] * receiver->length;
        break;

        case INTG_LOGNF:
            auxdata->simplified = 1;
            failure = hcubature(1, exp_logNf, auxdata, 1, tmin, tmax, max_eval,
                                req_abs_error, req_rel_error, error_norm,
                                result, error);
            double rbar;
            rbar = vector_norm(sender->middle_point, receiver->middle_point);
            _Complex double exp_gr = cexp(-gamma*rbar)*result[0];
            result[0] = creal(exp_gr) * receiver->length;
            result[1] = cimag(exp_gr) * receiver->length;
        break;
    }
    free(auxdata);
    return failure;
}

_Complex double internal_impedance(
    _Complex double s, double rho, double radius, double mur)
{
    _Complex double etapr = csqrt(s*mur*MU0/rho);
    _Complex double etapr_radius = etapr*radius;
    double zr = creal(etapr_radius);
    double zi = cimag(etapr_radius);
    double fnu0 = 0;
    double fnu1 = 1;
    int kode = 1;
    int n = 1;
    double cyr0, cyi0, cyr1, cyi1;
    int nz, ierr;
    zbesi_(&zr, &zi, &fnu0, &kode, &n, &cyr0, &cyi0, &nz, &ierr);
    if (ierr != 0)
    {
        printf("in function 'internal_impedance'\n");
        printf("error in call to zbesi(fnu = 0)\n");
        printf("error flag = %i\n", ierr);
        exit(ierr);
    }
    zbesi_(&zr, &zi, &fnu1, &kode, &n, &cyr1, &cyi1, &nz, &ierr);
    if (ierr != 0)
    {
        printf("in function 'internal_impedance'\n");
        printf("error in call to zbesi(fnu = 1)\n");
        printf("error flag = %i\n", ierr);
        exit(ierr);
    }
    return (etapr*rho*(cyr0 + cyi0*I))/(TWO_PI*radius*(cyr1 + cyi1*I));
}

// Longitudinal impedance
_Complex double longitudinal_self(const Electrode* electrode,
                                  _Complex double s, double mur)
{
    double ls = electrode->length;
    double k1 = electrode->radius/ls;
    double k2 = sqrt(1.0 + k1*k1);
    _Complex double zi = electrode->zi;
    return s*mur*MU0*ls/(TWO_PI)*(log( (k2 + 1.)/k1 ) - k2 + k1) + zi;
}

_Complex double longitudinal_mutual(const Electrode* sender,
                                    const Electrode* receiver,
                                    _Complex double s, double mur,
                                    _Complex double gamma, size_t max_eval,
                                    double req_abs_error, double req_rel_error,
                                    int error_norm, int integration_type,
                                    double result[2], double error[2])
{
    double cost = 0.0; //cosine between segments
    double k1, k2;
    for (size_t i = 0; i < 3; i++)
    {
        k1 = (sender->end_point[i] - sender->start_point[i]);
        k2 = (receiver->end_point[i] - receiver->start_point[i]);
        cost += k1*k2;
    }
    cost = abs(cost/(sender->length * receiver->length));
    if (fabs(cost - 0.0) < DBL_EPSILON)
    {
        return 0.0;
    }
    else
    {
        integral(sender, receiver, gamma, max_eval, req_abs_error,
                 req_rel_error, error_norm, integration_type, result, error);
        _Complex double intg = result[0] + I*result[1];
        return s*mur*MU0/(FOUR_PI)*intg*cost;
    }
}

// Transveral impedance
_Complex double transversal_self(const Electrode* electrode,
                                 _Complex double kappa)
{
    double ls = electrode->length;
    double k1 = electrode->radius/ls;
    double k2 = sqrt(1.0 + k1*k1);
    return 1.0/(TWO_PI*kappa*ls)*(log( (k2 + 1.)/k1 ) - k2 + k1);
}

_Complex double transversal_mutual(const Electrode* sender,
                                   const Electrode* receiver,
                                   _Complex double kappa,
                                   _Complex double gamma, size_t max_eval,
                                   double req_abs_error, double req_rel_error,
                                   int error_norm, int integration_type,
                                   double result[2], double error[2])
{
    double ls = sender->length;
    double lr = receiver->length;
    integral(sender, receiver, gamma, max_eval, req_abs_error,
             req_rel_error, error_norm, integration_type, result, error);
    _Complex double intg = result[0] + I*result[1];
    return 1.0/(FOUR_PI*kappa*ls*lr)*intg;
}

// TODO store ZT and ZL as upper/lower triangular matrices, as they are symmetric
// TODO use COLUMN_MAJOR ordering to avoid having transpose the matrices for LAPACK
int calculate_impedances(const Electrode* electrodes, size_t num_electrodes,
                         _Complex double* zl, _Complex double* zt,
                         _Complex double gamma, _Complex double s, double mur,
                         _Complex double kappa, size_t max_eval,
                         double req_abs_error, double req_rel_error,
                         int error_norm, int integration_type)
{
    double result[2], error[2], ls, lr, k1, k2, cost;
    _Complex double iwu_4pi = s*mur*MU0/(FOUR_PI);
    _Complex double one_4pik = 1.0/(FOUR_PI*kappa);
    _Complex double intg;
    size_t i, k, m;
    int failure;
    // _self and _mutual impedances are not used to reduce the number of
    // operations, as some of them would be done twice or more unnecessarily
    for (i = 0; i < num_electrodes; i++)
    {
        ls = electrodes[i].length;
        k1 = electrodes[i].radius/ls;
        k2 = sqrt(1.0 + k1*k1);
        cost = 2.0*(log( (k2 + 1.)/k1 ) - k2 + k1);
        zl[i*num_electrodes + i] = iwu_4pi*ls*cost + electrodes[i].zi;
        zt[i*num_electrodes + i] = one_4pik/ls*cost;
        for (k = i+1; k < num_electrodes; k++)
        {
            lr = electrodes[k].length;
            cost = 0.0;
            for (m = 0; m < 3; m++)
            {
                k1 = (electrodes[i].end_point[m] - electrodes[i].start_point[m]);
                k2 = (electrodes[k].end_point[m] - electrodes[k].start_point[m]);
                cost += k1*k2;
            }
            cost = abs(cost/(ls*lr));
            failure = integral(&(electrodes[i]), &(electrodes[k]), gamma,
                               max_eval, req_abs_error, req_rel_error,
                               error_norm, integration_type, result, error);
            if (failure) return failure;
            intg = result[0] + I*result[1];
            zl[i*num_electrodes + k] = iwu_4pi*intg*cost;
            zt[i*num_electrodes + k] = one_4pik/(ls*lr)*intg;

            zl[k*num_electrodes + i] = zl[i*num_electrodes + k];
            zt[k*num_electrodes + i] = zt[i*num_electrodes + k];
        }
    }
    return 0;
}

int impedances_images(const Electrode* electrodes, const Electrode* images,
                      size_t num_electrodes, _Complex double* zl,
                      _Complex double* zt, _Complex double gamma,
                      _Complex double s, double mur, _Complex double kappa,
                      _Complex double ref_l, _Complex double ref_t,
                      size_t max_eval, double req_abs_error,
                      double req_rel_error, int error_norm, int integration_type)
{
    double result[2], error[2], ls, lr, k1, k2, cost;
    _Complex double iwu_4pi = s*mur*MU0/(FOUR_PI);
    _Complex double one_4pik = 1.0/(FOUR_PI*kappa);
    _Complex double intg;
    int failure;
    // _mutual impedances are not used to reduce the number of
    // operations, as some of them would be done twice or more unnecessarily
    for (size_t i = 0; i < num_electrodes; i++)
    {
        ls = electrodes[i].length;
        for (size_t k = i; k < num_electrodes; k++)
        {
            lr = images[k].length;
            cost = 0.0;
            for (size_t m = 0; m < 3; m++)
            {
                k1 = (electrodes[i].end_point[m] - electrodes[i].start_point[m]);
                k2 = (images[k].end_point[m] - images[k].start_point[m]);
                cost += k1*k2;
            }
            cost = abs(cost/(ls*lr));
            failure = integral(&(electrodes[i]), &(images[k]), gamma, max_eval,
                               req_abs_error, req_rel_error, error_norm,
                               integration_type, result, error);
            if (failure) return failure;
            intg = result[0] + I*result[1];
            /*if (i != k)*/ zl[i*num_electrodes + k] += ref_l*iwu_4pi*intg*cost;
            zt[i*num_electrodes + k] += ref_t*one_4pik/(ls*lr)*intg;

            zl[k*num_electrodes + i] = zl[i*num_electrodes + k];
            zt[k*num_electrodes + i] = zt[i*num_electrodes + k];
        }
    }
    return 0;
}

int fill_incidence(_Complex double* we, const Electrode* electrodes,
                   size_t num_electrodes, const double nodes[][3],
                   size_t num_nodes)
{
    int e, n, condition, node_is_start, node_is_end, a, b, c, d, no_incidence;
    size_t ld = (2*num_electrodes + num_nodes)*num_electrodes;
    size_t ld_2 = 2*ld;
    //TODO use a more efficient search or fill-up method if possible
    for (n = 0; n < num_nodes; n++)
    {
        no_incidence = 1;
        for (e = 0; e < num_electrodes; e++)
        {
            node_is_start = equal_points(electrodes[e].start_point, nodes[n]);
            node_is_end = 0;
            /* it is assumed that start_point and end_point of an electrode
            are never equal */
            condition = 0;
            if (node_is_start)
            {
                condition = 1;
            }
            else
            {
                node_is_end = equal_points(electrodes[e].end_point, nodes[n]);
                if (node_is_end)
                {
                    condition = 2;
                }
            }
            //================================================
            a = 2*num_electrodes*(e + 1) + e*num_nodes + n;
            b = ld + a;
            c = ld_2 + n*(2*num_electrodes + num_nodes) + e;
            d = num_electrodes + c;
            if (condition == 1)
            {
                we[a] = -1.0; //A
                we[b] = -0.5; //B
                we[c] = 1.0; //C
                we[d] = 0.0; //D
            }
            else if (condition == 2)
            {
                we[a] = 1.0; //A
                we[b] = -0.5; //B
                we[c] = 0.0; //C
                we[d] = 1.0; //D
            }
            else
            {
                we[a] = 0.0; //A
                we[b] = 0.0; //B
                we[c] = 0.0; //C
                we[d] = 0.0; //D
            }
            //================================================
            if (condition > 0)
            {
                no_incidence = 0; //false
            }
        }
        if (no_incidence)
        {
            printf("No electrode is connected to node[%i]\n", n);
            return 1;
        }
    }
    return 0;
}

int fill_impedance(_Complex double* we, const Electrode* electrodes,
                   size_t num_electrodes, size_t num_nodes,
                   const _Complex double* zl, const _Complex double* zt,
                   const _Complex double* yn)
{
    size_t i, k;
    size_t row_size = num_electrodes*2 + num_nodes;
    size_t ld = row_size*num_electrodes;
    for (i = 0; i < num_electrodes; i++)
    {
        for (k = 0; k < num_electrodes; k++)
        {
            we[i*row_size + k] = zl[i*num_electrodes + k]/2.0;
            we[i*row_size + num_electrodes + k] = -zl[i*num_electrodes + k]/2.0;
            we[ld + i*row_size + k] = zt[i*num_electrodes + k];
            we[ld + i*row_size + num_electrodes + k] = zt[i*num_electrodes + k];
        }
    }
    size_t ld_2 = ld*2 + 2*num_electrodes;
    for (i = 0; i < num_nodes; i++)
    {
        for (k = 0; k < num_nodes; k++)
        {
            we[ld_2 + i*row_size + k] = yn[i*num_nodes + k];
        }
    }
    return 0;
}

int solve_electrodes(_Complex double* we, _Complex double* ie,
                     size_t num_electrodes, size_t num_nodes)
{
    MKL_INT n = num_electrodes*2 + num_nodes;
    MKL_INT ipiv[n]; //pivot indices
    int info = LAPACKE_zgesv(LAPACK_ROW_MAJOR, n, 1, we, n, ipiv, ie, 1);
    /* Check for the exact singularity */
    if(info > 0)
    {
        printf("The diagonal element of the triangular factor of WE,\n");
        printf("U(%i,%i) is zero, so that WE is singular;\n", info, info);
        printf("the solution could not be computed.\n");
        exit(info);
    }
    return info;
}

// Alternative approach using equivalent Ynodal
int incidence_alt(double* a, double* b, const Electrode* electrodes,
                  size_t num_electrodes, const double nodes[][3],
                  size_t num_nodes)
{
    int e, n, condition, node_is_start, node_is_end, no_incidence, pos;
    //TODO use a more efficient search or fill-up method if possible
    for (n = 0; n < num_nodes; n++)
    {
        no_incidence = 1;
        for (e = 0; e < num_electrodes; e++)
        {
            node_is_start = equal_points(electrodes[e].start_point, nodes[n]);
            node_is_end = 0;
            /* it is assumed that start_point and end_point of an electrode
            are never equal */
            condition = 0;
            if (node_is_start)
            {
                condition = 1;
            }
            else
            {
                node_is_end = equal_points(electrodes[e].end_point, nodes[n]);
                if (node_is_end)
                {
                    condition = 2;
                }
            }
            //================================================
            pos = num_nodes*e + n;
            if (condition == 1)
            {
                a[pos] = 0.5; //A
                b[pos] = 1.0; //B
            }
            else if (condition == 2)
            {
                a[pos] = 0.5; //A
                b[pos] = -1.0; //B
            }
            else
            {
                a[pos] = 0.0; //A
                b[pos] = 0.0; //B
            }
            //================================================
            if (condition > 0)
            {
                no_incidence = 0; //false
            }
        }
        if (no_incidence)
        {
            printf("No electrode is connected to node[%i]\n", n);
            return 1;
        }
    }
    return 0;
}

int ynodal_eq(_Complex double* yn, const double* a, const double* b,
	          _Complex double* zl, _Complex double* zt, size_t num_electrodes,
              size_t num_nodes)
{
    // TODO use COLUMN_MAJOR matrices, to avoid LAPACKE doing the matrix
    // copying when they are ROW_MAJOR
    // FIXME copy arrays a and b into complex arr, else the zgemm will give
    // wrong results.
    //Force a and b to be (_Complex double*) arguments instead of (double*)?
    lapack_complex_double* arr = malloc(
                   sizeof(lapack_complex_double)*(num_electrodes*num_nodes));

    // yn = aT*(zt^-1)*a + bT*(zl^-1)*b
    lapack_int* ipiv = malloc(sizeof(lapack_int)*num_electrodes);
    // invert zl and zt taking advantage of its symmetry
    /* inv(ZT) is "exploding" (~10e280)
    LAPACKE_zsytrf(LAPACK_ROW_MAJOR, 'U', num_electrodes, zt, num_electrodes, ipiv);
    LAPACKE_zsytri(LAPACK_ROW_MAJOR, 'U', num_electrodes, zt, num_electrodes, ipiv);
    LAPACKE_zsytrf(LAPACK_ROW_MAJOR, 'U', num_electrodes, zl, num_electrodes, ipiv);
    LAPACKE_zsytri(LAPACK_ROW_MAJOR, 'U', num_electrodes, zl, num_electrodes, ipiv);*/
    // invert zl and zt by general matrix LU factorization
    // as both zl and zt are symmetric, it does not matter if COL or ROW MAJOR.
    // Use COL Major so LAPACKE does not transpose them.
    LAPACKE_zgetrf(LAPACK_COL_MAJOR, num_electrodes, num_electrodes,
                   zt, num_electrodes, ipiv);
    LAPACKE_zgetri(LAPACK_COL_MAJOR, num_electrodes, zt, num_electrodes, ipiv);
    LAPACKE_zgetrf(LAPACK_COL_MAJOR, num_electrodes, num_electrodes,
                   zl, num_electrodes, ipiv);
    LAPACKE_zgetri(LAPACK_COL_MAJOR, num_electrodes, zl, num_electrodes, ipiv);
    free(ipiv);
    const double alpha = 1.0;
    const double beta = 0.0;
    lapack_complex_double* c = malloc(
                   sizeof(lapack_complex_double)*(num_electrodes*num_nodes));
    // c = yt*a
    for (size_t i = 0; i < num_electrodes*num_nodes; i++)
    {
        arr[i] = (lapack_complex_double) a[i];
    }
    /*cblas_zsymm(CblasRowMajor, CblasLeft, CblasUpper,
                num_electrodes, num_nodes,
                &alpha, zt, num_electrodes, a, num_nodes,
                &beta, c, num_nodes);*/
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                num_electrodes, num_nodes, num_electrodes,
                &alpha, zt, num_electrodes, arr, num_nodes,
                &beta, c, num_nodes);
    // yn = aT*c
    cblas_zgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                num_nodes, num_nodes, num_electrodes,
                &alpha, arr, num_nodes, c, num_nodes,
                &beta, yn, num_nodes);
    // c = yl*b
    for (size_t i = 0; i < num_electrodes*num_nodes; i++)
    {
        arr[i] = (lapack_complex_double) b[i];
    }
    /*cblas_zsymm(CblasRowMajor, CblasLeft, CblasUpper,
                num_electrodes, num_nodes,
                &alpha, zl, num_electrodes, b, num_nodes,
                &beta, c, num_nodes);*/
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                num_electrodes, num_nodes, num_electrodes,
                &alpha, zl, num_electrodes, arr, num_nodes,
                &beta, c, num_nodes);
    // yn = bT*c + yn
    cblas_zgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                num_nodes, num_nodes, num_electrodes,
                &alpha, arr, num_nodes, c, num_nodes,
                &alpha, yn, num_nodes);
    free(c);
    free(arr);
    return 0;
}

int harmonic_impedance1(size_t ns, const _Complex double* s,
                        const _Complex double* kappa1,
                        const _Complex double* kappa2,
                        const _Complex double* gamma1,
                        const Electrode* electrodes, const Electrode* images,
                        size_t num_electrodes, const double nodes[][3],
                        size_t num_nodes, size_t max_eval, double req_abs_error,
                        double req_rel_error, int error_norm, double rsource,
                        _Complex double* zh)
{
    //TODO extract nodes from electrodes instead of receiving them as an argument?
    //build system to be solved
    double mur1 = 1.0; //TODO take mur1 an argument
    size_t ne2 = num_electrodes*num_electrodes;
    size_t nn2 = num_nodes*num_nodes;
    size_t ss1 = (2*num_electrodes + num_nodes);
    size_t ss2 = ss1*ss1;
    _Complex double zinternal, ref_l, ref_t;
    _Complex double* zl = malloc(sizeof(_Complex double)*ne2);
    _Complex double* zt = malloc(sizeof(_Complex double)*ne2);
    _Complex double* yn = calloc(nn2, sizeof(_Complex double));
    _Complex double* ie = malloc(sizeof(_Complex double)*ss1);
    _Complex double* we = malloc(sizeof(_Complex double)*ss2);
    _Complex double* we_incidence = malloc(sizeof(_Complex double)*ss2);
    //yn[0] = 1.0/rsource; FIXME why did I comment this line?
    fill_incidence(we_incidence, electrodes, num_electrodes, nodes, num_nodes);
    // solve for each frequency: WE*VE = IE
    for (size_t i = 0; i < ns; i++)
    {
        //reset IE
        for (size_t k = 0; k < ss1; k++)
        {
            ie[k] = 0.0;
        }
        ie[ss1 - num_nodes] = 1.0;
        ref_l = (kappa1[i] - kappa2[i])/(kappa1[i] + kappa2[i]);
        ref_t = ref_l;
        //TODO specialized impedance calculation taking advantage of symmetry
        calculate_impedances(electrodes, num_electrodes, zl, zt, gamma1[i],
                             s[i], mur1, kappa1[i], max_eval, req_abs_error,
                             req_rel_error, error_norm, INTG_DOUBLE);
        for (size_t k = 0; k < num_electrodes; k++)
        {
            //FIXME crashing when interfaced with Mathematica because it does not see zbesi_
            zinternal = internal_impedance(s[i], RHO_CU, electrodes[k].radius, 1.0)
                        *electrodes[k].length;
            //zinternal = 0.0;
            zl[k*num_electrodes + k] += zinternal;
        }
        impedances_images(electrodes, images, num_electrodes, zl, zt, gamma1[i],
                          s[i], mur1, kappa1[i], ref_l, ref_t, max_eval,
                          req_abs_error, req_rel_error, error_norm, INTG_DOUBLE);
        //The matrices are pivoted in-place. To avoid overwriting them, copy.
        //That way the filling of the incidence, which is costly, needs to be
        //done only once
        copy_array(we_incidence, we, ss2);
        fill_impedance(we, electrodes, num_electrodes, num_nodes, zl, zt, yn);
        solve_electrodes(we, ie, num_electrodes, num_nodes);
        zh[i] = ie[ss1 - num_nodes];
    }
    free(zl);
    free(zt);
    free(yn);
    free(ie);
    free(we);
    free(we_incidence);
    return 0;
}

int harmonic_impedance1_alt(size_t ns, const _Complex double* s,
                            const _Complex double* kappa1,
                            const _Complex double* kappa2,
                            const _Complex double* gamma1,
                            const Electrode* electrodes,
                            const Electrode* images, size_t num_electrodes,
                            const double nodes[][3], size_t num_nodes,
                            size_t max_eval, double req_abs_error,
                            double req_rel_error, int error_norm,
                            double rsource, _Complex double* zh)
{
    //TODO extract nodes from electrodes instead of receiving it as an argument?
    //build system to be solved
    double mur1 = 1.0; //TODO take mur1 as argument
    size_t ne2 = num_electrodes*num_electrodes;
    size_t nn2 = num_nodes*num_nodes;
    _Complex double zinternal, ref_l, ref_t;
    _Complex double* zl = malloc(sizeof(_Complex double)*ne2);
    _Complex double* zt = malloc(sizeof(_Complex double)*ne2);
    _Complex double* yn = calloc(nn2, sizeof(_Complex double));
    _Complex double* ie = malloc(sizeof(_Complex double)*num_nodes);
    double* a = malloc(sizeof(double)*(num_electrodes*num_nodes));
    double* b = malloc(sizeof(double)*(num_electrodes*num_nodes));
    incidence_alt(a, b, electrodes, num_electrodes, nodes, num_nodes);
    //yn[0] = 1.0/rsource; FIXME why did I comment this line?
    // solve for each frequency: YN*VN = IN
    MKL_INT n = num_electrodes*2 + num_nodes;
    MKL_INT ipiv[n]; //pivot indices
    int info;
    for (size_t i = 0; i < ns; i++)
    {
        //reset IN
        ie[0] = 1.0;
        for (size_t k = 1; k < num_nodes; k++)
        {
            ie[k] = 0.0;
        }
        ref_l = (kappa1[i] - kappa2[i])/(kappa1[i] + kappa2[i]);
        ref_t = ref_l;
        //TODO specialized impedance calculation taking advantage of symmetry
        calculate_impedances(electrodes, num_electrodes, zl, zt, gamma1[i],
                             s[i], mur1, kappa1[i], max_eval, req_abs_error,
                             req_rel_error, error_norm, INTG_DOUBLE);
        for (size_t k = 0; k < num_electrodes; k++)
        {
            zinternal = internal_impedance(s[i], RHO_CU, electrodes[k].radius, 1.0)
                        *electrodes[k].length;
            //zinternal = 0.0;
            zl[k*num_electrodes + k] += zinternal;
        }
        impedances_images(electrodes, images, num_electrodes, zl, zt, gamma1[i],
                          s[i], mur1, kappa1[i], ref_l, ref_t, max_eval,
                          req_abs_error, req_rel_error, error_norm, INTG_DOUBLE);
        ynodal_eq(yn, a, b, zl, zt, num_electrodes, num_nodes);
        info = LAPACKE_zgesv(LAPACK_ROW_MAJOR, num_nodes, 1, yn, num_nodes, ipiv, ie, 1);
        // Check for the exact singularity
        if(info > 0)
        {
            printf("The diagonal element of the triangular factor of YN,\n");
            printf("U(%i,%i) is zero, so that YN is singular;\n", info, info);
            printf("the solution could not be computed.\n");
            exit(info);
        }
        zh[i] = ie[0];
    }
    free(zl);
    free(zt);
    free(yn);
    free(ie);
    free(a);
    free(b);
    return 0;
}
