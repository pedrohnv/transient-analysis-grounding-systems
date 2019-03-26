#include "electrode.h"
#include "linalg.h"
#include "mkl.h"
#include "mkl_lapacke.h"
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <omp.h>

int
solve_electrodes (_Complex double *we, _Complex double *ie, size_t num_electrodes,
                  size_t num_nodes)
{
    MKL_INT n = num_electrodes*2 + num_nodes;
    MKL_INT ipiv[n]; //pivot indices
    int info = LAPACKE_zgesv(LAPACK_ROW_MAJOR, n, 1, we, n, ipiv, ie, 1);
    /* Check for the exact singularity */
    if (info > 0) {
        printf("The diagonal element of the triangular factor of WE,\n");
        printf("U(%i,%i) is zero, so that WE is singular;\n", info, info);
        printf("the solution could not be computed.\n");
        exit(info);
    }
    return info;
}

int
ynodal_eq (_Complex double *yn, const double *a, const double *b,
           _Complex double *zl, _Complex double *zt, size_t num_electrodes,
           size_t num_nodes)
{
    // TODO use COLUMN_MAJOR matrices, to avoid LAPACKE doing the matrix
    // copying when they are ROW_MAJOR
    // FIXME copy arrays a and b into complex arr, else the zgemm will give
    // wrong results.
    //Force a and b to be (_Complex double*) arguments instead of (double*)?
    lapack_complex_double *arr = malloc(
                   sizeof(lapack_complex_double)*(num_electrodes*num_nodes));
    // yn = aT*(zt^-1)*a + bT*(zl^-1)*b
    lapack_int *ipiv = malloc(sizeof(lapack_int)*num_electrodes);
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
    lapack_complex_double *c = malloc(
                   sizeof(lapack_complex_double)*(num_electrodes*num_nodes));
    // c = yt*a
    for (size_t i = 0; i < num_electrodes*num_nodes; i++) {
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
    for (size_t i = 0; i < num_electrodes*num_nodes; i++) {
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

int
harmonic_impedance1 (size_t ns, const _Complex double *s,
                     const _Complex double *kappa1,
                     const _Complex double *kappa2,
                     const _Complex double *gamma1,
                     const Electrode *electrodes, const Electrode *images,
                     size_t num_electrodes, const double nodes[][3],
                     size_t num_nodes, size_t max_eval, double req_abs_error,
                     double req_rel_error, int error_norm, double rsource,
                     _Complex double *zh)
{
    //TODO extract nodes from electrodes instead of receiving them as an argument?
    //build system to be solved
    double mur1 = 1.0; //TODO take mur1 an argument
    size_t ne2 = num_electrodes*num_electrodes;
    size_t nn2 = num_nodes*num_nodes;
    size_t ss1 = (2*num_electrodes + num_nodes);
    size_t ss2 = ss1*ss1;
    _Complex double zinternal, ref_l, ref_t;
    _Complex double *zl = malloc(sizeof(_Complex double)*ne2);
    _Complex double *zt = malloc(sizeof(_Complex double)*ne2);
    _Complex double *yn = calloc(nn2, sizeof(_Complex double));
    _Complex double *ie = malloc(sizeof(_Complex double)*ss1);
    _Complex double *we = malloc(sizeof(_Complex double)*ss2);
    _Complex double *we_incidence = malloc(sizeof(_Complex double)*ss2);
    //yn[0] = 1.0/rsource; FIXME why did I comment this line?
    fill_incidence(we_incidence, electrodes, num_electrodes, nodes, num_nodes);
    // solve for each frequency: WE*VE = IE
    for (size_t i = 0; i < ns; i++) {
        //reset IE
        for (size_t k = 0; k < ss1; k++) ie[k] = 0.0;
        ie[ss1 - num_nodes] = 1.0;
        ref_l = (kappa1[i] - kappa2[i])/(kappa1[i] + kappa2[i]);
        ref_t = ref_l;
        //TODO specialized impedance calculation taking advantage of symmetry
        calculate_impedances(electrodes, num_electrodes, zl, zt, gamma1[i],
                             s[i], mur1, kappa1[i], max_eval, req_abs_error,
                             req_rel_error, error_norm, INTG_DOUBLE);
        for (size_t k = 0; k < num_electrodes; k++) {
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
        for (size_t i = 0; i < ss2; i++) {
            we[i] = we_incidence[i];
        }
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

int
harmonic_impedance1 (size_t ns, const _Complex double *s,
                     const _Complex double *kappa1,
                     const _Complex double *kappa2,
                     const _Complex double *gamma1,
                     const Electrode *electrodes, const Electrode *images,
                     size_t num_electrodes, const double nodes[][3],
                     size_t num_nodes, size_t max_eval, double req_abs_error,
                     double req_rel_error, int error_norm, double rsource,
                     _Complex double *zh)
{
    //TODO extract nodes from electrodes instead of receiving it as an argument?
    //build system to be solved
    double mur1 = 1.0; //TODO take mur1 as argument
    size_t ne2 = num_electrodes*num_electrodes;
    size_t nn2 = num_nodes*num_nodes;
    _Complex double zinternal, ref_l, ref_t;
    _Complex double *zl = malloc(sizeof(_Complex double)*ne2);
    _Complex double *zt = malloc(sizeof(_Complex double)*ne2);
    _Complex double *yn = calloc(nn2, sizeof(_Complex double));
    _Complex double *ie = malloc(sizeof(_Complex double)*num_nodes);
    double *a = malloc(sizeof(double)*(num_electrodes*num_nodes));
    double *b = malloc(sizeof(double)*(num_electrodes*num_nodes));
    incidence_alt(a, b, electrodes, num_electrodes, nodes, num_nodes);
    //yn[0] = 1.0/rsource; FIXME why did I comment this line?
    // solve for each frequency: YN*VN = IN
    MKL_INT n = num_electrodes*2 + num_nodes;
    MKL_INT ipiv[n]; //pivot indices
    int info;
    for (size_t i = 0; i < ns; i++) {
        //reset IN
        ie[0] = 1.0;
        for (size_t k = 1; k < num_nodes; k++) ie[k] = 0.0;
        ref_l = (kappa1[i] - kappa2[i])/(kappa1[i] + kappa2[i]);
        ref_t = ref_l;
        //TODO specialized impedance calculation taking advantage of symmetry
        calculate_impedances(electrodes, num_electrodes, zl, zt, gamma1[i],
                             s[i], mur1, kappa1[i], max_eval, req_abs_error,
                             req_rel_error, error_norm, INTG_DOUBLE);
        for (size_t k = 0; k < num_electrodes; k++) {
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
        if(info > 0) {
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
