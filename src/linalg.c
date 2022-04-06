#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "electrode.h"
#include "auxiliary.h"
#include "cubature.h"
#include "linalg.h"
#include "blas.h"
#include "lapack.h"

int
fill_incidence_imm (_Complex float* wg, const Electrode* electrodes,
                    size_t num_electrodes, const float* nodes,
                    size_t num_nodes)
{
    enum {
        NO_INCIDENCE = 0,
        NODE_IS_START = 1,
        NODE_IS_END = 2
    };
    int a, b, c, d;
    int condition;
    bool no_incidence;
    size_t ld = (2 * num_electrodes + num_nodes);
    float tol = 1e-5;  // using DBL_EPSILON leads to wrong incidence matrix
    long int suma = 0;
    float sumb = 0.0;
    for (size_t n = 0; n < num_nodes; n++) {
        no_incidence = true;
        for (size_t e = 0; e < num_electrodes; e++) {
            // it is assumed that start_point and end_point of an electrode
            // are never equal
            condition = NO_INCIDENCE;
            if ( equal_points_tol(electrodes[e].start_point, nodes + n*3, tol) ) {
                condition = NODE_IS_START;
            } else if ( equal_points_tol(electrodes[e].end_point, nodes + n*3, tol) ) {
                condition = NODE_IS_END;
            }
            // ================================================
            a = (e + num_nodes) + ld * n;
            b = (e + num_nodes + num_electrodes) + ld * n;
            c = n + ld * (e + num_nodes);
            d = n + ld * (e + num_nodes + num_electrodes);
            switch (condition) {
                case NODE_IS_START:
                    wg[a] = 1.0; // A
                    wg[b] = 0.5; // B
                    wg[c] = 1.0; // A^T
                    wg[d] = 0.5; // B^T
                    break;

                case NODE_IS_END:
                    wg[a] = -1.0; // A
                    wg[b] = 0.5; // B
                    wg[c] = -1.0; // A^T
                    wg[d] = 0.5; // B^T
                    break;

                default:  // NO_INCIDENCE
                    wg[a] = 0.0; // A
                    wg[b] = 0.0; // B
                    wg[c] = 0.0; // A^T
                    wg[d] = 0.0; // B^T
                    break;
            }
            suma += wg[a];
            sumb += wg[b];
            // ================================================
            if (condition != NO_INCIDENCE) no_incidence = false;
        }
        if (no_incidence) {
            printf("No electrode is connected to node[%i]\n", (int) n);
            return 1;
        }
    }
    if (suma != 0) {
        printf("incidence matrix A is wrong. sum(A) = %li != 0\n", suma);
        return -suma;
    }
    if (sumb != num_electrodes) {
        printf("incidence matrix B is wrong. sum(B) = %.1f != %.1f\n", sumb, (float) num_electrodes);
        return (sumb - num_electrodes);
    }
    return 0;
}

int
fill_impedance_imm (_Complex float* wg, const _Complex float* zl,
                    const _Complex float* zt, size_t num_electrodes,
                    size_t num_nodes)
{
    size_t n = num_nodes;
    size_t m = num_electrodes;
    size_t ld = n + 2*m;
    for (size_t k = 0; k < num_electrodes; k++) {
        for (size_t i = 0; i < num_electrodes; i++) {
            wg[(i + n) + ld*(k + n)] = -zl[i + m*k];
            wg[(i + n) + ld*(k + n + m)] = 0.0;
            wg[(i + n + m) + ld*(k + n)] = 0.0;
            wg[(i + n + m) + ld*(k + n + m)] = -zt[i + m*k];
        }
    }
    for (size_t k = 0; k < num_nodes; k++) {
        for (size_t i = 0; i < num_nodes; i++) {
            wg[i + ld*k] = 0.0;
        }
    }
    return 0;
}

int
solve_immittance (_Complex float* wg, _Complex float* ie,
                  size_t num_electrodes, size_t num_nodes)
{
    int n = num_electrodes*2 + num_nodes;
    int ipiv[n]; //pivot indices
    int info = 0;
    int nrhs = 1;
    //zgesv_(&n, &nrhs, wg, &n, ipiv, ie, &n, &info);
    char uplo = 'L';
    int lwmax = 100;
    _Complex float* work = malloc(lwmax * sizeof(_Complex float));
    //Query the optimal workspace.
    int lwork = -1;
    csysv_(&uplo, &n, &nrhs, wg, &n, ipiv, ie, &n, work, &lwork, &info);
    lwork = creal(work[0]);
    if (lwork > lwmax) {lwork = lwmax;}
    work = realloc(work, lwmax * sizeof(_Complex float));
    csysv_(&uplo, &n, &nrhs, wg, &n, ipiv, ie, &n, work, &lwork, &info);
    // Check for the exact singularity
    if (info > 0) {
        printf("The diagonal element of the triangular factor of YN,\n");
        printf("U(%i,%i) is zero, so that YN is singular;\n", info, info);
        printf("the solution could not be computed.\n");
    }
    free(work);
    return info;
}


int
fill_incidence_adm (_Complex float* a, _Complex float* b,
                    const Electrode* electrodes, size_t num_electrodes,
                    const float* nodes, size_t num_nodes)
{
    enum {
        NO_INCIDENCE = 0,
        NODE_IS_START = 1,
        NODE_IS_END = 2
    };
    int condition, pos;
    bool no_incidence;
    float tol = 1e-5;  // using DBL_EPSILON leads to wrong incidence matrix
    long int suma = 0;
    float sumb = 0.0;
    for (size_t n = 0; n < num_nodes; n++) {
        no_incidence = true;
        for (size_t e = 0; e < num_electrodes; e++) {
            // it is assumed that start_point and end_point of an electrode
            // are never equal
            condition = NO_INCIDENCE;
            if ( equal_points_tol(electrodes[e].start_point, nodes + n*3, tol) ) {
                condition = NODE_IS_START;
            } else if ( equal_points_tol(electrodes[e].end_point, nodes + n*3, tol) ) {
                condition = NODE_IS_END;
            }
            // ================================================
            pos = num_electrodes * n + e;
            switch (condition) {
                case NODE_IS_START:
                    b[pos] = 0.5;
                    a[pos] = 1.0;
                    break;

                case NODE_IS_END:
                    b[pos] = 0.5;
                    a[pos] = -1.0;
                    break;

                default:  // NO_INCIDENCE
                    a[pos] = 0.0;
                    b[pos] = 0.0;
                    break;
            }
            if (condition != NO_INCIDENCE) no_incidence = false;
            suma += a[pos];
            sumb += b[pos];
        }
        if (no_incidence) {
            printf("No electrode is connected to node[%i]\n", (int) n);
            return 1;
        }
    }
    if (suma != 0) {
        printf("incidence matrix A is wrong. sum(A) = %li != 0\n", suma);
        return -suma;
    }
    if (sumb != num_electrodes) {
        printf("incidence matrix B is wrong. sum(B) = %.1f != %.1f\n", sumb, (float) num_electrodes);
        return (sumb - num_electrodes);
    }
    return 0;
}

int
fill_impedance_adm (_Complex float* yn, _Complex float* zl, _Complex float* zt,
                    _Complex float* a, _Complex float* b, size_t num_electrodes,
                    size_t num_nodes)
{
    // yn = aT*(zl^-1)*a + bT*(zt^-1)*b
    int ne = num_electrodes;
    int nn = num_nodes;
    char uplo = 'L';
    char notrans = 'N';
    char trans = 'T';
    char side = 'L';
    _Complex float one = 1.0;
    _Complex float zero = 0.0;
    int lwork = -1;
    int info;
    _Complex float* c = malloc( (ne*nn) * sizeof(_Complex float) );
    int* ipiv = malloc(ne * sizeof(size_t));
    _Complex float* work = malloc(100 * sizeof(_Complex float));
    // Query the optimal workspace.
    csytrf_(&uplo, &ne, zl, &ne, ipiv, work, &lwork, &info);
    if (info != 0) return info;
    lwork = creal(work[0]);
    work = realloc(work, lwork * sizeof(_Complex float));
    // inv(ZL)
    csytrf_(&uplo, &ne, zl, &ne, ipiv, work, &lwork, &info);
    if (info != 0) return info;
    csytri_(&uplo, &ne, zl, &ne, ipiv, work, &info);
    // inv(ZT)
    csytrf_(&uplo, &ne, zt, &ne, ipiv, work, &lwork, &info);
    if (info != 0) return info;
    csytri_(&uplo, &ne, zt, &ne, ipiv, work, &info);
    // c = yl*a + c*0
    csymm_(&side, &uplo, &ne, &nn, &one, zl, &ne, a, &ne, &zero, c, &ne);
    // yn = aT*c + yn*0
    cgemm_(&trans, &notrans, &nn, &nn, &ne, &one, a, &ne, c, &ne, &zero, yn, &nn);
    // c = yt*b + c*0
    csymm_(&side, &uplo, &ne, &nn, &one, zt, &ne, b, &ne, &zero, c, &ne);
    // yn = bT*c + yn
    cgemm_(&trans, &notrans, &nn, &nn, &ne, &one, b, &ne, c, &ne, &one, yn, &nn);
    // if using Intel MKL, replace the above BLAS calls by:
    /*cblas_zsymm(CblasColMajor, CblasLeft, CblasLower,
                ne, nn, &one, zl, ne, a, ne, &zero, c, ne);
    cblas_zgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                nn, nn, ne, &one, a, ne, c, ne, &zero, yn, nn);
    cblas_zsymm(CblasColMajor, CblasLeft, CblasLower,
                ne, nn, &one, zt, ne, b, ne, &zero, c, ne);
    cblas_zgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                nn, nn, ne, &one, b, ne, c, ne, &one, yn, nn);*/
    free(c);
    free(ipiv);
    free(work);
    return 0;
}

int
calculate_yla_ytb (_Complex float* yla, _Complex float* ytb,
                   _Complex float* zl, _Complex float* zt,
                   _Complex float* a, _Complex float* b,
                   size_t num_electrodes, size_t num_nodes)
{
    int ne = num_electrodes;
    int nn = num_nodes;
    char uplo = 'L';
    char side = 'L';
    _Complex float one = 1.0;
    _Complex float zero = 0.0;
    int lwork = -1;
    int info;
    int* ipiv = malloc(ne * sizeof(size_t));
    _Complex float* work = malloc(100 * sizeof(_Complex float));
    // Query the optimal workspace.
    csytrf_(&uplo, &ne, zl, &ne, ipiv, work, &lwork, &info);
    if (info != 0) return info;
    lwork = creal(work[0]);
    work = realloc(work, lwork * sizeof(_Complex float));
    // inv(ZL)
    csytrf_(&uplo, &ne, zl, &ne, ipiv, work, &lwork, &info);
    if (info != 0) return info;
    csytri_(&uplo, &ne, zl, &ne, ipiv, work, &info);
    // inv(ZT)
    csytrf_(&uplo, &ne, zt, &ne, ipiv, work, &lwork, &info);
    if (info != 0) return info;
    csytri_(&uplo, &ne, zt, &ne, ipiv, work, &info);
    csymm_(&side, &uplo, &ne, &nn, &one, zl, &ne, a, &ne, &zero, yla, &ne);
    csymm_(&side, &uplo, &ne, &nn, &one, zt, &ne, b, &ne, &zero, ytb, &ne);
    // if using Intel MKL, replace the above BLAS calls by:
    /*cblas_zsymm(CblasColMajor, CblasLeft, CblasLower,
                ne, nn, &one, zl, ne, a, ne, &zero, yla, ne);
    cblas_zsymm(CblasColMajor, CblasLeft, CblasLower,
                ne, nn, &one, zt, ne, b, ne, &zero, ytb, ne);*/
    free(ipiv);
    free(work);
    return 0;
}

int
fill_impedance_adm2 (_Complex float* yn, _Complex float* yla,
                     _Complex float* ytb, _Complex float* a,
                     _Complex float* b, size_t num_electrodes,
                     size_t num_nodes)
{
    int ne = num_electrodes;
    int nn = num_nodes;
    char notrans = 'N';
    char trans = 'T';
    _Complex float one = 1.0;
    _Complex float zero = 0.0;
    // yn = aT*(zl^-1)*a + bT*(zt^-1)*b
    // yn = aT*yla + yn*0
    cgemm_(&trans, &notrans, &nn, &nn, &ne, &one, a, &ne, yla, &ne, &zero, yn, &nn);
    // yn = bT*ytb + yn
    cgemm_(&trans, &notrans, &nn, &nn, &ne, &one, b, &ne, ytb, &ne, &one, yn, &nn);
    // if using Intel MKL, replace the above BLAS calls by:
    /*cblas_zgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                nn, nn, ne, &one, a, ne, yla, ne, &zero, yn, nn);
    cblas_zgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                nn, nn, ne, &one, b, ne, ytb, ne, &one, yn, nn);*/
    return 0;
}

int
solve_admittance (_Complex float* yn, _Complex float* ic, size_t num_nodes)
{
    int n = (int) num_nodes;
    int ipiv[n]; //pivot indices
    int info;
    int nrhs = 1;
    char uplo = 'L';
    int lwork = -1;
    _Complex float* work = malloc(100 * sizeof(_Complex float));
    // Query the optimal workspace.
    csysv_(&uplo, &n, &nrhs, yn, &n, ipiv, ic, &n, work, &lwork, &info);
    //if (info != 0) return info;
    lwork = creal(work[0]);
    work = realloc(work, lwork * sizeof(_Complex float));
    csysv_(&uplo, &n, &nrhs, yn, &n, ipiv, ic, &n, work, &lwork, &info);
    // Check for the exact singularity
    if (info > 0) {
        printf("The diagonal element of the triangular factor of YN,\n");
        printf("U(%i,%i) is zero, so that YN is singular;\n", info, info);
        printf("the solution could not be computed.\n");
    }
    return info;
}
