#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <omp.h>
#include "electrode.h"
#include "auxiliary.h"
#include "cubature.h"
#include "linalg.h"

int
fill_incidence_imm (_Complex double *we, const Electrode *electrodes,
                    size_t num_electrodes, const double *nodes,
                    size_t num_nodes)
{
    int condition, node_is_start, node_is_end, a, b, c, d, no_incidence;
    size_t ld = (2*num_electrodes + num_nodes);
    //TODO use a more efficient search or fill-up method if possible
    for (size_t n = 0; n < num_nodes; n++) {
        no_incidence = 1;
        for (size_t e = 0; e < num_electrodes; e++) {
            node_is_start = equal_points(electrodes[e].start_point, nodes + n*3);
            node_is_end = 0;
            /* it is assumed that start_point and end_point of an electrode
            are never equal */
            condition = 0;
            if (node_is_start) {
                condition = 1;
            } else {
                node_is_end = equal_points(electrodes[e].end_point, nodes + n*3);
                if (node_is_end) condition = 2;
            }
            //================================================
            a = (e + num_nodes) + ld*n;
            b = (e + num_nodes + num_electrodes) + ld*n;
            c = n + ld*(e + num_nodes);
            d = n + ld*(e + num_nodes + num_electrodes);
            if (condition == 1) {
                we[a] = -1.0; //A
                we[b] = -0.5; //B
                we[c] = 1.0; //C
                we[d] = 0.0; //D
            } else if (condition == 2) {
                we[a] = 1.0; //A
                we[b] = -0.5; //B
                we[c] = 0.0; //C
                we[d] = 1.0; //D
            } else {
                we[a] = 0.0; //A
                we[b] = 0.0; //B
                we[c] = 0.0; //C
                we[d] = 0.0; //D
            }
            //================================================
            if (condition > 0) no_incidence = 0; //false
        }
        if (no_incidence) {
            printf("No electrode is connected to node[%i]\n", (int) n);
            return 1;
        }
    }
    return 0;
}

int
fill_impedance_imm (_Complex double *we, size_t num_electrodes, size_t num_nodes,
                    const _Complex double *zl, const _Complex double *zt,
                    const _Complex double *ye)
{
    _Complex double zl2;
    size_t n = num_nodes;
    size_t m = num_electrodes;
    size_t ld = n + 2*m;
    for (size_t k = 0; k < num_electrodes; k++) {
        for (size_t i = 0; i < num_electrodes; i++) {
            zl2 = zl[i + m*k]/2.0;
            we[(i + n) + ld*(k + n)] = zl2;
            we[(i + n) + ld*(k + n + m)] = -zl2;
            we[(i + n + m) + ld*(k + n)] = zt[i + m*k];
            we[(i + n + m) + ld*(k + n + m)] = zt[i + m*k];
        }
    }
    for (size_t k = 0; k < num_nodes; k++) {
        for (size_t i = 0; i < num_nodes; i++) {
            we[i + ld*k] = ye[i + n*k];
        }
    }
    return 0;
}

int
solve_immittance (_Complex double *we, _Complex double *ie,
                  size_t num_electrodes, size_t num_nodes)
{
    int n = num_electrodes*2 + num_nodes;
    int ipiv[n]; //pivot indices
    int info;
    int nrhs = 1;
    zgesv_(&n, &nrhs, we, &n, ipiv, ie, &n, &info);
    //print_matrix("IE", num_nodes, 1, ie, 1);
    // Check for the exact singularity
    if (info > 0) {
        printf("The diagonal element of the triangular factor of WE,\n");
        printf("U(%i,%i) is zero, so that WE is singular;\n", info, info);
        printf("the solution could not be computed.\n");
    }
    return info;
}


int
fill_incidence_adm (double *a, double *b, const Electrode *electrodes,
                    size_t num_electrodes, const double *nodes,
                    size_t num_nodes)
{
    int condition, node_is_start, node_is_end, no_incidence, pos;
    //TODO use a more efficient search or fill-up method if possible
    for (size_t n = 0; n < num_nodes; n++) {
        no_incidence = 1;
        for (size_t e = 0; e < num_electrodes; e++) {
            node_is_start = equal_points(electrodes[e].start_point, nodes + n*3);
            node_is_end = 0;
            /* it is assumed that start_point and end_point of an electrode
            are never equal */
            condition = 0;
            if (node_is_start) {
                condition = 1;
            } else {
                node_is_end = equal_points(electrodes[e].end_point, nodes + n*3);
                if (node_is_end) condition = 2;
            }
            //================================================
            pos = num_nodes*e + n;
            if (condition == 1) {
                a[pos] = 1.0;
                b[pos] = 0.5;
            } else if (condition == 2) {
                a[pos] = -1.0;
                b[pos] = 0.5;
            } else {
                a[pos] = 0.0;
                b[pos] = 0.0;
            }
            //================================================
            if (condition > 0) no_incidence = 0; //false
        }
        if (no_incidence) {
            printf("No electrode is connected to node[%i]\n", (int) n);
            return 1;
        }
    }
    return 0;
}

int
fill_impedance_adm (_Complex double *yn, const double *a, const double *b,
                    _Complex double *zl, _Complex double *zt, size_t num_electrodes,
                    size_t num_nodes, const _Complex double *ye)
{
    int ne = (int) num_electrodes;
    int nn = (int) num_nodes;
    _Complex double *arr = malloc((ne*nn)*sizeof(_Complex double));
    // yn = aT*(zl^-1)*a + bT*(zt^-1)*b
    int *ipiv = malloc(ne * sizeof(int));
    int info;
    /*
    char uplo = 'U';
    int lwmax = 100;
    _Complex double *work = malloc(lwmax * sizeof(_Complex double));
    // inv(ZL)
    //Query the optimal workspace.
    int lwork = -1;
    zsytrf_(&uplo, &ne, zl, &ne, ipiv, work, &lwork, &info);
    lwork = creal(work[0]);
    if (lwork > lwmax) {lwork = lwmax;}
    zsytrf_(&uplo, &ne, zl, &ne, ipiv, work, &lwork, &info);
    zsytri_(&uplo, &ne, zl, &ne, ipiv, work, &info);
    // inv(ZT)
    lwork = -1;
    zsytrf_(&uplo, &ne, zt, &ne, ipiv, work, &lwork, &info);
    lwork = creal(work[0]);
    if (lwork > lwmax) {lwork = lwmax;}
    zsytrf_(&uplo, &ne, zt, &ne, ipiv, work, &lwork, &info);
    zsytri_(&uplo, &ne, zt, &ne, ipiv, work, &info);
    */
    _Complex double *work = malloc(100 * sizeof(_Complex double));
    zgetrf_(&ne, &ne, zt, &ne, ipiv, &info);
    zgetri_(&ne, zt, &ne, ipiv, work, &ne, &info);
    zgetrf_(&ne, &ne, zl, &ne, ipiv, &info);
    zgetri_(&ne, zl, &ne, ipiv, work, &ne, &info);

    _Complex double alpha = 1.0;
    _Complex double beta = 0.0;
    _Complex double *c = malloc((num_electrodes*num_nodes)*sizeof(_Complex double));
    // c = yl*a
    for (size_t i = 0; i < num_electrodes*num_nodes; i++) {
        arr[i] = (_Complex double) a[i];
    }
    char notrans = 'N';
    char trans = 'T';
    //char side = 'U';
    zgemm_(&notrans, &notrans, &ne, &nn, &ne, &alpha, zl, &ne, arr, &nn,
            &beta, c, &nn);
    //zsymm_(&side, &uplo, &ne, &nn, &alpha, zt, &ne, arr, &nn, &beta, c, &nn);
    // yn = aT*c
    zgemm_(&trans, &notrans, &nn, &nn, &ne, &alpha, arr, &nn, c, &nn, &beta,
           yn, &nn);
    // c = yt*b
    for (size_t i = 0; i < num_electrodes*num_nodes; i++) {
        arr[i] = (_Complex double) b[i];
    }
    zgemm_(&notrans, &notrans, &ne, &nn, &ne, &alpha, zt, &nn, arr, &nn,
            &beta, c, &nn);
    //zsymm_(&side, &uplo, &ne, &nn, &alpha, zl, &ne, arr, &nn, &beta, c, &nn);
    // yn = bT*c + yn
    zgemm_(&trans, &notrans, &nn, &nn, &ne, &alpha, arr, &nn, c, &nn, &alpha,
           yn, &nn);
    free(c);
    free(arr);
    free(ipiv);
    free(work);
    return 0;
}

int
solve_admittance (_Complex double *yn, _Complex double *ic, size_t num_nodes)
{
    int n = (int) num_nodes;
    int ipiv[n]; //pivot indices
    int info;
    int nrhs = 1;
    char uplo = 'U';
    int lwmax = 100;
    _Complex double *work = malloc(lwmax * sizeof(_Complex double));
    //Query the optimal workspace.
    int lwork = -1;
    zsysv_(&uplo, &n, &nrhs, yn, &n, ipiv, ic, &n, work, &lwork, &info);
    lwork = creal(work[0]);
    if (lwork > lwmax) {lwork = lwmax;}
    zsysv_(&uplo, &n, &nrhs, yn, &n, ipiv, ic, &n, work, &lwork, &info);
    // Check for the exact singularity
    if (info > 0) {
        printf("The diagonal element of the triangular factor of YN,\n");
        printf("U(%i,%i) is zero, so that YN is singular;\n", info, info);
        printf("the solution could not be computed.\n");
    }
    return info;
    return 0;
}

int
zh_immittance (size_t ns, const _Complex double *s, double sigma, double epsr,
               double mur, const Electrode *electrodes, const Electrode *images,
               size_t num_electrodes, const double *nodes, size_t num_nodes,
               size_t max_eval, double req_abs_error, double req_rel_error,
               _Complex double *zh)
{
    int ne2 = num_electrodes*num_electrodes;
    int nn2 = num_nodes*num_nodes;
    int ss1 = (2*num_electrodes + num_nodes);
    int ss2 = ss1*ss1;
    int nrhs = num_nodes;
    int lengthie = ss1*nrhs;
    _Complex double *zl = malloc(ne2*sizeof(_Complex double));
    _Complex double *zt = malloc(ne2*sizeof(_Complex double));
    _Complex double *ye = calloc(nn2, sizeof(_Complex double));
    _Complex double *inj = malloc(lengthie*sizeof(_Complex double));
    _Complex double *ie = malloc(lengthie*sizeof(_Complex double));
    _Complex double *incidence = malloc(ss2*sizeof(_Complex double));
    _Complex double *we = malloc(ss2*sizeof(_Complex double));
    _Complex double ref_l, ref_t;
    _Complex double kappa, gamma, kappa_air, zi;
    int ipiv[ss2]; //pivot indices
    int info;
    for (size_t k = 0; k < nrhs; k++) {
        inj[k*ss1 + k] = 1.0;
    }
    fill_incidence_imm(incidence, electrodes, num_electrodes, nodes, num_nodes);
    for (size_t i = 0; i < ns; i++) {
        // The matrices are pivoted in-place. To recover them, copy
        for (size_t k = 0; k < ss2; k++) we[k] = incidence[k];
        for (size_t k = 0; k < lengthie; k++) ie[k] = inj[k];
        kappa = (sigma + s[i]*epsr*EPS0); //soil complex conductivity
        gamma = csqrt(s[i]*mur*MU0*kappa); //soil propagation constant
        kappa_air = (s[i]*EPS0);
        ref_l = (kappa - kappa_air)/(kappa + kappa_air);
        ref_t = ref_l;
        calculate_impedances(electrodes, num_electrodes, zl, zt, gamma, s[i],
                             mur, kappa, max_eval, req_abs_error, req_rel_error,
                             ERROR_PAIRED, INTG_DOUBLE);
        for (size_t k = 0; k < num_electrodes; k++) {
            zi = internal_impedance(s[i], RHO_CU, electrodes[k].radius, 1.0, &info);
            zl[k*num_electrodes + k] += zi*electrodes[k].length + electrodes[k].zi;
        }
        impedances_images(electrodes, images, num_electrodes, zl, zt, gamma,
                          s[i], mur, kappa, ref_l, ref_t, max_eval, req_abs_error,
                          req_rel_error, ERROR_PAIRED, INTG_DOUBLE);
        fill_impedance_imm(we, num_electrodes, num_nodes, zl, zt, ye);
        zgesv_(&ss1, &nrhs, we, &ss1, ipiv, ie, &ss1, &info);
        // Check for the exact singularity
        if (info > 0) {
            printf("The diagonal element of the triangular factor of WE,\n");
            printf("U(%i,%i) is zero, so that WE is singular;\n", info, info);
            printf("the solution could not be computed.\n");
        }
        for (size_t k = 0; k < nrhs; k++) {
            zh[k + i*nrhs] = ie[k*ss1 + k];
        }
    }
    free(zl);
    free(zt);
    free(ye);
    free(inj);
    free(ie);
    free(incidence);
    free(we);
    return 0;
}

int
sim_immittance (size_t ns, const _Complex double *s, double sigma, double epsr,
                double mur, const Electrode *electrodes, const Electrode *images,
                size_t num_electrodes, const double *nodes, size_t num_nodes,
                size_t max_eval, double req_abs_error, double req_rel_error,
                size_t inj_node, const _Complex double *inj_current,
                const _Complex double *inj_adm, _Complex double *u,
                _Complex double *il, _Complex double *it)
{
    int ierr; //error code
    int ne2 = num_electrodes*num_electrodes;
    int nn2 = num_nodes*num_nodes;
    int ss1 = (2*num_electrodes + num_nodes);
    int ss2 = ss1*ss1;
    _Complex double *zl = malloc(ne2*sizeof(_Complex double));
    _Complex double *zt = malloc(ne2*sizeof(_Complex double));
    _Complex double *ye = calloc(nn2, sizeof(_Complex double));
    _Complex double *ie = malloc(ss1*sizeof(_Complex double));
    _Complex double *incidence = malloc(ss2*sizeof(_Complex double));
    _Complex double *we = malloc(ss2*sizeof(_Complex double));
    _Complex double ref_l, ref_t;
    _Complex double kappa, gamma, kappa_air, zi, i1, i2;
    fill_incidence_imm(incidence, electrodes, num_electrodes, nodes, num_nodes);
    for (size_t i = 0; i < ns; i++) {
        // The matrices are pivoted in-place. To recover them, copy
        for (size_t k = 0; k < ss2; k++) we[k] = incidence[k];
        for (size_t k = 0; k < num_nodes; k++) ie[k] = 0.0;
        ie[inj_node] = inj_current[i];
        ye[inj_node*(1 + inj_node)] = inj_adm[i];
        kappa = (sigma + s[i]*epsr*EPS0); //soil complex conductivity
        gamma = csqrt(s[i]*mur*MU0*kappa); //soil propagation constant
        kappa_air = (s[i]*EPS0);
        ref_l = (kappa - kappa_air)/(kappa + kappa_air);
        ref_t = ref_l;
        calculate_impedances(electrodes, num_electrodes, zl, zt, gamma, s[i],
                             mur, kappa, max_eval, req_abs_error, req_rel_error,
                             ERROR_PAIRED, INTG_DOUBLE);
        for (size_t k = 0; k < num_electrodes; k++) {
            zi = internal_impedance(s[i], RHO_CU, electrodes[k].radius, 1.0, &ierr);
            zl[k*num_electrodes + k] += zi*electrodes[k].length + electrodes[k].zi;
        }
        impedances_images(electrodes, images, num_electrodes, zl, zt, gamma,
                          s[i], mur, kappa, ref_l, ref_t, max_eval, req_abs_error,
                          req_rel_error, ERROR_PAIRED, INTG_DOUBLE);
        fill_impedance_imm(we, num_electrodes, num_nodes, zl, zt, ye);
        solve_immittance(we, ie, num_electrodes, num_nodes);
        for (size_t k = 0; k < num_nodes; k++) {
            u[k + i*num_nodes] = ie[k];
        }
        for (size_t k = 0; k < num_electrodes; k++) {
            i1 = ie[num_nodes + k];
            i2 = ie[num_nodes + k + num_electrodes];
            il[k + i*num_electrodes] = (i1 - i2)/2;
            it[k + i*num_electrodes] = i1 + i2;
        }
    }
    free(zl);
    free(zt);
    free(ye);
    free(ie);
    free(incidence);
    free(we);
    return 0;
}
