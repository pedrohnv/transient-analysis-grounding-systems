/*
Calculates the harmonic impedance reading the electrodes and frequencies from files.
Uses mHEM formulation considering a Alipio-Visacro soil model with mean values.

Soil is considered to be at z=0.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include <omp.h>
#include "auxiliary.h"
#include "electrode.h"
#include "linalg.h"
#include "cubature.h"
#include "grid.h"

/**
Runs the simulations.

Parameters
----------
    sigma0 : soil conductivity in low frequency, in S
    freq : array of frequencies, in Hz
    ns : number of frequencies
    inj_node : injection node number
    ne : number of electrodes
    nn : number of nodes
*/
int
run_case (double sigma0, double* freq, size_t ns, size_t inj_node, size_t ne, size_t nn)
{
      const unsigned int NRHS = 1;
    // numerical integration parameters
    size_t max_eval = 0;
    double req_abs_error = 1e-6;
    double req_rel_error = 1e-4;

    // soil model (Alipio) parameters
    double mur = 1.0;  // soil rel. magnetic permeability
    // parameters that I fitted:
    double h_soil = 1.26 * pow(sigma0 * 1e3, -0.73);
    double g_soil = 0.54;
    double eps_ratio = 12;  // soil rel. permitivitty ratio

    // frequencies
    _Complex double *s = malloc(ns * sizeof(_Complex double));
    for (int k = 0; k < ns; k++) {
        s[k] = freq[k] * TWO_PI * I;
    }
    // electrodes definition ===================================================
    Electrode* electrodes = malloc(sizeof(Electrode) * ne);
    electrodes_file("electrodes.csv", electrodes, ne);
    double* nodes = malloc(nn * 3 * sizeof(double));
    nodes_file("nodes.csv", nodes, nn);
    Electrode* images = malloc(sizeof(Electrode) * ne);
    for (size_t m = 0; m < ne; m++) {
        populate_electrode(images + m, electrodes[m].start_point,
                           electrodes[m].end_point, electrodes[m].radius);
        images[m].start_point[2] = -images[m].start_point[2];
        images[m].end_point[2] = -images[m].end_point[2];
        images[m].middle_point[2] = -images[m].middle_point[2];
    }
    // malloc matrices =========================================================
    size_t ne2 = ne * ne;
    size_t nn2 = nn * nn;
    _Complex double* potzl = malloc(ne2 * sizeof(_Complex double));
    _Complex double* potzt = malloc(ne2 * sizeof(_Complex double));
    _Complex double* potzli = malloc(ne2 * sizeof(_Complex double));
    _Complex double* potzti = malloc(ne2 * sizeof(_Complex double));
    _Complex double* a = malloc((ne * nn) * sizeof(_Complex double));
    _Complex double* b = malloc((ne * nn) * sizeof(_Complex double));
    _Complex double* zh = malloc(NRHS * ns * sizeof(_Complex double));
    int err;
    err = fill_incidence_adm(a, b, electrodes, ne, nodes, nn);
    if (err != 0) printf("Could not build incidence matrices\n");
    err = calculate_impedances(potzl, potzt, electrodes, ne, 0.0, 0.0, 0.0, 0.0,
                               max_eval, req_abs_error, req_rel_error, INTG_MHEM);
    if (err != 0) printf("integration error\n");
    err = impedances_images(potzli, potzti, electrodes, images, ne,
                            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, max_eval,
                            req_abs_error, req_rel_error, INTG_MHEM);
    if (err != 0) printf("integration error\n");
    double begin = omp_get_wtime();  // to estimate time until completion
    #pragma omp parallel private(err)
    {
        #pragma omp single
        {
            printf("avaible threads: %i\n", omp_get_num_threads());
        }
        _Complex double* zl = malloc(ne2 * sizeof(_Complex double));
        _Complex double* zt = malloc(ne2 * sizeof(_Complex double));
        _Complex double* yn = malloc(nn2 * sizeof(_Complex double));
        _Complex double* ie = malloc(NRHS * nn * sizeof(_Complex double));
        _Complex double sigma, epsr, kappa, gamma;  // soil parameters
        _Complex double ref_l, ref_t;  // reflection coefficients (images)
        _Complex double iwu_4pi, one_4pik, exp_gr;
        double rbar;
        // use raw BLAS and LAPACK to solve system with multiple RHS
        int nn1 = (int) nn;
        int ldie = nn1;  // leading dimension of ie
        int* ipiv = malloc(nn1 * sizeof(int));  // pivot indices
        int info;
        int nrhs = NRHS;
        char uplo = 'L';  // matrices are symmetric, only Lower Half is set
        int lwork = -1;  // signal to Query the optimal workspace
        _Complex double wkopt;
        zsysv_(&uplo, &nn1, &nrhs, yn, &nn1, ipiv, ie, &ldie, &wkopt, &lwork, &info);
        lwork = creal(wkopt);
        _Complex double* work = malloc(lwork * sizeof(_Complex double));
        #pragma omp for
        for (size_t i = 0; i < ns; i++) {
            //printf("i = %li from thread %d\n", i, omp_get_thread_num());
            alipio_soil(&sigma, &epsr, sigma0, s[i], h_soil, g_soil, eps_ratio);
            kappa = (sigma + s[i] * epsr * EPS0);  // soil complex conductivity
            gamma = csqrt(s[i] * MU0 * kappa);  // soil propagation constant
            iwu_4pi = s[i] * mur * MU0 / (FOUR_PI);
            one_4pik = 1.0 / (FOUR_PI * kappa);
            // reflection coefficient, soil / air
            ref_t = (kappa - s[i] * EPS0) / (kappa + s[i] * EPS0);
            ref_l = 1.0;  // longitudinal current has +1 image
            // modified HEM (mHEM):
            for (size_t m = 0; m < ne; m++) {
                for (size_t k = m; k < ne; k++) {
                    rbar = vector_length(electrodes[k].middle_point,
                                         electrodes[m].middle_point);
                    exp_gr = cexp(-gamma * rbar);
                    zl[m * ne + k] = exp_gr * potzl[m * ne + k];
                    zt[m * ne + k] = exp_gr * potzt[m * ne + k];
                    rbar = vector_length(electrodes[k].middle_point,
                                         images[m].middle_point);
                    exp_gr = cexp(-gamma * rbar);
                    zl[m * ne + k] += ref_l * exp_gr * potzli[m * ne + k];
                    zt[m * ne + k] += ref_t * exp_gr * potzti[m * ne + k];
                    zl[m * ne + k] *= iwu_4pi;
                    zt[m * ne + k] *= one_4pik;
                }
            }
            // Traditional HEM (highly discouraged):
            /*calculate_impedances(zl, zt, electrodes, ne, gamma, s[i], 1.0, kappa,
                                 max_eval, req_abs_error, req_rel_error, INTG_DOUBLE);
            impedances_images(zl, zt, electrodes, images, ne, gamma, s[i], mur,
                              kappa, ref_l, ref_t, max_eval, req_abs_error,
                              req_rel_error, INTG_DOUBLE);*/
            fill_impedance_adm(yn, zl, zt, a, b, ne, nn);
            for (size_t m = 0; m < (NRHS * nn); m++) {
                ie[m] = 0.0;
            }
            ie[inj_node] = 1.0;
            // solve
            zsysv_(&uplo, &nn1, &nrhs, yn, &nn1, ipiv, ie, &ldie, work, &lwork, &info);
            // Check for the exact singularity
            if (info > 0) {
                printf("The diagonal element of the triangular factor of YN,\n");
                printf("U(%i,%i) is zero, so that YN is singular;\n", info, info);
                printf("the solution could not be computed.\n");
                printf("  in the %li-th frequency\n", i);
                exit(info);
            } else if (info < 0) {
                printf("the %i-th parameter to zsysv had an illegal value.\n", info);
            }
            zh[i] = ie[inj_node];
            if (i == 0) {
                printf("Expected more time until completion of the frequency loop: %.2f min.\n",
                       (omp_get_wtime() - begin) * ns / 60.0 / omp_get_num_threads());
            }
        }
        free(zl);
        free(zt);
        free(yn);
        free(ie);
        free(ipiv);
        free(work);
    } // end parallel
    printf("Frequency loop ended. Saving results.\n");
    char zh_file_name[60];
    sprintf(zh_file_name, "zh.csv");
    FILE* zh_file = fopen(zh_file_name, "w");
    for (size_t i = 0; i < ns; i++) {
        fprintf(zh_file, "%e%+e%s\n", creal(zh[i]), cimag(zh[i]), "im");
    }
    fclose(zh_file);
    free(s);
    free(potzl);
    free(potzt);
    free(potzli);
    free(potzti);
    free(a);
    free(b);
    free(zh);
    free(electrodes);
    free(images);
    free(nodes);
    return 0;
}

/** Returns the number of lines in file. */
int count_lines(FILE* fileptr)
{
      int ns = 1;
    char chr = getc(fileptr);
    while (chr != EOF) {
        if (chr == '\n') {
            ns++;
        }
        chr = getc(fileptr);
    }
    return ns;
}

int
main (int argc, char *argv[])
{
    if (argc != 3) {
        printf("Wrong number of arguments, the following are needed:\n");
        printf("  rho0 : soil low frequency resistivity in (ohm*m)\n");
        printf("  inj_node : the injection node number (from 1 to N)\n");
        exit(argc);
    }
    char *p;
    double sigma0 = 1.0 / strtod(argv[1], &p);
    int inj_node = atoi(argv[2]) - 1;
    double start_time = omp_get_wtime();
    // read files
    char fname_electrodes[] = "electrodes.csv";
    FILE* file_electrodes = fopen(fname_electrodes, "r");
    if (file_electrodes == NULL) {
        printf("Cannot open file %s\n", fname_electrodes);
        return -11;
    }
    size_t ne = count_lines(file_electrodes);
    printf("num. electrodes: %li\n", ne);
    fclose(file_electrodes);
    char fname_nodes[] = "nodes.csv";
    FILE* file_nodes = fopen(fname_nodes, "r");
    if (file_nodes == NULL) {
        printf("Cannot open file %s\n", fname_nodes);
        return -12;
    }
    size_t nn = count_lines(file_nodes);
    printf("num. nodes: %li\n", nn);
    fclose(file_nodes);
    char fname_freq[] = "frequencies.csv";
    FILE* file_freq = fopen(fname_freq, "r");
    if (file_freq == NULL) {
        printf("Cannot open file %s\n", fname_freq);
        return -13;
    }
    size_t ns = count_lines(file_freq);
    printf("num. frequencies: %li\n", ns);
    fclose(file_freq);
    file_freq = fopen(fname_freq, "r");
    double* freq = malloc(ns * sizeof(double));
    for (size_t i = 0; i < ns; i++) {
        fscanf(file_freq, "%lf", freq + i);
    }
    fclose(file_freq);
    run_case(sigma0, freq, ns, inj_node, ne, nn);
    free(freq);
    double end_time = omp_get_wtime();
    printf("Elapsed time: %.2f minutes\n", (end_time - start_time) / 60.0);
    return 0;
}
