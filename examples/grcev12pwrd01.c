/*
Reproducing the results in [1] for a square grounding grid.
This example is single threaded (no parallelism).

[1] L. D. Grcev and M. Heimbach, "Frequency dependent and transient
characteristics of substation grounding systems," in IEEE Transactions on
Power Delivery, vol. 12, no. 1, pp. 172-178, Jan. 1997.
doi: 10.1109/61.568238
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "auxiliary.h"
#include "electrode.h"
#include "linalg.h"
#include "cubature.h"
#include "grid.h"

/** Run example using admittance formulation
@param gs grid size in meters
@param div how many segments each edge of the grid has (that is, each segment
will have length $L = 10 / div$)
@param nf number of frequencies logspaced between $10^{2}$ and $10^{6.4}$ Hz
*/
int
run_case_adm (int gs, int div, size_t nf)
{
    // integration parameters
    size_t max_eval = 0;
    double req_abs_error = 1e-4;
    double req_rel_error = 1e-5;
    // soil parameters, constant with frequency
    double sigma = 1.0 / 1000.0;  // soil conductivity
    double er = 10.0;  // soil rel. permitivitty
    // frequencies of interest
    double freq[nf];
    logspace(2, 6.4, nf, freq);
    // electrode definition
    int err;
    Grid grid = {gs/10 + 1, gs/10 + 1, gs, gs, div, div, 7e-3, -0.5};
    size_t num_electrodes = number_segments(grid);
    size_t num_nodes = number_nodes(grid);
    printf("Num. segments = %li\nNum. nodes = %li\n", num_electrodes, num_nodes);
    Electrode* electrodes = malloc(sizeof(Electrode)*num_electrodes);
    Electrode* images = malloc(sizeof(Electrode)*num_electrodes);
    double* nodes = malloc(num_nodes * 3 * sizeof(double));
    electrode_grid(grid, electrodes, nodes);
    electrode_grid(grid, images, nodes);
    for (size_t m = 0; m < num_electrodes; m++) {
        images[m].start_point[2]  = -images[m].start_point[2];
        images[m].end_point[2]    = -images[m].end_point[2];
        images[m].middle_point[2] = -images[m].middle_point[2];
    }
    size_t ne = num_electrodes;
    size_t nn = num_nodes;
    size_t ne2 = ne * ne;
    size_t nn2 = nn * nn;
    _Complex double* potzl = malloc(ne2 * sizeof(_Complex double));
    _Complex double* potzt = malloc(ne2 * sizeof(_Complex double));
    _Complex double* potzli = malloc(ne2 * sizeof(_Complex double));
    _Complex double* potzti = malloc(ne2 * sizeof(_Complex double));
    _Complex double* a = malloc((ne * nn) * sizeof(_Complex double));
    _Complex double* b = malloc((ne * nn) * sizeof(_Complex double));
    // Incidence and "z-potential" (mHEM) matrices =============================
    err = fill_incidence_adm(a, b, electrodes, ne, nodes, nn);
    if (err != 0) printf("Could not build incidence matrices\n");
    err = calculate_impedances(potzl, potzt, electrodes, ne, 0.0, 0.0, 0.0, 0.0,
                               max_eval, req_abs_error, req_rel_error, INTG_MHEM);
    if (err != 0) printf("integration error\n");
    err = impedances_images(potzli, potzti, electrodes, images, ne,
                            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, max_eval,
                            req_abs_error, req_rel_error, INTG_MHEM);
    if (err != 0) printf("integration error\n");
    // pre-compute distance matrices
    double* rbar = malloc(ne2 * sizeof(double));
    double* rbari = malloc(ne2 * sizeof(double));
    for (size_t m = 0; m < ne; m++) {
        for (size_t k = m; k < ne; k++) {
            rbar[m * ne + k] = vector_length(electrodes[k].middle_point,
                                             electrodes[m].middle_point);
            rbari[m * ne + k] = vector_length(electrodes[k].middle_point,
                                              images[m].middle_point);
        }
    }
    _Complex double* zh = malloc(nf * sizeof(_Complex double));
    openblas_set_num_threads(1);  // this is important! We want to multithread the frequency loop, not the matrix operations
    #pragma omp parallel
    {
        _Complex double* restrict zl = malloc(ne2 * sizeof(_Complex double));
        _Complex double* restrict zt = malloc(ne2 * sizeof(_Complex double));
        _Complex double* restrict ie = malloc(nn2 * sizeof(_Complex double));
        _Complex double* restrict yn = malloc(nn2 * sizeof(_Complex double));
        _Complex double s, kappa, gamma, ref_l, ref_t, exp_gr, iwu_4pi, one_4pik;
        #pragma omp for
        for (size_t i = 0; i < nf; i++) {
            ie[0] = 1.0;
            for (size_t m = 1; m < num_nodes; m++) {
                ie[m] = 0.0;
            }
            s = I * TWO_PI * freq[i];
            kappa = (sigma + s * er * EPS0);  // soil complex conductivity
            gamma = csqrt(s * MU0 * kappa);  // soil propagation constant
            ref_t = (kappa - s * EPS0) / (kappa + s * EPS0);
            ref_l = 1.0;  // longitudinal current has positive image
            // modified HEM (mHEM):
            iwu_4pi = s * MU0 / (FOUR_PI);
            one_4pik = 1.0 / (FOUR_PI * kappa);
            for (size_t m = 0; m < ne; m++) {
                for (size_t k = m; k < ne; k++) {
                    exp_gr = cexp(-gamma * rbar[m * ne + k]);
                    zl[m * ne + k] = exp_gr * potzl[m * ne + k];
                    zt[m * ne + k] = exp_gr * potzt[m * ne + k];
                    exp_gr = cexp(-gamma * rbari[m * ne + k]);
                    zl[m * ne + k] += ref_l * exp_gr * potzli[m * ne + k];
                    zt[m * ne + k] += ref_t * exp_gr * potzti[m * ne + k];
                    zl[m * ne + k] *= iwu_4pi;
                    zt[m * ne + k] *= one_4pik;
                }
            }
            // Traditional HEM (highly discouraged):
            /*calculate_impedances(zl, zt, electrodes, ne, gamma, s, 1.0, kappa,
                                max_eval, req_abs_error, req_rel_error, INTG_DOUBLE);
            impedances_images(zl, zt, electrodes, images, ne, gamma, s, 1.0,
                            kappa, ref_l, ref_t, max_eval, req_abs_error,
                            req_rel_error, INTG_DOUBLE);*/
            fill_impedance_adm(yn, zl, zt, a, b, num_electrodes, num_nodes);
            err = solve_admittance(yn, ie, num_nodes);
            if (err != 0) exit(err);
            zh[i] = ie[0];
        }
        free(zl);
        free(zt);
        free(ie);
        free(yn);
    }  // end parallel
    free(electrodes);
    free(images);
    free(nodes);
    free(potzl);
    free(potzt);
    free(potzli);
    free(potzti);
    free(rbar);
    free(rbari);
    free(a);
    free(b);
    char file_name[50];
    sprintf(file_name, "gs%d.csv", gs);
    FILE *save_file = fopen(file_name, "w");
    for (int i = 0; i < nf; i++) {
        fprintf(save_file, "%f + %f%s\n", creal(zh[i]), cimag(zh[i]), "im");
    }
    free(zh);
    fclose(save_file);
    return 0;
}

int
run_case (int gs, int div, size_t nf)
{
    return run_case_adm(gs, div, nf);
}

/**
Run cases GS10, GS20, GS30, GS60 and GS120, all with segments of length
(10/Lmax) meters and nf frequencies.
*/
int
sweep ()
{
    printf("Making a sweep of many cases.\n");
    size_t nf = 100;
    double begin;  // to estimate time until completion
    double len = 1.0;  // Lmax
    int div = ceil(10.0 / len);
    // ===================================================
    printf("computing GS10\n");
    begin = omp_get_wtime();
    run_case(10, div, nf);
    printf("GS10: end; elapsed time: %f s\n", omp_get_wtime() - begin);
    // ===================================================
    printf("computing GS20\n");
    begin = omp_get_wtime();
    run_case(20, div, nf);
    printf("GS20: end; elapsed time: %f s\n", omp_get_wtime() - begin);
    // ===================================================
    printf("computing GS30\n");
    begin = omp_get_wtime();
    run_case(30, div, nf);
    printf("GS30: end; elapsed time: %f s\n", omp_get_wtime() - begin);
    // ===================================================
    printf("computing GS60\n");
    begin = omp_get_wtime();
    run_case(60, div, nf);
    printf("GS60: end; elapsed time: %f s\n", omp_get_wtime() - begin);
    // ===================================================
    printf("computing GS120\n");
    begin = omp_get_wtime();
    run_case(120, div, nf);
    printf("GS120: end; elapsed time: %f s\n", omp_get_wtime() - begin);
    // ===================================================
    return 0;
}

int
main (int argc, char *argv[])
{
    char *p;
    if (argc == 2 && strtol(argv[1], &p, 10) == 0) return sweep();
    if (argc < 3) {
        printf("wrong number of arguments, must be 3:\n");
        printf("  gs  : grid size in meters\n");
        printf("  len : segments' maximum length in meters\n");
        printf("  nf  : number of frequencies (logspace from 10^2 to 10^(6.4) Hz)\n");
        exit(argc);
    }
    int gs  = strtol(argv[1], &p, 10);
    float len = strtof(argv[2], &p);
    int nf  = strtol(argv[3], &p, 10);
    printf("gs  = %i x %i m^2\n", gs, gs);
    printf("len = %.2f m\n", len);
    printf("nf  = %i\n", nf);
    double time_spent;
    double begin = omp_get_wtime();
    int div = ceil(10/len);
    run_case(gs, div, nf);
    time_spent = (double) (omp_get_wtime() - begin);
    printf("elapsed time: %f s\n", time_spent);
    return 0;
}
