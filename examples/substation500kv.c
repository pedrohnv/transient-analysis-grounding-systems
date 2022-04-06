/*
A grid buried at a depth of 0.5 m in a two layer soil (second layer at a depth
of 5.0 m) is simulated. It is calculated the harmonic impedance for every node
of the grid, defined in an external file.

The interfaces are accounted for using the image method. The electromagnetic
wave reflects at each interface, creating a new image with each successive
reflection. The effect of each image is added to the impedance matrices until
the smallest relative change (to an element in the matrices) is less than 0.01.

Soil parameters:
- first:
    ρ1 = 2100 Ω⋅m
    ϵr1 = 15
- second:
    ρ2 = 900 Ω⋅m
    ϵr2 = 15
    H = 5.0 m

To use less memory, the matrices of the system are float instead of float.
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
//#include <omp.h>
#include "auxiliary.h"
#include "electrode.h"
#include "linalg.h"
#include "cubature.h"
#include "grid.h"

int set_image_dist (Electrode *images, size_t num_electrodes, float h)
{
    for (size_t m = 0; m < num_electrodes; m++) {
        images[m].start_point[2] = h;
        images[m].end_point[2] = h;
        images[m].middle_point[2] = h;
    }
    return 0;
}

int
run_case ()
{
    char elec_file_name[] = "examples/substation500kv.txt";
    char zh_file_name[] = "substation500kv_zh.csv";
    //char gpr_file_name[] = "substation500kv_gpr.csv";
    //char stepv_file_name[] = "substation500kv_stepv.csv";
    // frequencies of interest
    size_t nf = 150;
    float freq[nf];
    logspace(2, 6.4, nf, freq);
    // Integration parameters
    size_t max_eval = 0;
    double req_abs_error = 1e-6;
    double req_rel_error = 1e-6;
    int intg_type = INTG_DOUBLE;
    // soil parameters
    double mur = 1.0;
    double sigma1 = 1.0 / 2100.0;  // soil conductivity
    double er1 = 15.0;  // soil rel. permitivitty
    double sigma2 = 1.0 / 900.0;  // soil conductivity
    double er2 = 15.0;  // soil rel. permitivitty
    float H = 5.0;  // second layer depth
    float h = 0.5;  // grid burial depth
    int c1, c2, k;
    // read electrodes from file (already segmented)
    int read;
    size_t num_electrodes = lines_in_file(elec_file_name);
    Electrode *electrodes = malloc(num_electrodes * sizeof(Electrode));
    read = electrodes_file(elec_file_name, electrodes, num_electrodes);
    if (read == -10) exit(read);
    Electrode *images = malloc(num_electrodes * sizeof(Electrode));
    for (size_t m = 0; m < num_electrodes; m++) {
        populate_electrode(&(images[m]), electrodes[m].start_point,
                           electrodes[m].end_point, electrodes[m].radius);
    }
    // malloc the maximum possible number of nodes:
    float *nodes = malloc(2 * num_electrodes * 3 * sizeof(float));
    size_t num_nodes = nodes_from_elecs(nodes, electrodes, num_electrodes);
    nodes = realloc(nodes, num_nodes * 3 * sizeof(float));
    size_t ne2 = num_electrodes * num_electrodes;
    size_t nn2 = num_nodes * num_nodes;
    // malloc matrices
    size_t nbytes = sizeof(_Complex float) * (ne2 * 4 + nn2 * 2 + (num_electrodes * num_nodes) * 2);
    printf("Expected memory needed for matrices: %f GB\n", nbytes * 1e-9);
    printf("Enter 0 if you want to continue: ");
    scanf("%i", &read);
    if (read != 0) exit(0);
    _Complex float ref_l10, ref_t10, ref_l12, ref_t12;
    _Complex float* zl = malloc(ne2 * sizeof(_Complex float));
    _Complex float* zt = malloc(ne2 * sizeof(_Complex float));
    _Complex float* zli = malloc(ne2 * sizeof(_Complex float));
    _Complex float* zti = malloc(ne2 * sizeof(_Complex float));
    _Complex float* ie = malloc(nn2 * sizeof(_Complex float));
    _Complex float* yn = malloc(nn2 * sizeof(_Complex float));
    _Complex float* a = malloc((num_electrodes * num_nodes) * sizeof(_Complex float));
    _Complex float* b = malloc((num_electrodes * num_nodes) * sizeof(_Complex float));
    fill_incidence_adm(a, b, electrodes, num_electrodes, nodes, num_nodes);
    FILE *zh_file = fopen(zh_file_name, "w");
    _Complex double kappa1, kappa2, gamma, s;
    int converged;
    printf("entering loop\n");
    fflush(stdout);
    // use raw BLAS and LAPACK to solve system with multiple RHS
    int nn1 = (int) num_nodes;
    int ldie = nn1;  // leading dimension of ie
    int* ipiv = malloc(nn1 * sizeof(int));  // pivot indices
    int info;
    int nrhs = num_nodes;
    char uplo = 'L';  // matrices are symmetric, only Lower Half is set
    int lwork = -1;  // signal to Query the optimal workspace
    _Complex float wkopt;
    csysv_(&uplo, &nn1, &nrhs, yn, &nn1, ipiv, ie, &ldie, &wkopt, &lwork, &info);
    lwork = crealf(wkopt);
    _Complex float* work = malloc(lwork * sizeof(_Complex float));
    for (size_t i = 0; i < nf; i++) {
        s = I * TWO_PI * freq[i];
        kappa1 = (sigma1 + s * er1 * EPS0);  // soil complex conductivity
        kappa2 = (sigma2 + s * er2 * EPS0);  // soil complex conductivity
        gamma = csqrt(s * MU0 * kappa1);  //soil propagation constant
        calculate_impedances(zl, zt, electrodes, num_electrodes, gamma, s, mur,
                             kappa1, max_eval, req_abs_error, req_rel_error, intg_type);
        // Images =====================
        // reflection coefficient, air
        ref_t10 = (kappa1 - s * EPS0) / (kappa1 + s * EPS0);
        ref_l10 = 1.0;
        // reflection coefficient, 2nd soil
        ref_t12 = (kappa1 - kappa2)/(kappa1 + kappa2);
        ref_l12 = 1.0;
        k = 0;
        converged = 0;
        while (!converged) {
            k += 1;
            for (size_t m = 0; m < ne2; m++) {
                zli[m] = 0.0;
                zti[m] = 0.0;
            }
            // second image group in Air
            set_image_dist(images, num_electrodes, 2*k*H);
            impedances_images(zli, zti, electrodes, images, num_electrodes, gamma,
                              s, mur, kappa1, 1, 1, max_eval, req_abs_error,
                              req_rel_error, intg_type);
            // first image group in 2nd Soil
            // has the same distance as second group in air
            for (size_t m = 0; m < ne2; m++) {
                zli[m] *= (ref_l10 + ref_l12);
                zti[m] *= (ref_t10 + ref_t12);
            }
            // first image group in Air
            set_image_dist(images, num_electrodes, 2*(k - 1)*H + 2*h);
            impedances_images(zli, zti, electrodes, images, num_electrodes, gamma,
                              s, mur, kappa1, ref_l10, ref_t10, max_eval,
                              req_abs_error, req_rel_error, intg_type);
            // second image group in 2nd Soil
            set_image_dist(images, num_electrodes, 2*k*H - 2*h);
            impedances_images(zli, zti, electrodes, images, num_electrodes, gamma,
                              s, mur, kappa1, ref_l12, ref_t12, max_eval,
                              req_abs_error, req_rel_error, intg_type);
            c1 = (cabsf(zli[0]) < req_rel_error * cabsf(zl[0]));
            c2 = (cabsf(zti[0]) < req_rel_error * cabsf(zt[0]));
            for (size_t m = 0; m < ne2; m++) {
                zl[m] += zli[m];
                zt[m] += zti[m];
            }
            converged = (c1 && c2);
        }
        printf("f = %.1f Hz, number of images: %i\n", freq[i], 2*k);
        fill_impedance_adm(yn, zl, zt, a, b, num_electrodes, num_nodes);
        for (size_t m = 0; m < num_nodes; m++) {
            for (size_t k = 0; m < num_nodes; m++) {
                if (m == k) {
                    ie[m * num_nodes + k] = 1.0;
                } else {
                    ie[m * num_nodes + k] = 0.0;
                }
            }
        }
        // solve
        csysv_(&uplo, &nn1, &nrhs, yn, &nn1, ipiv, ie, &ldie, work, &lwork, &info);
        // Check for the exact singularity
        if (info > 0) {
            printf("The diagonal element of the triangular factor of YN,\n");
            printf("U(%i,%i) is zero, so that YN is singular;\n", info, info);
            printf("the solution could not be computed.\n");
        }

        fprintf(zh_file, "%f + %f%s\n", crealf(ie[0]), cimagf(ie[0]), "im");
        fflush(zh_file);
    }
    fclose(zh_file);
    free(electrodes);
    free(images);
    free(nodes);
    free(zl);
    free(zt);
    free(zli);
    free(zti);
    free(ie);
    free(yn);
    free(a);
    free(b);
    return 0;
}

int
main ()
{
    clock_t begin, end;
    float time_spent;
    begin = clock();
    run_case();
    end = clock();
    time_spent = (float) (end - begin) / CLOCKS_PER_SEC;
    printf("elapsed time: %.2f hours\n", time_spent / 3600.0);
    return 0;
}
