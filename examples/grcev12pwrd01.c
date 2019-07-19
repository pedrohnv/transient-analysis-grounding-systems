/*
Test case grcev12pwrd01-a

Reproducing the results in [1] for a grounding grid.
Execution time is expected to be up to 20 min.

[1] L. D. Grcev and M. Heimbach, "Frequency dependent and transient
characteristics of substation grounding systems," in IEEE Transactions on
Power Delivery, vol. 12, no. 1, pp. 172-178, Jan. 1997.
doi: 10.1109/61.568238
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

int
run_case (int gs, int num_electrodes, int num_nodes)
{
    size_t max_eval = 0;
    double req_abs_error = 1e-3;
    double req_rel_error = 1e-4;
    char file_name[50];
    sprintf(file_name, "gs%d.dat", gs);
    FILE *save_file = fopen(file_name, "w");
    // parameters
    double sigma = 1.0/1000.0; // soil conductivity
    double er = 10.0; // soil rel. permitivitty
    double rho_c = 1.9e-6; // copper resistivity
    // frequencies of interest
    size_t nf = 150;
    double freq[nf];
    logspace(2, 6.4, nf, freq);
    // electrode definition
    int read;
    sprintf(file_name, "examples/grcev12pwrd01_auxfiles/elec_gs%d.txt", gs);
    Electrode *electrodes = malloc(sizeof(Electrode)*num_electrodes);
    read = electrodes_file(file_name, electrodes, num_electrodes);
    if (read == -10) exit(read);
    Electrode *images = malloc(sizeof(Electrode)*num_electrodes);
    electrodes_file(file_name, images, num_electrodes);
    //double nodes[num_nodes][3];
    sprintf(file_name, "examples/grcev12pwrd01_auxfiles/nodes_gs%d.txt", gs);
    double *nodes = malloc(num_nodes*3*sizeof(double));
    read = nodes_file(file_name, nodes, num_nodes);
    if (read == -10) exit(read);
    for (size_t m = 0; m < num_electrodes; m++) {
        images[m].start_point[2] = -images[m].start_point[2];
        images[m].end_point[2] = -images[m].end_point[2];
        images[m].middle_point[2] = -images[m].middle_point[2];
    }
    size_t ne2 = num_electrodes*num_electrodes;
    size_t nn2 = num_nodes*num_nodes;
    size_t ss1 = (2*num_electrodes + num_nodes);
    size_t ss2 = ss1*ss1;
    _Complex double kappa, gamma, zinternal;
    _Complex double *zl = malloc(sizeof(_Complex double)*ne2);
    _Complex double *zt = malloc(sizeof(_Complex double)*ne2);
    _Complex double *yn = calloc(nn2, sizeof(_Complex double));
    _Complex double *ie = calloc(ss1, sizeof(_Complex double));
    _Complex double *ie_cp = malloc(sizeof(_Complex double)*ss1);
    _Complex double *we = malloc(sizeof(_Complex double)*ss2);
    _Complex double *we_cp = malloc(sizeof(_Complex double)*ss2);
    ie[0] = 1.0;
    fill_incidence_imm(we, electrodes, num_electrodes, nodes, num_nodes);
    // solve for each frequency: WE*VE = IE
    _Complex double s;
    _Complex double ref_l, ref_t;
    for (size_t k = 0; k < nf; k++) {
        s = I*TWO_PI*freq[k];
        kappa = (sigma + s*er*EPS0); //soil complex conductivity
        gamma = csqrt(s*MU0*kappa); //soil propagation constant
        ref_t = (kappa - s*EPS0)/(kappa + s*EPS0);
        //ref_l = (1 - ref_t);
        ref_l = ref_t;
        //TODO especialized impedance calculation taking advantage of symmetry
        calculate_impedances(electrodes, num_electrodes, zl, zt, gamma, s, 1.0,
                             kappa, max_eval, req_abs_error, req_rel_error,
                             ERROR_PAIRED, INTG_DOUBLE);
        zinternal = internal_impedance(s, rho_c, electrodes[0].radius, 1.0, &read)
                    *electrodes[0].length;
        for (size_t m = 0; m < num_electrodes; m++) {
            zl[m*num_electrodes + m] += zinternal;
        }
        impedances_images(electrodes, images, num_electrodes, zl, zt, gamma,
                          s, 1.0, kappa, ref_l, ref_t, max_eval, req_abs_error,
                          req_rel_error, ERROR_PAIRED, INTG_DOUBLE);
        fill_impedance_imm(we, num_electrodes, num_nodes, zl, zt, yn);
        //The matrices are pivoted in-place. To recover them, copy
        for (size_t i = 0; i < ss2; i++) {
            we_cp[i] = we[i];
        }
        for (size_t i = 0; i < ss1; i++) {
            ie_cp[i] = ie[i];
        }
        read = solve_immittance(we_cp, ie_cp, num_electrodes, num_nodes);
        if (read != 0) exit(read);
        fprintf(save_file, "%f %f\n", creal(ie[0]), cimag(ie[0]));
    }
    fclose(save_file);
    free(electrodes);
    free(images);
    free(nodes);
    free(zl);
    free(zt);
    free(yn);
    free(ie);
    free(ie_cp);
    free(we);
    free(we_cp);
    return 0;
}

int
main ()
{
    clock_t begin, end;
    double time_spent;
    // ===================================================
    printf("computing GS10\n");
    begin = clock();
    run_case(10, 12, 12);
    end = clock();
    time_spent = (double) (end - begin)/CLOCKS_PER_SEC;
    printf("GS10: end; elapsed time: %f s\n", time_spent);
    // ===================================================
    printf("computing GS20\n");
    begin = clock();
    run_case(20, 36, 33);
    end = clock();
    time_spent = (double) (end - begin)/CLOCKS_PER_SEC;
    printf("GS20: end; elapsed time: %f s\n", time_spent);
    // ===================================================
    printf("computing GS30\n");
    begin = clock();
    run_case(30, 72, 64);
    end = clock();
    time_spent = (double) (end - begin)/CLOCKS_PER_SEC;
    printf("GS30: end; elapsed time: %f s\n", time_spent);
    // ===================================================
    printf("computing GS60\n");
    begin = clock();
    run_case(60, 252, 217);
    end = clock();
    time_spent = (double) (end - begin)/CLOCKS_PER_SEC;
    printf("GS60: end; elapsed time: %f s\n", time_spent);
    // ===================================================
    printf("computing GS120\n");
    begin = clock();
    run_case(120, 936, 793);
    end = clock();
    time_spent = (double) (end - begin)/CLOCKS_PER_SEC;
    printf("GS120: end; elapsed time: %f s\n", time_spent);
    // ===================================================
    return 0;
}
