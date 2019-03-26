/*
Test case 51emc03grcev

Reproducing the results in [1] for an horizontal electrode buriend in ground.

[1] L. Grcev and S. Grceva. “On HF Circuit Models of Horizontal Grounding
Electrodes”. In: IEEE Transactions on Electromagnetic Compatibility 51.3
(2009), pp. 873–875. ISSN: 0018-9375. DOI: 10.1109/TEMC.2009.2023330.
*/
#include "auxiliary.h"
#include "electrode.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
//#include <omp.h>

int
run_case (double length, double rho, char file_name[])
{
    // parameters
    double h = 0.8; //burial depth
    double r = 7e-3; //radius
    double sigma = 1/rho; // soil conductivity
    double er = 10.0; //soil rel. permitivitty
    double rho_c = 1.9e-6; //copper resistivity
    //char file_name[] = "20pwrd02grcev_L<>rho<>.dat";
    // frequencies of interest
    size_t nf = 250;
    double freq[nf];
    double start_point[3] = {0., 0., -h};
    double end_point[3] = {length, 0., -h};
    remove(file_name);
    FILE *save_file = fopen(file_name, "w");
    if (save_file == NULL) {
        printf("Cannot open file %s\n",  file_name);
        exit(1);
    }
    logspace(2, 7, nf, freq);
    // electrode definition and segmentation
    double lambda = wave_length(freq[nf - 1], sigma, er*EPS0, 1.0); //smallest
    size_t num_electrodes = ceil( length/(lambda/6.0) ) ;
    size_t num_nodes = num_electrodes + 1;
    double nodes[num_nodes][3];
    Electrode *electrodes = (Electrode*) malloc(sizeof(Electrode)*num_electrodes);
    // the internal impedance is added "outside" later
    segment_electrode(
        electrodes, nodes, num_electrodes, start_point, end_point, r, 0.0);
    // create images
    start_point[2] = h;
    end_point[2] = h;
    double nodes_images[num_nodes][3];
    //Electrode images[num_electrodes];
    Electrode *images = (Electrode*) malloc(sizeof(Electrode)*num_electrodes);
    segment_electrode(
        images, nodes_images, num_electrodes, start_point, end_point, r, 0.0);
    //build system to be solved
    size_t ne2 = num_electrodes*num_electrodes;
    size_t nn2 = num_nodes*num_nodes;
    size_t ss1 = (2*num_electrodes + num_nodes);
    size_t ss2 = ss1*ss1;
    _Complex double s;
    _Complex double kappa, gamma, zinternal;
    _Complex double *zl = malloc(sizeof(_Complex double)*ne2);
    _Complex double *zt = malloc(sizeof(_Complex double)*ne2);
    _Complex double *yn = calloc(nn2, sizeof(_Complex double));
    _Complex double *ie = calloc(ss1, sizeof(_Complex double));
    _Complex double *ie_cp = malloc(sizeof(_Complex double)*ss1);
    _Complex double *we = malloc(sizeof(_Complex double)*ss2);
    _Complex double *we_cp = malloc(sizeof(_Complex double)*ss2);
    ie[ss1 - num_nodes] = 1.0;
    fill_incidence(we, electrodes, num_electrodes, nodes, num_nodes);
    // solve for each frequency: WE*VE = IE
    for (size_t i = 0; i < nf; i++) {
        s = I*TWO_PI*freq[i];
        kappa = (sigma + s*er*EPS0); //soil complex conductivity
        gamma = csqrt(s*MU0*kappa); //soil propagation constant
        //TODO especialized impedance calculation taking advantage of symmetry
        calculate_impedances(electrodes, num_electrodes, zl, zt, gamma, s, 1.0,
                             kappa, 200, 1e-3, 1e-4, ERROR_PAIRED, INTG_DOUBLE);
        zinternal = internal_impedance(s, rho_c, r, 1.0)*electrodes[0].length;
        for (size_t k = 0; k < num_electrodes; k++) {
            zl[k*num_electrodes + k] += zinternal;
        }
        impedances_images(electrodes, images, num_electrodes, zl, zt, gamma,
                          s, 1.0, kappa, 0.0, 1.0, 200, 1e-3, 1e-4, ERROR_PAIRED,
                          INTG_DOUBLE);
        fill_impedance(we, electrodes, num_electrodes, num_nodes, zl, zt, yn);
        //The matrices are pivoted in-place. To recover them, copy
        for (size_t i = 0; i < ss2; i++) {
            we_cp[i] = we[i];
        }
        for (size_t i = 0; i < ss1; i++) {
            ie_cp[i] = ie[i];
        }
        solve_electrodes(we_cp, ie_cp, num_electrodes, num_nodes);
        fprintf(save_file, "%f %f\n",
                creal(ie_cp[ss1 - num_nodes]), cimag(ie_cp[ss1 - num_nodes]));
    }
    fclose(save_file);
    free(electrodes);
    free(images);
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
main (int argc, char *argv[])
{
    printf("test case 51emc03grcev\n=== START ===\n");
    run_case(10.0, 10.0, "examples/51emc03grcev_L10rho10.dat");
    run_case(10.0, 100.0, "examples/51emc03grcev_L10rho100.dat");
    run_case(10.0, 1000.0, "examples/51emc03grcev_L10rho1000.dat");
    run_case(100.0, 10.0, "examples/51emc03grcev_L100rho10.dat");
    run_case(100.0, 100.0, "examples/51emc03grcev_L100rho100.dat");
    run_case(100.0, 1000.0, "examples/51emc03grcev_L100rho1000.dat");
    printf("==== END ====\n");
    return 0;
}
