/*
Test case grcev51emc03

Reproducing the results in [1] for an horizontal electrode buriend in ground.

[1] L. Grcev and S. Grceva. “On HF Circuit Models of Horizontal Grounding
Electrodes”. In: IEEE Transactions on Electromagnetic Compatibility 51.3
(2009), pp. 873–875. ISSN: 0018-9375. DOI: 10.1109/TEMC.2009.2023330.
*/
#include "auxiliary.h"
#include "electrode.h"
#include "linalg.h"
#include "cubature.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
//#include <omp.h>

int
run_case_imm (double length, double rho, char file_name[])
{
    size_t max_eval = 0;
    double req_abs_error = 1e-3;
    double req_rel_error = 1e-4;
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
    segment_electrode(electrodes, nodes, num_electrodes, start_point, end_point,
                      r, 0.0);
    // create images
    start_point[2] = h;
    end_point[2] = h;
    double nodes_images[num_nodes][3];
    //Electrode images[num_electrodes];
    Electrode *images = (Electrode*) malloc(sizeof(Electrode)*num_electrodes);
    segment_electrode(images, nodes_images, num_electrodes, start_point,
                      end_point, r, 0.0);
    //build system to be solved
    size_t ne2 = num_electrodes*num_electrodes;
    size_t nn2 = num_nodes*num_nodes;
    size_t ss1 = (2*num_electrodes + num_nodes);
    size_t ss2 = ss1*ss1;
    _Complex double s, ref_l, ref_t;
    _Complex double kappa, gamma, zinternal, kappa_air;
    _Complex double *zl = malloc(ne2*sizeof(_Complex double));
    _Complex double *zt = malloc(ne2*sizeof(_Complex double));
    _Complex double *ye = calloc(nn2, sizeof(_Complex double));
    _Complex double *ie = calloc(ss1, sizeof(_Complex double));
    _Complex double *ie_cp = malloc(ss1*sizeof(_Complex double));
    _Complex double *we = malloc(ss2*sizeof(_Complex double));
    _Complex double *we_cp = malloc(ss2*sizeof(_Complex double));
    fill_incidence_imm(we, electrodes, num_electrodes, nodes, num_nodes);
    // solve for each frequency: WE*VE = IE
    ie[0] = 1.0;
    for (size_t i = 0; i < nf; i++) {
        s = I*TWO_PI*freq[i];
        kappa = (sigma + s*er*EPS0); //soil complex conductivity
        gamma = csqrt(s*MU0*kappa); //soil propagation constant
        kappa_air = (s*EPS0);
        ref_l = (kappa - kappa_air)/(kappa + kappa_air);
        ref_t = ref_l;
        //TODO especialized impedance calculation taking advantage of symmetry
        calculate_impedances(electrodes, num_electrodes, zl, zt, gamma, s, 1.0,
                             kappa, max_eval, req_abs_error, req_rel_error,
                             ERROR_PAIRED, INTG_DOUBLE);
        zinternal = internal_impedance(s, rho_c, r, 1.0)*electrodes[0].length;
        for (size_t k = 0; k < num_electrodes; k++) {
            zl[k*num_electrodes + k] += zinternal;
        }
        impedances_images(electrodes, images, num_electrodes, zl, zt, gamma,
                          s, 1.0, kappa, ref_l, ref_t, max_eval, req_abs_error,
                          req_rel_error, ERROR_PAIRED, INTG_DOUBLE);
        fill_impedance_imm(we, num_electrodes, num_nodes, zl, zt, ye);
        //The matrices are pivoted in-place. To recover them, copy
        for (size_t i = 0; i < ss2; i++) {
            we_cp[i] = we[i];
        }
        for (size_t i = 0; i < ss1; i++) {
            ie_cp[i] = ie[i];
        }
        solve_immittance(we_cp, ie_cp, num_electrodes, num_nodes);
        fprintf(save_file, "%f %f\n", creal(ie_cp[0]), cimag(ie_cp[0]));
    }
    fclose(save_file);
    free(electrodes);
    free(images);
    free(zl);
    free(zt);
    free(ye);
    free(ie);
    free(ie_cp);
    free(we);
    free(we_cp);
    return 0;
}

int
main (int argc, char *argv[])
{
    printf("test case grcev51emc03\n=== START ===\n");
    run_case_imm(10.0, 10.0, "grcev51emc03_L10rho10.dat");
    run_case_imm(10.0, 100.0, "grcev51emc03_L10rho100.dat");
    run_case_imm(10.0, 1000.0, "grcev51emc03_L10rho1000.dat");
    run_case_imm(100.0, 10.0, "grcev51emc03_L100rho10.dat");
    run_case_imm(100.0, 100.0, "grcev51emc03_L100rho100.dat");
    run_case_imm(100.0, 1000.0, "grcev51emc03_L100rho1000.dat");
    printf("==== END ====\n");
    return 0;
}
