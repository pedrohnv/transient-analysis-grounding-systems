/*
Test case 12pwrd01grcev-a

Reproducing the results in [1] for a grounding grid.

[1] L. D. Grcev and M. Heimbach, "Frequency dependent and transient
characteristics of substation grounding systems," in IEEE Transactions on
Power Delivery, vol. 12, no. 1, pp. 172-178, Jan. 1997.
doi: 10.1109/61.568238
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <cubature.h>
#include <Electrode.h>
#include <auxiliary.h>
//#include <omp.h>

/*
argv[0] = grid size
argv[1] = num_electrodes
argv[2] = num_nodes
*/
int main(int argc, char *argv[])
{
    int gs = strtol(argv[1], NULL, 10);;
    int num_electrodes = strtol(argv[2], NULL, 10);;
    int num_nodes = strtol(argv[3], NULL, 10);;
    char file_name[50];
    sprintf(file_name, "examples/ex12pwrd01grcev/gs%d.dat", gs);
    FILE* save_file = fopen(file_name, "w");
    // parameters
    double sigma = 1.0/1000.0; // soil conductivity
    double er = 10.0; // soil rel. permitivitty
    double rho_c = 1.9e-6; // copper resistivity
    // frequencies of interest
    int nf = 150;
    double freq[nf];
    logspace(2, 6.4, nf, freq);

    // electrode definition
    Electrode* electrodes = (Electrode*) malloc(sizeof(Electrode)*num_electrodes);
    sprintf(file_name, "examples/ex12pwrd01grcev/elec_gs%d.csv", gs);
    int read;
    read = electrodes_file(file_name, electrodes, num_electrodes);
    if (read > 0)
    {
        exit(read);
    }
    Electrode* images = (Electrode*) malloc(sizeof(Electrode)*num_electrodes);
    electrodes_file(file_name, images, num_electrodes);

    double nodes[num_nodes][3];
    sprintf(file_name, "examples/ex12pwrd01grcev/nodes_gs%d.csv", gs);
    nodes_file(file_name, nodes, num_nodes);
    for (int m = 0; m < num_electrodes; m++)
    {
        images[m].start_point[2] = -images[m].start_point[2];
        images[m].end_point[2] = -images[m].end_point[2];
        images[m].middle_point[2] = -images[m].middle_point[2];
    }
    int ne2 = num_electrodes*num_electrodes;
    int nn2 = num_nodes*num_nodes;
    int ss1 = (2*num_electrodes + num_nodes);
    int ss2 = ss1*ss1;
    _Complex double kappa, gamma, zinternal;
    _Complex double* zl = (_Complex double*) malloc(sizeof(_Complex double)*ne2);
    _Complex double* zt = (_Complex double*) malloc(sizeof(_Complex double)*ne2);
    _Complex double* yn = (_Complex double*) calloc(nn2, sizeof(_Complex double));
    _Complex double* ie = (_Complex double*) calloc(ss1, sizeof(_Complex double));
    _Complex double* ie_cp = (_Complex double*) calloc(ss1, sizeof(_Complex double));
    _Complex double* we = (_Complex double*) malloc(sizeof(_Complex double)*ss2);
    _Complex double* we_cp = (_Complex double*) malloc(sizeof(_Complex double)*ss2);
    int i, k;
    ie[ss1 - num_nodes] = 1.0;
    fill_incidence(we, electrodes, num_electrodes, nodes, num_nodes);
    // solve for each frequency: WE*VE = IE
    _Complex double s;
    _Complex double ref_l, ref_t;
    for (i = 0; i < nf; i++)
    {
        s = I*TWO_PI*freq[i];
        kappa = (sigma + s*er*EPS0); //soil complex conductivity
        gamma = csqrt(s*MU0*kappa); //soil propagation constant
        ref_t = (kappa - s*EPS0)/(kappa + s*EPS0);
        ref_l = (1 - ref_t);
        //TODO especialized impedance calculation taking advantage of symmetry
        calculate_impedances(
            electrodes, num_electrodes, zl, zt, gamma, s, MU0, kappa,
            200, 1e-3, 1e-4, ERROR_PAIRED, INTG_DOUBLE);
        zinternal = internal_impedance(s, rho_c,
            electrodes[0].radius, MU0)*electrodes[0].length;
        for (k = 0; k < num_electrodes; k++)
        {
            zl[k*num_electrodes + k] += zinternal;
        }
        impedances_images(electrodes, images, num_electrodes, zl, zt, gamma,
            s, MU0, kappa, ref_l, ref_t, 200, 1e-3, 1e-4, ERROR_PAIRED, INTG_DOUBLE);
        fill_impedance(we, electrodes, num_electrodes, num_nodes, zl, zt, yn);
        //The matrices are pivoted in-place. To recover them, copy
        copy_array(we, we_cp, ss2);
        copy_array(ie, ie_cp, ss1);
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
