/*
400m^2 grounding mesh.
*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cubature.h>
#include <Electrode.h>
#include <auxiliary.h>
//#include <omp.h>

int main()
{
    /*FILE* zl_file = fopen("examples/zl.dat", "w");
    FILE* zt_file = fopen("examples/zt.dat", "w");
    FILE* yl_file = fopen("examples/yl.dat", "w");
    FILE* yt_file = fopen("examples/yt.dat", "w");
    FILE* a_file = fopen("examples/a.dat", "w");
    FILE* b_file = fopen("examples/b.dat", "w");
    FILE* yn_file = fopen("examples/yn.dat", "w");*/
    FILE* save_file = fopen("examples/malha01.dat", "w");
    // parameters
    double h = 1.0; // burial depth
    double radius = 0.7e-3;//12.5e-3; // radius
    double sigma1 = 1e-3/5.0;//0.001; // soil conductivity
    double er1 = 30.0; //soil rel. permitivitty
    double rho_c = 1.9e-6; // copper resistivity
    // frequencies of interest
    int nf = 150;
    double freq[nf];
    logspace(2, 7, nf, freq);
    // electrode definition and segmentation
    int num_electrodes = 200;
    int num_nodes = 185;
    Electrode* electrodes = (Electrode*) malloc(sizeof(Electrode)*num_electrodes);
    electrodes_file("examples/malha01_electrodes.txt", electrodes, num_electrodes);
    double nodes[num_nodes][3];
    nodes_file("examples/malha01_nodes.txt", nodes, num_nodes);
    Electrode* images = (Electrode*) malloc(sizeof(Electrode)*num_electrodes);
    electrodes_file("examples/malha01_electrodes.txt", images, num_electrodes);
    for (int m = 0; m < num_electrodes; m++)
    {
        images[m].start_point[2] = -images[m].start_point[2];
        images[m].end_point[2] = -images[m].end_point[2];
        images[m].middle_point[2] = -images[m].middle_point[2];
    }
    int ne2 = num_electrodes*num_electrodes;
    int nn2 = num_nodes*num_nodes;
    _Complex double* zl = (_Complex double*) malloc(sizeof(_Complex double)*ne2);
    _Complex double* zt = (_Complex double*) malloc(sizeof(_Complex double)*ne2);
    _Complex double* yn = (_Complex double*) malloc(sizeof(_Complex double)*nn2);
    _Complex double ref_l = 0.0; //reflection coefficient, longitudinal
    _Complex double ref_t; //reflection coefficient, transversal
    double* a = (double*) malloc(sizeof(double)*(num_electrodes*num_nodes));
    double* b = (double*) malloc(sizeof(double)*(num_electrodes*num_nodes));
    incidence_alt(a, b, electrodes, num_electrodes, nodes, num_nodes);
    // for each frequency
    _Complex double kappa1, gamma, zinternal;
    _Complex double s;
    int i, k;
    for (i = 0; i < nf; i++)
    {
        printf("i = %i\n", i);
        s = I*TWO_PI*freq[i];
        kappa1 = (sigma1 + s*er1*EPS0); //soil complex conductivity
        gamma = csqrt(s*MU0*kappa1); //soil 1 propagation constant
        calculate_impedances(
            electrodes, num_electrodes, zl, zt, gamma, w, MU0, kappa1,
            200, 1e-3, 1e-4, ERROR_PAIRED, INTG_DOUBLE);
        zinternal = internal_impedance(w, rho_c, radius, MU0);
        for (k = 0; k < num_electrodes; k++)
        {
            zl[k*num_electrodes + k] += zinternal*electrodes[k].length;
        }
        ref_t = (kappa1 - s*EPS0)/(kappa1 + s*EPS0);
        impedances_images(electrodes, images, num_electrodes, zl, zt, gamma,
            s, MU0, kappa1, ref_l, ref_t, 200, 1e-3, 1e-4, ERROR_PAIRED, INTG_DOUBLE);
        /*ynodal_eq(yn, a, b, zl, zt, num_electrodes, num_nodes);
        print_zmatrix_file(num_nodes, num_nodes, yn, num_nodes, save_file);
        fprintf(save_file, "f = %f\n", freq[i]);*/

        /*print_dmatrix_file(num_electrodes, num_nodes, a, num_nodes, a_file);
        print_dmatrix_file(num_electrodes, num_nodes, b, num_nodes, b_file);
        print_zmatrix_file(num_electrodes, num_electrodes, zl, num_electrodes, zl_file);
        print_zmatrix_file(num_electrodes, num_electrodes, zt, num_electrodes, zt_file);*/
        ynodal_eq(yn, a, b, zl, zt, num_electrodes, num_nodes);
        print_zmatrix_file(num_nodes, num_nodes, yn, num_nodes, save_file);
        fprintf(save_file, "f = %f\n", freq[i]);
        /*print_zmatrix_file(num_nodes, num_nodes, yn, num_nodes, yn_file);
        print_zmatrix_file(num_electrodes, num_electrodes, zl, num_electrodes, yl_file);
        print_zmatrix_file(num_electrodes, num_electrodes, zt, num_electrodes, yt_file);*/
    }
    free(electrodes);
    free(images);
    free(zl);
    free(zt);
    free(yn);
    free(a);
    free(b);
    fclose(save_file);
    /*fclose(zl_file);
    fclose(zt_file);
    fclose(a_file);
    fclose(b_file);*/
    return 0;
}
