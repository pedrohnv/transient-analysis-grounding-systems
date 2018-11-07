/*
Dissertação Miranda, caso 6.4.

Malha de aterramento em solo de duas camadas.
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
    FILE* save_file = fopen("examples/miranda.dat", "w");

    // parameters
    double h; // image distance
    double mesh_depth = 1.0; //burial depth
    double layer_depth = 10.0; //first layer size
    double radius = 1.25e-2; //radius FIXME
    double sigma1 = 0.001; // soil conductivity
    double sigma2 = 0.010; // soil conductivity
    double er1 = 30.0; //soil rel. permitivitty
    double er2 = 50.0; //soil rel. permitivitty
    double rho_c = 1.9e-6; //copper resistivity
    // frequencies of interest
    int nf = 150;
    double freq[nf];
    logspace(2, 7, nf, freq);
    // electrode definition and segmentation
    int num_electrodes = 118;
    Electrode* electrodes = (Electrode*) malloc(sizeof(Electrode)*num_electrodes);
    electrodes_file("examples/miranda64_electrodes.txt", electrodes, num_electrodes);
    int num_nodes = 116;
    double nodes[num_nodes][3];
    nodes_file("examples/miranda64_nodes.txt", nodes, num_nodes);
    //make images as copy of electrodes, change the points coordinates later
    Electrode* images = (Electrode*) malloc(sizeof(Electrode)*num_electrodes);
    electrodes_file("examples/miranda64_electrodes.txt", images, num_electrodes);

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
    _Complex double kappa1, kappa2, gamma, zinternal;
    double w;
    int i, k, m;
    for (i = 0; i < nf; i++)
    {
        printf("i = %i\n", i);
        w = TWO_PI*freq[i];
        kappa1 = (sigma1 + I*w*er1*EPS0); //soil complex conductivity
        kappa2 = (sigma2 + I*w*er2*EPS0);
        gamma = csqrt(I*w*MU0*kappa1); //soil 1 propagation constant
        calculate_impedances(
            electrodes, num_electrodes, zl, zt, gamma, w, MU0, kappa1,
            200, 1e-3, 1e-4, ERROR_PAIRED, INTG_DOUBLE);
        zinternal = internal_impedance(w, rho_c, radius, MU0)*electrodes[0].length;
        for (k = 0; k < num_electrodes; k++)
        {
            zl[k*num_electrodes + k] += zinternal;
        }
        // considering up to 28 reflections, which gives 28*2 image groups
        // they can be divided into 4 airthmetic sequences of images
        // hence, 14 of each
        for (k = 1; k < 15; k++)
        {
            // === AIR ===
            ref_t = (kappa1 - I*w*EPS0)/(kappa1 + I*w*EPS0);
            // first group, Air
            h = 2*(k - 1)*layer_depth + 2*mesh_depth;
            for (m = 0; m < num_electrodes; m++)
            {
                images[m].start_point[2] = h;
                images[m].end_point[2] = h;
            }
            impedances_images(electrodes, images, num_electrodes, zl, zt, gamma,
                w, MU0, kappa1, ref_l, ref_t, 200, 1e-3, 1e-4, ERROR_PAIRED, INTG_DOUBLE);
            // second group, Air
            h = 2*k*layer_depth;
            for (m = 0; m < num_electrodes; m++)
            {
                images[m].start_point[2] = h;
                images[m].end_point[2] = h;
            }
            impedances_images(electrodes, images, num_electrodes, zl, zt, gamma,
                w, MU0, kappa1, ref_l, ref_t, 200, 1e-3, 1e-4, ERROR_PAIRED, INTG_DOUBLE);

            // === 2nd SOIL ===
            ref_t = (kappa1 - kappa2)/(kappa1 + kappa2);
            // third group, Soil
            // this group has the same distance as the second one (in air)
            impedances_images(electrodes, images, num_electrodes, zl, zt, gamma,
                w, MU0, kappa1, ref_l, ref_t, 200, 1e-3, 1e-4, ERROR_PAIRED, INTG_DOUBLE);
            // fourth group, Soil
            h = 2*k*layer_depth - 2*mesh_depth;
            for (m = 0; m < num_electrodes; m++)
            {
                images[m].start_point[2] = h;
                images[m].end_point[2] = h;
            }
            impedances_images(electrodes, images, num_electrodes, zl, zt, gamma,
                w, MU0, kappa1, ref_l, ref_t, 200, 1e-3, 1e-4, ERROR_PAIRED, INTG_DOUBLE);
        }
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
