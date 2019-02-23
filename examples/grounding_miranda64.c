/*
Dissertação Miranda, caso 6.4.

Malha de aterramento em solo de duas camadas.
*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <cubature.h>
#include <Electrode.h>
#include <auxiliary.h>
//#include <omp.h>
#include <mkl_lapacke.h>

int main()
{
    clock_t begin, end;
    double time_spent;
    begin = clock();
    char file_name[] = "miranda.dat";
    FILE* save_file = fopen(file_name, "w");
    if (save_file == NULL)
    {
        printf("Cannot open file %s\n",  file_name);
        exit(1);
    }
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
    Electrode* electrodes = malloc(sizeof(Electrode)*num_electrodes);
    electrodes_file("miranda64_electrodes.txt", electrodes, num_electrodes);
    int num_nodes = 116;
    double nodes[num_nodes][3];
    nodes_file("miranda64_nodes.txt", nodes, num_nodes);
    //make images as copy of electrodes, change the points coordinates later
    Electrode* images = malloc(sizeof(Electrode)*num_electrodes);
    electrodes_file("miranda64_electrodes.txt", images, num_electrodes);

    int ne2 = num_electrodes*num_electrodes;
    int nn2 = num_nodes*num_nodes;
    _Complex double* zl = malloc(sizeof(_Complex double)*ne2);
    _Complex double* zt = malloc(sizeof(_Complex double)*ne2);
    _Complex double* yn = malloc(sizeof(_Complex double)*nn2);
    _Complex double* ie = malloc(sizeof(_Complex double)*num_nodes);
    double* a = malloc(sizeof(double)*(num_electrodes*num_nodes));
    double* b = malloc(sizeof(double)*(num_electrodes*num_nodes));
    incidence_alt(a, b, electrodes, num_electrodes, nodes, num_nodes);
    // for each frequency
    _Complex double kappa1, kappa2, gamma, zinternal, s;
    _Complex double ref_l = 0.0; //reflection coefficient, longitudinal
    _Complex double ref_t; //reflection coefficient, transversal
    int i, k, m, info;
    MKL_INT n = num_electrodes*2 + num_nodes;
    MKL_INT ipiv[n]; //pivot indices
    for (i = 0; i < nf; i++)
    {
        //reset IN
        ie[0] = 1.0;
        for (size_t k = 1; k < num_nodes; k++)
        {
            ie[k] = 0.0;
        }
        printf("i = %i\n", i);
        s = I*TWO_PI*freq[i];
        kappa1 = (sigma1 + s*er1*EPS0); //soil complex conductivity
        kappa2 = (sigma2 + s*er2*EPS0);
        gamma = csqrt(s*MU0*kappa1); //soil 1 propagation constant
        calculate_impedances(
            electrodes, num_electrodes, zl, zt, gamma, s, 1.0, kappa1,
            200, 1e-3, 1e-4, ERROR_PAIRED, INTG_DOUBLE);
        zinternal = internal_impedance(s, rho_c, radius, 1.0)*electrodes[0].length;
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
            ref_t = (kappa1 - s*EPS0)/(kappa1 + s*EPS0);
            // first group, Air
            h = 2*(k - 1)*layer_depth + 2*mesh_depth;
            for (m = 0; m < num_electrodes; m++)
            {
                images[m].start_point[2] = h;
                images[m].end_point[2] = h;
            }
            impedances_images(electrodes, images, num_electrodes, zl, zt, gamma,
                s, 1.0, kappa1, ref_l, ref_t, 200, 1e-3, 1e-4, ERROR_PAIRED, INTG_DOUBLE);
            // second group, Air
            h = 2*k*layer_depth;
            for (m = 0; m < num_electrodes; m++)
            {
                images[m].start_point[2] = h;
                images[m].end_point[2] = h;
            }
            impedances_images(electrodes, images, num_electrodes, zl, zt, gamma,
                s, 1.0, kappa1, ref_l, ref_t, 200, 1e-3, 1e-4, ERROR_PAIRED, INTG_DOUBLE);

            // === 2nd SOIL ===
            ref_t = (kappa1 - kappa2)/(kappa1 + kappa2);
            // third group, Soil
            // this group has the same distance as the second one (in air)
            //TODO avoid this calculation, as it is done twice unnecessarily
            impedances_images(electrodes, images, num_electrodes, zl, zt, gamma,
                s, 1.0, kappa1, ref_l, ref_t, 200, 1e-3, 1e-4, ERROR_PAIRED, INTG_DOUBLE);
            // fourth group, Soil
            h = 2*k*layer_depth - 2*mesh_depth;
            for (m = 0; m < num_electrodes; m++)
            {
                images[m].start_point[2] = h;
                images[m].end_point[2] = h;
            }
            impedances_images(electrodes, images, num_electrodes, zl, zt, gamma,
                s, 1.0, kappa1, ref_l, ref_t, 200, 1e-3, 1e-4, ERROR_PAIRED, INTG_DOUBLE);
        }
        ynodal_eq(yn, a, b, zl, zt, num_electrodes, num_nodes);
        info = LAPACKE_zgesv(LAPACK_ROW_MAJOR, num_nodes, 1, yn, num_nodes, ipiv, ie, 1);
        // Check for the exact singularity
        if(info > 0)
        {
            printf("The diagonal element of the triangular factor of YN,\n");
            printf("U(%i,%i) is zero, so that YN is singular;\n", info, info);
            printf("the solution could not be computed.\n");
            exit(info);
        }
        fprintf(save_file, "%f %f\n", creal(ie[0]), cimag(ie[0]));
    }
    free(electrodes);
    free(images);
    free(zl);
    free(zt);
    free(yn);
    free(a);
    free(b);
    free(ie);
    fclose(save_file);
    end = clock();
    time_spent = (double) (end - begin)/CLOCKS_PER_SEC;
    printf("elapsed time: %f s\n", time_spent);
    return 0;
}
