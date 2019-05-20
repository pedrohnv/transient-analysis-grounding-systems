/** HP_HEM
High Performance implementation of the Hybrid Electromagnetic Model

Main program
*/
#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "linalg.h"
#include "electrode.h"
#include "auxiliary.h"

int
soil_param (double *sigma, double *epsr, double *mur)
{
    printf("Enter the soil parameters:\n");
    printf("rho [Ohm.m]: ");
    scanf("%lf", sigma);
    printf("relative permittivity, epsr: ");
    scanf("%lf", epsr);
    printf("relative permeability, mur: ");
    scanf("%lf", mur);
    printf("Values read:\n");
    printf("rho = %lf\n", *sigma);
    printf("epsr = %lf\n", *epsr);
    printf("mur = %lf\n", *mur);
    *sigma = 1.0/(*sigma);
    return 0;
}

int
freq_param (int nf, double *freq)
{
    int v1 = 0;
    while(v1 != 1 && v1 != 2) {
        printf("Do you want them to be linearly (1) or logarithmic (2) spaced? ");
        scanf("%i", &v1);
    }
    float f1, f2;
    if (v1 == 1) {
        printf("Enter minimum frequency [Hz]: ");
        scanf("%f", &f1);
        printf("Enter maximum frequency [Hz]: ");
        scanf("%f", &f2);
        linspace(f1, f2, nf, freq);
    } else {
        printf("Enter minimum frequency expoent: ");
        scanf("%f", &f1);
        printf("Enter maximum frequency expoent: ");
        scanf("%f", &f2);
        logspace(f1, f2, nf, freq);
    }
    printf("Frequencies [Hz]:\n");
    for (int i = 0; i < nf; i++) printf("  %lf\n", freq[i]);
    return 0;
}

int
grid_param (double *radius, double *depth, double *frac, int *num_electrodes,
            int *num_nodes, int *ve1, double *length1, int *le1, int *ve2,
            double *length2, int *le2, double maxfreq, double sigma, double epsr,
            double mur)
{
    printf("Enter the conductors' radius: ");
    scanf("%lf", radius);
    printf("radius = %lf\n", *radius);
    printf("Enter the grid burial depth: ");
    scanf("%lf", depth);
    printf("depth = %lf\n", *radius);
    double lambda = wave_length(maxfreq, sigma, epsr*EPS0, mur);
    double len1, len2, lmax;
    int l1, l2;
    char r = 'n';
    int v1, v2;
    while (r == 'N' || r == 'n') {
        printf("Segmentation will be done such that each segment has length\n");
        printf("at most lambda/n, where lambda is the wave length for the\n");
        printf("maximum frequency of interest.\n");
        printf("lambda [m] = %lf\n", lambda);
        printf("Enter `n` (6 is recommended): ");
        scanf("%lf", frac);
        printf("Enter the number of vertice rows `v1`: ");
        scanf("%i", &v1);
        printf("Enter the number of vertice columns `v2`: ");
        scanf("%i", &v2);
        printf("Enter the total length `L1`: ");
        scanf("%lf", &len1);
        printf("Enter the total length `L2`: ");
        scanf("%lf", &len2);
        lmax = (lambda/(*frac));
        printf("Lmax = %lf\n", lmax);
        l1 = ceil((len1/(v1 - 1))/lmax);
        l2 = ceil((len2/(v2 - 1))/lmax);
        *num_electrodes = l1*v2*(v1 - 1) + l2*v1*(v2 - 1);
        *num_nodes = v1*v2 + v1*(v2 - 1)*(l2 - 1) + v2*(v1 - 1)*(l1 - 1);
        printf("\nThis is a (%.2f x %.2f) grid that is going to have, for the\n",
               len1, len2);
        printf("chosen frequencies and segmentation, %i segments\n", *num_electrodes);
        printf("and %i nodes. Is that correct?\n", *num_nodes);
        printf("[y/n] ");
        r = getchar();
        if (r == '\n') r = getchar();
        while(r != 'n' && r != 'N' && r != 'y' && r != 'Y') {
            printf("invalid input, enter the choice(y/Y/n/N) again : ");
            r = getchar();
            if (r == '\n') r = getchar();
        }
    }
    *ve1 = v1;
    *ve2 = v2;
    *le1 = l1;
    *le2 = l2;
    *length1 = len1;
    *length2 = len2;
    return 0;
}

int
intg_param (size_t *max_eval, double *req_abs_error, double *req_rel_error)
{
    printf("Enter the maximum number of evaluations to be done\n");
    printf("during integration (0 for no limit): ");
    int me;
    scanf("%i", &me);
    *max_eval = (size_t) me;
    printf("Enter the requested absolute error for integration\n");
    printf("(0 to ignore): ");
    scanf("%lf", req_abs_error);
    printf("Enter the requested relative error for integration\n");
    printf("(0 to ignore): ");
    scanf("%lf", req_rel_error);
    return 0;
}

int
debug ()
{
    printf("DEBUG\n");
    double sigma = 1.0/1000;
    double epsr = 10.0;
    double mur = 1.0;
    int nf = 150;
    double *freq = malloc(nf * sizeof(double));
    logspace(1.0, 6.4, nf, freq);
    double radius = 7e-3;
    double depth = 0.5;
    int v1 = 3;
    int v2 = 3;
    int l1 = 3;
    int l2 = 3;
    double length1 = 10.0*(v1 - 1);
    double length2 = 10.0*(v2 - 1);
    int num_electrodes = l1*v2*(v1 - 1) + l2*v1*(v2 - 1);
    int num_nodes = v1*v2 + v1*(v2 - 1)*(l2 - 1) + v2*(v1 - 1)*(l1 - 1);
    Electrode *electrodes = malloc(num_electrodes * sizeof(Electrode));
    double *nodes = malloc(3 * num_nodes * sizeof(double));
    _Complex double zi = 0.0;
    electrode_grid(electrodes, nodes, num_nodes, radius, depth, zi, v1, length1,
                   l1, v2, length2, l2);
    Electrode *images = malloc(num_electrodes * sizeof(Electrode));
    for(size_t i = 0; i < num_electrodes; i++) {
        populate_electrode(&(images[i]), electrodes[i].start_point,
                           electrodes[i].end_point, electrodes[i].radius,
                           electrodes[i].zi);
        images[i].start_point[2] = -images[i].start_point[2];
        images[i].end_point[2] = -images[i].end_point[2];
    }
    size_t max_eval = 0;
    double req_abs_error = 1E-3;
    double req_rel_error = 1E-4;
    _Complex double *s = malloc(nf*sizeof(_Complex double));
    for(size_t i = 0; i < nf; i++) {
        s[i] = (_Complex double) I*TWO_PI*freq[i];
    }
    _Complex double *zh = malloc(num_nodes * nf * sizeof(_Complex double));
    zh_immittance(nf, s, sigma, epsr, mur, electrodes, images, num_electrodes,
                  nodes, num_nodes, max_eval, req_abs_error, req_rel_error, zh);
    for (int i = 0; i < nf; i++) {
        printf("%f\n", cabs(zh[i*num_nodes]));
    }
    // TODO output IT and IL
    free(freq);
    free(electrodes);
    free(images);
    free(s);
    free(zh);
    printf("END DEBUG\n");
    return 0;
}

int
zh_grid ()
{
    printf("\nSimulating a rectangular grounding grid with (v1 x v2) vertices.\n");
    printf("The harmonic impedance will be calculated and output saved\n");
    char file_name[50] = "gs_L1xL2.dat";
    printf("to file '%s' where each line is a frequency point\n", file_name);
    printf("and each column an injection node.\n");
    printf("     1   .....   v1\n\n");
    printf("     o---o---o---o  1\n");
    printf("     |   |   |   |  .\n");
    printf(" L2  o---o---o---o  .\n");
    printf("     |   |   |   |  .\n");
    printf("     o---o---o---o  v2\n");
    printf("          L1\n");
    double sigma, epsr, mur;
    soil_param(&sigma, &epsr, &mur);
    int nf;
    printf("Enter the number of frequency points: ");
    scanf("%i", &nf);
    double *freq = malloc(nf * sizeof(double));
    if (nf > 1) {
        freq_param(nf, freq);
    } else {
        printf("Enter frequency [Hz]: ");
        scanf("%lf", freq);
    }
    double radius, frac, length1, length2, depth;
    int num_electrodes, num_nodes, v1, l1, v2, l2;
    grid_param(&radius, &depth, &frac, &num_electrodes, &num_nodes,
               &v1, &length1, &l1, &v2, &length2, &l2,
               freq[nf-1], sigma, epsr, mur);
    Electrode *electrodes = malloc(num_electrodes * sizeof(Electrode));
    double *nodes = malloc(3 * num_nodes * sizeof(double));
    _Complex double zi = 0.0;
    electrode_grid(electrodes, nodes, num_nodes, radius, depth, zi, v1, length1,
                   l1, v2, length2, l2);
    Electrode *images = malloc(num_electrodes * sizeof(Electrode));
    for(size_t i = 0; i < num_electrodes; i++) {
        populate_electrode(&(images[i]), electrodes[i].start_point,
                           electrodes[i].end_point, electrodes[i].radius,
                           electrodes[i].zi);
        images[i].start_point[2] = -images[i].start_point[2];
        images[i].end_point[2] = -images[i].end_point[2];
    }
    size_t max_eval;
    double req_abs_error, req_rel_error;
    intg_param(&max_eval, &req_abs_error, &req_rel_error);
    printf("Simulating...\n");
    _Complex double *s = malloc(nf*sizeof(_Complex double));
    for(size_t i = 0; i < nf; i++) {
        s[i] = (_Complex double) I*TWO_PI*freq[i];
    }
    _Complex double *zh = malloc(num_nodes * nf * sizeof(_Complex double));
    // TODO sort nodes
    zh_immittance(nf, s, sigma, epsr, mur, electrodes, images, num_electrodes,
                  nodes, num_nodes, max_eval, req_abs_error, req_rel_error, zh);
    // TODO specialized calculation: identify 1/4 of the grid and inject
    // currents only on them (symmetry)
    // save results to file
    sprintf(file_name, "gs_%.2fx%.2f.dat", (v1 - 1)*length1, (v2 - 1)*length2);
    remove(file_name);
    FILE *save_file = fopen(file_name, "w");
    if (save_file == NULL) {
        printf("Cannot open file %s\n",  file_name);
        exit(1);
    }
    for (int i = 0; i < nf; i++) {
        v1 = i*num_nodes;
        for (int k = 0; k < (num_nodes - 1); k++) {
            fprintf(save_file, "%f, ", cabs(zh[k + v1]));
        }
        fprintf(save_file, "%f\n", cabs(zh[(num_nodes - 1) + v1]));
    }
    fclose(save_file);
    free(freq);
    free(electrodes);
    free(images);
    free(s);
    free(zh);
    printf("Simulation ended.\n");
    return 0;
}

int main (int argc, char *argv[])
{
    printf("HP_HEM %s\n", VERSION);
    printf("High Performance implementation of the Hybrid Electromagnetic Model\n");
    printf("https://github.com/pedrohnv/hp_hem\n");
    printf("Author: Pedro Henrique Nascimento Vieira\n");
    printf("pedrohnv@hotmail.com\n\n");
    if (argc == 5) {
        // params, elec, ne, out
        printf("Sorry, not implemented yet.\n");
    } else if (argc == 11) {
        // nf, (1 | 2), fmin, fmax, rho,
        // epsr, mur, max_eval, abs_err, rel_err,
        // ne, elec.txt, nn, nodes.txt, inj_node,
        // out
        int ns = atoi(argv[1]);
        int spac = atoi(argv[2]);
        double fmin = atof(argv[3]);
        double fmax = atof(argv[4]);
        double *freq = malloc(ns*sizeof(double));
        if (spac == 1) {
            linspace(fmin, fmax, ns, freq);
        } else if (spac == 2) {
            logspace(fmin, fmax, ns, freq);
        } else {
            printf("Invalid frequency spacing argument: %i\n", spac);
            printf("Aborting...\n");
            exit(2);
        }
        _Complex double *s = malloc(ns*sizeof(_Complex double));
        for(size_t i = 0; i < ns; i++) {
            s[i] = (_Complex double) I*TWO_PI*freq[i];
        }
        free(freq);
        double sigma = 1.0/atof(argv[5]);
        double epsr = atof(argv[6]);
        double mur = atof(argv[7]);
        size_t max_eval = atoi(argv[8]);
        double req_abs_error = atof(argv[9]);
        double req_rel_error = atof(argv[10]);
        size_t num_electrodes = (size_t) atoi(argv[11]);
        Electrode *electrodes = malloc(num_electrodes*sizeof(Electrode));
        electrodes_file(argv[12], electrodes, num_electrodes);
        size_t num_nodes = atoi(argv[13]);
        double *nodes = malloc(3 * num_nodes * sizeof(double));
        nodes_file(argv[14], nodes, num_nodes);
        size_t inj_node = atoi(argv[15]);
        _Complex double *inj_current = malloc(ns*sizeof(_Complex double));
        for (size_t i = 0; i < ns; i++) inj_current[i] = 1.0;
        char *out = argv[17];
        FILE *save_file = fopen(out, "w");
        if (save_file == NULL) {
            printf("Cannot open file %s\n", out);
            exit(1);
        }
        Electrode *images = malloc(num_electrodes * sizeof(Electrode));
        for(size_t i = 0; i < num_electrodes; i++) {
            populate_electrode(&(images[i]), electrodes[i].start_point,
                               electrodes[i].end_point, electrodes[i].radius,
                               electrodes[i].zi);
            images[i].start_point[2] = -images[i].start_point[2];
            images[i].end_point[2] = -images[i].end_point[2];
        }
        _Complex double *u = malloc((ns * num_nodes) * sizeof(_Complex double));
        _Complex double *il = malloc((ns * num_electrodes) * sizeof(_Complex double));
        _Complex double *it = malloc((ns * num_electrodes) * sizeof(_Complex double));
        _Complex double *inj_adm = calloc(ns, sizeof(_Complex double));
        sim_immittance(ns, s, sigma, epsr, mur, electrodes, images,
                       num_electrodes, nodes, num_nodes, max_eval, req_abs_error,
                       req_rel_error, inj_node, inj_current, inj_adm, u, il, it);
        free(electrodes);
        fclose(save_file);
    } else {
        printf("No correct input detected. argc = %i\n", argc);
        printf("Expected 3 files: 'params.txt, electrodes.txt, output.dat'\n");
        printf("or parameters list and 2 files:\n");
        printf("'nf, (1 | 2), fmin, fmax, rho, epsr, mur, electrodes.txt, output.dat'\n");
        zh_grid();
        //debug();
    }
    return 0;
}
