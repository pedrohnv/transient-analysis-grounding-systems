/*
Functions for testing the code
*/

#include "cubature.h"
#include "electrode.h"
#include "auxiliary.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
//#include <omp.h>

/* Auxiliary routine: printing a matrix */
void
print_matrix (char *desc, MKL_INT m, MKL_INT n, _Complex double *a, MKL_INT lda)
{
    MKL_INT i, j;
    printf( "\n %s\n", desc );
    for( i = 0; i < m; i++ ) {
        for( j = 0; j < n; j++ )
            printf( " (%6.2f,%6.2f)", creal(a[i*lda+j]), cimag(a[i*lda+j]) );
        printf( "\n" );
    }
}

int
we_building ()
{
    /* electrode system:
    3---4
    | / |
    0---1---2
    */
    double radius = 1e-3;
    _Complex double zi = 0.0;
    size_t ne = 6;
    Electrode electrodes[ne];
    size_t nn = 5;
    double nodes[nn][3];
    nodes[0][0] = 0.0;
    nodes[0][1] = 0.0;
    nodes[0][2] = 0.0;

    nodes[1][0] = 1.0;
    nodes[1][1] = 0.0;
    nodes[1][2] = 0.0;

    nodes[2][0] = 2.0;
    nodes[2][1] = 0.0;
    nodes[2][2] = 0.0;

    nodes[3][0] = 0.0;
    nodes[3][1] = 1.0;
    nodes[3][2] = 0.0;

    nodes[4][0] = 1.0;
    nodes[4][1] = 1.0;
    nodes[4][2] = 0.0;
    populate_electrode(&(electrodes[0]), nodes[0], nodes[1], radius, zi);
    populate_electrode(&(electrodes[1]), nodes[1], nodes[2], radius, zi);
    populate_electrode(&(electrodes[2]), nodes[0], nodes[3], radius, zi);
    populate_electrode(&(electrodes[3]), nodes[1], nodes[4], radius, zi);
    populate_electrode(&(electrodes[4]), nodes[0], nodes[4], radius, zi);
    populate_electrode(&(electrodes[5]), nodes[3], nodes[4], radius, zi);
    _Complex double we[(2*ne + nn)*(2*ne + nn)];
    fill_incidence(we, electrodes, ne, nodes, nn);
    _Complex double zl[ne*ne], zt[ne*ne], yn[nn*nn];
    for (size_t i = 0; i < (ne*ne); i++) {
        zl[i] = 3.0;
        zt[i] = 2.0;
    }
    for (size_t i = 0; i < (nn*nn); i++) {
        yn[i] = 4.0;
    }
    fill_impedance(we, electrodes, ne, nn, zl, zt, yn);
    _Complex double we_target[] =
        {1.5, 1.5, 1.5, 1.5, 1.5, 1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5,
        -1., 1., 0., 0., 0., 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, -1.5, -1.5, -1.5,
        -1.5, -1.5, -1.5, 0., -1., 1., 0., 0., 1.5, 1.5, 1.5, 1.5, 1.5, 1.5,
        -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1., 0., 0., 1., 0., 1.5, 1.5,
        1.5, 1.5, 1.5, 1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, 0., -1., 0.,
        0., 1., 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, -1.5, -1.5, -1.5, -1.5, -1.5,
        -1.5, -1., 0., 0., 0., 1., 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, -1.5, -1.5,
        -1.5, -1.5, -1.5, -1.5, 0., 0., 0., -1., 1., 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, -0.5, -0.5, 0., 0., 0., 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 0., -0.5, -0.5, 0., 0., 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, -0.5,
        0., 0., -0.5, 0., 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0., -0.5, 0.,
        0., -0.5, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, -0.5, 0., 0., 0., -0.5,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0., 0., 0., -0.5, -0.5, 1., 0.,
        1., 0., 1., 0., 0., 0., 0., 0., 0., 0., 4., 4., 4., 4., 4., 0., 1.,
        0., 1., 0., 0., 1., 0., 0., 0., 0., 0., 4., 4., 4., 4., 4., 0., 0.,
        0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 4., 4., 4., 4., 4., 0., 0.,
        0., 0., 0., 1., 0., 0., 1., 0., 0., 0., 4., 4., 4., 4., 4., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 4., 4., 4., 4., 4.};
    //=====================================================
    assert(sizeof(we) == sizeof(we_target));
    // print WE ?
    if ( 0 ) {
        printf("WE\n");
        for (size_t i = 0; i < (2*ne + nn); i++) {
            for (size_t k = 0; k < (2*ne + nn); k++) {
                printf("% .1f  ", creal(we[i*(2*ne + nn) + k]));
            }
            printf("\n");
        }
        printf("\nWE_target\n");
        for (size_t i = 0; i < (2*ne + nn); i++) {
            for (size_t k = 0; k < (2*ne + nn); k++) {
                printf("% .1f  ", creal(we_target[i*(2*ne + nn) + k]));
            }
            printf("\n");
        }
    }
    //=====================================================
    size_t a, b, c, d;
    size_t jump = (2*ne + nn)*ne;
    size_t jump_2 = 2*jump;
    for (size_t e = 0; e < ne; e++) {
        for (size_t n = 0; n < nn; n++) {
            a = 2*ne*(e + 1) + e*nn + n;
            b = jump + a;
            c = jump_2 + n*(2*ne + nn) + e;
            d = ne + c;
            assert( creal(we_target[a]) == creal(we[a]) );
            assert( creal(we_target[b]) == creal(we[b]) );
            assert( creal(we_target[c]) == creal(we[c]) );
            assert( creal(we_target[d]) == creal(we[d]) );
        }
    }
    printf("we building: incidence matrices: passed\n");
    //=====================================================
    for (size_t i = 0; i < (2*ne + nn)*(2*ne + nn); i++) {
        assert(we[i] == we_target[i]);
    }
    printf("we building: impedance and admitance: passed\n");
    return 0;
}

int
segmentation_example ()
{
    // TODO write on doc
    printf("3D linspace\n");
    double length = 35.0;
    double pf[] = {21.0025, 25.41, 11.757};
    double radius = 1e-2;
    size_t n = round(length/(10*radius));
    double xx[n];
    linspace(0.0, pf[0], n, xx);
    double yy[n];
    linspace(0.0, pf[1], n, yy);
    double zz[n];
    linspace(0.0, pf[2], n, zz);
    for (size_t i = 0; i < n; i++) {
        printf("x = %f | y = %f | z = %f\n", xx[i], yy[i], zz[i]);
    }
    return 0;
}

int
test_case (double rho, double length, double frac)
{
    // setup
    size_t nf = 250;
    double freq[nf];
    logspace(2.0, 7.5, nf, freq);
    double h = 0.001;
    _Complex double zi = 0.0;
    double radius = 1.25e-2;
    double mur = 1.0;
    double eps = 10.0*EPS0;
    double sigma = 1/rho;
    double rho_c = 1.9e-6;
    // execute
    remove("test_case.dat");
    FILE *ftest_case = fopen("test_case.dat", "w");
    if (ftest_case == NULL) {
        printf("Cannot open file %s\n", "test_case.dat");
        exit(1);
    }
    size_t nn2, ne2, n2en1, n2en2, num_nodes, num_electrodes;
    num_nodes = round(length/(frac*radius));
    num_electrodes = num_nodes - 1;
    ne2 = num_electrodes*num_electrodes;
    nn2 = num_nodes*num_nodes;
    n2en1 = (2*num_electrodes + num_nodes);
    n2en2 = n2en1*n2en1;
    _Complex double k1, c, zl[ne2], zt[ne2], yn[nn2];
    // segment electrode
    double nodes[num_nodes][3], nodes_image[2][3];
    double segment[num_nodes];
    linspace(-h, -h-length, num_nodes, segment);
    for (size_t i = 0; i < num_nodes; i++) {
        nodes[i][0] = 0.0;
        nodes[i][1] = 0.0;
        nodes[i][2] = segment[i];
    }
    for ( size_ti = 0; i < nn2; i++) {
        yn[i] = 0.0;
    }
    Electrode electrodes[num_electrodes];
    Electrode images[num_electrodes];
    for (size_t i = 0; i < num_electrodes; i++) {
        populate_electrode(
            &(electrodes[i]), nodes[i], nodes[i + 1], radius, zi);
        nodes_image[0][0] = nodes[i][0];
        nodes_image[0][1] = nodes[i][1];
        nodes_image[0][2] = -nodes[i][2];
        nodes_image[1][0] = nodes[i + 1][0];
        nodes_image[1][1] = nodes[i + 1][1];
        nodes_image[1][2] = -nodes[i + 1][2];
        populate_electrode(
            &(images[i]), nodes_image[0], nodes_image[1], radius, zi);
    }
    _Complex double rhs[n2en1], incidence[n2en2], we[n2en2], ie[n2en1];
    for (size_t i = 0; i < n2en1; i++) {
        rhs[i] = 0.0;
    }
    rhs[n2en1 - num_nodes] = 1.0;
    fill_incidence(incidence, electrodes, num_electrodes, nodes, num_nodes);
    _Complex double s;
    _Complex double rt, rl; // reflection coefficients
    rt = 1.0; // TODO calculate for each frequency
    rl = 1.0;
    _Complex double zinternal;
    double l = electrodes[0].length;
    // solve for each frequency
    for (size_t i = 0; i < nf; i++) {
        //printf("\nf = %.2f Hz\n", freq[i]);
        s = I*TWO_PI*freq[i];
        memcpy(ie, rhs, sizeof(ie));
        memcpy(we, incidence, sizeof(we));
        zinternal = internal_impedance(s, rho_c, radius, mur)*l;
        for (size_t k = 0; k < num_electrodes; k++) {
            electrodes[k].zi = zinternal;
        }
        c = (sigma + s*eps);
        k1 = csqrt(s*mur*MU0*c); //gamma
        calculate_impedances(electrodes, num_electrodes, zl, zt, k1, s, mur, c,
            200, 1e-3, 1e-4, ERROR_PAIRED, INTG_DOUBLE);
        impedances_images(
            electrodes, images, num_electrodes, zl, zt, k1, s, mur, c, rt, rl,
            200, 1e-3, 1e-4, ERROR_PAIRED, INTG_DOUBLE);
        fill_impedance(we, electrodes, num_electrodes, num_nodes, zl, zt, yn);
        solve_electrodes(we, ie, num_electrodes, num_nodes);
        fprintf(ftest_case, "%f %f\n",
            creal(ie[n2en1 - num_nodes]), cimag(ie[n2en1 - num_nodes]));
    }
    //print_matrix( "\nWE", n2en1, n2en1, we, n2en1 );
    //print_matrix( "\nZL", num_electrodes, num_electrodes, zl, num_electrodes );
    //print_matrix( "\nZT", num_electrodes, num_electrodes, zt, num_electrodes );
    fclose(ftest_case);
    return 0;
}

int
zi_bessel ()
{
    _Complex double s = I*377;
    double rho = 1.9e-6;
    double radius = 7e-3;
    _Complex double zi = internal_impedance(s, rho, radius, 1.0);
    assert(abs(creal(zi) - 0.0123426) < 1e-6);
    assert(abs(cimag(zi) - 0.00001885) < 1e-6);
    printf("Internal impedance: passed\n");
    return 0;
}

int
integration ()
{
    double h = 0.3;
    double len = 3.0;
    double r1 = 12.5e-3;
    double sigma = 1.0/100.0;
    double er = 10.0;
    // frequencies of interest
    size_t nf = 150;
    double freq[nf];
    logspace(2, 7, nf, freq);
    char file_name[] = "integration_test.dat";
    remove(file_name);
    FILE *save_file = fopen(file_name, "w");
    if (save_file == NULL) {
        printf("Cannot open file %s\n",  file_name);
        exit(1);
    }
    double start_point[3] = {0., 0., h};
    double end_point[3] = {len, 0., h};
    Electrode *sender = (Electrode*) malloc(sizeof(Electrode));
    populate_electrode(sender, start_point, end_point, r1, 0.0);
    start_point[1] = len;
    start_point[2] = -h;
    end_point[2] = -h;
    Electrode *receiver = (Electrode*) malloc(sizeof(Electrode));
    populate_electrode(receiver, start_point, end_point, r1, 0.0);
    double result[2], error[2];
    _Complex double s;
    _Complex double kappa, k1;
    for (size_t i = 0; i < nf; i++) {
        s = I*TWO_PI*freq[i];
        kappa = (sigma + s*er*EPS0); //soil complex conductivity
        k1 = csqrt(s*MU0*kappa); //soil propagation constant
        integral(sender, receiver, k1, 200, 1e-3, 1e-4,
            ERROR_PAIRED, INTG_DOUBLE, result, error);
        fprintf(save_file, "%f %f\n", result[0], result[1]);
    }
    fclose(save_file);
    free(sender);
    free(receiver);
    return 0;
}

int
main (int argc, char *argv[])
{
    we_building();
    zi_bessel();
    integration();
    //test_case(30., 3., 40.); //prtl-gnd-test10, grcevL3rho30
}
