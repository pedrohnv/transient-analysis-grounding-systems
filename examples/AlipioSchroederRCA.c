/*
Test case AlipioSchroederRCA

Reproducing the results in [1].

[1] ALÍPIO, Rafael Silva et al. Modelagem de Aterramentos Elétricos para Fenômenos de Alta Frequência e Comparação com Resultados Experimentais. Revista Controle & Automação (SBA), v. 22, n. 1, p. 89-102, 2011.
*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cubature.h>
#include <Electrode.h>
#include <auxiliary.h>
//#include <omp.h>

int run_case(double length, double rho, char file_name[],
    _Complex double* ifonte, _Complex double* w, int nf)
{
    // parameters
    double h = 0.6; //burial depth
    double r = 0.00607651; //radius
    double sigma = 1/rho; // soil conductivity
    double er = 15.0; //soil rel. permitivitty
    double rho_c = 1.9e-6; //copper resistivity
    // frequencies of interest
    //int nf = 250;
    //double freq[nf];
    double start_point[3] = {0., 0., -h};
    double end_point[3] = {length, 0., -h};

    remove(file_name);
    FILE* save_file = fopen(file_name, "w");
    if (save_file == NULL)
    {
        printf("Cannot open file %s\n",  file_name);
        exit(1);
    }
    //logspace(2, 7, nf, freq);

    // electrode definition and segmentation
    double lambda = wave_length(-cimag(w[nf - 1])/TWO_PI, sigma, er*EPS0, MU0); //smallest
    int num_electrodes = ceil( length/(lambda/6.0) ) ;
    int num_nodes = num_electrodes + 1;
    double nodes[num_nodes][3];
    Electrode* electrodes = (Electrode*) malloc(sizeof(Electrode)*num_electrodes);
    // the internal impedance is added "outside" later
    segment_electrode(
        electrodes, nodes, num_electrodes, start_point, end_point, r, 0.0);

    // create images
    start_point[2] = h;
    end_point[2] = h;
    double nodes_images[num_nodes][3];
    //Electrode images[num_electrodes];
    Electrode* images = (Electrode*) malloc(sizeof(Electrode)*num_electrodes);
    segment_electrode(
        images, nodes_images, num_electrodes, start_point, end_point, r, 0.0);

    //build system to be solved
    int ne2 = num_electrodes*num_electrodes;
    int nn2 = num_nodes*num_nodes;
    int ss1 = (2*num_electrodes + num_nodes);
    int ss2 = ss1*ss1;
    //double w;
    _Complex double kappa, gamma, zinternal;
    _Complex double* zl = (_Complex double*) malloc(sizeof(_Complex double)*ne2);
    _Complex double* zt = (_Complex double*) malloc(sizeof(_Complex double)*ne2);
    _Complex double* yn = (_Complex double*) malloc(sizeof(_Complex double)*nn2);
    _Complex double* ie = (_Complex double*) malloc(sizeof(_Complex double)*ss1);
    _Complex double* ie_cp = (_Complex double*) malloc(sizeof(_Complex double)*ss1);
    _Complex double* we = (_Complex double*) malloc(sizeof(_Complex double)*ss2);
    _Complex double* we_cp = (_Complex double*) malloc(sizeof(_Complex double)*ss2);
    int i, k;
    for (i = 0; i < nn2; i++)
    {
        yn[i] = 0.0; //external nodal admittance
    }
    for (i = 0; i < ss1; i++)
    {
        ie[i] = 0.0;
    }
    fill_incidence(we, electrodes, num_electrodes, nodes, num_nodes);
    // solve for each frequency: WE*VE = IE
    for (i = 0; i < nf; i++)
    {
        ie[ss1 - num_nodes] = ifonte[i];
        kappa = (sigma + I*w[i]*er*EPS0); //soil complex conductivity
        gamma = csqrt(I*w[i]*MU0*kappa); //soil propagation constant
        //TODO especialized impedance calculation taking advantage of symmetry
        calculate_impedances(
            electrodes, num_electrodes, zl, zt, gamma, w[i], MU0, kappa,
            200, 1e-3, 1e-4, ERROR_PAIRED, INTG_DOUBLE);
        zinternal = internal_impedance(w[i], rho_c, r, MU0)*electrodes[0].length;
        for (k = 0; k < num_electrodes; k++)
        {
            zl[k*num_electrodes + k] += zinternal;
        }
        impedances_images(electrodes, images, num_electrodes, zl, zt, gamma,
            w[i], MU0, kappa, 0.0, 1.0, 200, 1e-3, 1e-4, ERROR_PAIRED, INTG_DOUBLE);
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

int main()
{
    double t[] = {9.22133e-10, 2.80886e-8, 4.45412e-8, 5.37883e-8,
            6.85252e-8, 7.96766e-8, 9.26036e-8, 1.01911e-7,
            1.13054e-7, 1.33255e-7, 1.51577e-7, 1.69994e-7,
            1.93652e-7, 2.17422e-7, 2.46467e-7, 2.75426e-7,
            3.18771e-7, 3.54996e-7, 3.82102e-7, 4.01951e-7,
            4.30808e-7, 4.7423e-7, 5.14117e-7, 5.53919e-7,
            5.8643e-7, 6.28093e-7, 6.6967e-7, 7.02138e-7,
            7.32907e-7, 7.92756e-7, 8.57914e-7, 8.92286e-7,
            9.30355e-7, 9.79207e-7, 9.79207e-5};
    double inj[] = {0.0473934, 0.236967, 2.03791, 4.21801, 7.06161, 10.2844,
            13.128, 15.9716, 19.0995, 22.3223, 24.7867, 28.2938, 29.7156,
            32.3697, 33.3175, 33.3175, 32.2749, 32.5592, 32.0853, 31.4218,
            30.2844, 30.0948, 30.8531, 30.6635, 29.9052, 30.2844, 29.7156,
            28.4834, 28.4834, 29.8104, 29.8104, 29.6209, 30.2844, 30.0948,
            30.0948};
    int nt = sizeof(t)/sizeof(double);
    int ns = 1024/2; //half the frequency sample (because symmetry)
    _Complex double* s = malloc(sizeof(_Complex double)*ns);
    double dt = t[nt - 1]/(2*ns);
    double dw = TWO_PI/(2*ns*dt);
    double sigma = log(0.001)/t[nt - 1];
    for (int k = 0; k < ns; k++)
    {
        s[k] = (sigma - I*dw*k);
    }
    _Complex double* ifonte = malloc(sizeof(_Complex double)*ns);
    laplace_transform(inj, t, nt, s, ns, ifonte);
    run_case(8.0, 65.0, "examples/AlipioSchroederRCA_1.dat", ifonte, s, ns);
    //run_case(15.0, 70.0, "examples/AlipioSchroederRCA_2.dat", ifonte);
    free(s);
    free(ifonte);
    return 0;
}
