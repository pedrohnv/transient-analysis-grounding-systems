/*
Reproducing some results for the electric field on ground level for an
horizontal electrode with frequency independent of the soil parameters.

[1] R.S. Alipio, M.A.O. Schroeder, M.M. Afonso, T.A.S. Oliveira, S.C. Assis,
Electric fields of grounding electrodes with frequency dependent soil parameters,
Electric Power Systems Research,
Volume 83, Issue 1, 2012, Pages 220-226, ISSN 0378-7796,
https://doi.org/10.1016/j.epsr.2011.11.011.
(http://www.sciencedirect.com/science/article/pii/S0378779611002781)
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include <omp.h>
#include "auxiliary.h"
#include "electrode.h"
#include "linalg.h"
#include "cubature.h"
#include "grid.h"

int
run_case (double lmax, double sigma0)
{
    const unsigned NRHS = 1;
    double mur = 1.0;  // soil rel. magnetic permeability
    // numerical integration parameters
    size_t max_eval = 0;
    double req_abs_error = 1e-5;
    double req_rel_error = 1e-6;

    // electrodes
    double r = 7e-3;
    double length = 15.0;
    double h = 1.0;
    double x0 = 2.5;
    double start_point[3] = {x0, 0., -h};
    double end_point[3] = {x0 + length, 0., -h};
    int ne = ceil(length / lmax) + 1;
    int nn = ne + 1;
    printf("Num. segments = %i\n", ne);
    printf("Num. nodes = %i\n", nn);
    double* nodes = malloc(3 * nn * sizeof(double));
    Electrode* electrodes = malloc(ne * sizeof(Electrode));
    segment_electrode(electrodes, nodes, ne, start_point, end_point, r);
    Electrode* images =malloc(ne * sizeof(Electrode));
    for (size_t m = 0; m < ne; m++) {
        populate_electrode(images + m, electrodes[m].start_point,
                           electrodes[m].end_point, electrodes[m].radius);
        images[m].start_point[2] = -images[m].start_point[2];
        images[m].end_point[2] = -images[m].end_point[2];
        images[m].middle_point[2] = -images[m].middle_point[2];
    }
    size_t inj_node = 0;

    // frequencies
    double freq[] = {100.0, 500e3, 1e6, 2e6};
    unsigned ns = 4;
    _Complex double s[ns];
    for (unsigned i = 0; i < ns; i++) s[i] = I * TWO_PI * freq[i];

   // define an array of points where to calculate quantities
   // line along X axis from 0 to 20
   double dr = 0.1;
   size_t num_points = ceil(20.0 / dr) + 1;
   printf("Num. points to calculate GPD = %li\n", num_points);
   double* points = malloc(num_points * 3 * sizeof(double));
   for (size_t i = 0; i < num_points; i++) {
       points[3*i + 0] = dr * i;
       points[3*i + 1] = 0.0;
       points[3*i + 1] = 0.0;
   }

    // malloc matrices =========================================================
    size_t ne2 = ne * ne;
    size_t nn2 = nn * nn;
    _Complex double* potzl = malloc(ne2 * sizeof(_Complex double));
    _Complex double* potzt = malloc(ne2 * sizeof(_Complex double));
    _Complex double* potzli = malloc(ne2 * sizeof(_Complex double));
    _Complex double* potzti = malloc(ne2 * sizeof(_Complex double));
    _Complex double* a = malloc((ne * nn) * sizeof(_Complex double));
    _Complex double* b = malloc((ne * nn) * sizeof(_Complex double));

    _Complex double* gpr_s = malloc(NRHS * ns * sizeof(_Complex double));
    _Complex double* ground_pot_s = calloc(NRHS * num_points * ns, sizeof(_Complex double));
    size_t np_field = num_points;
    _Complex double* efield_s = malloc(NRHS * 6 * np_field * ns * sizeof(_Complex double));
    // efield[i,j,k,m] => efield[m + (k * 6) + (j * 6 * NRHS) + (i * 6 * NRHS * np_field)]
    // dimensions: i-frequency, j-point, k-injection, m-field
    _Complex double* voltage_s = malloc(NRHS * np_field * ns * sizeof(_Complex double));

    // Incidence and "z-potential" (mHEM) matrices =============================
    int err;
    err = fill_incidence_adm(a, b, electrodes, ne, nodes, nn);
    if (err != 0) printf("Could not build incidence matrices\n");
    err = calculate_impedances(potzl, potzt, electrodes, ne, 0.0, 0.0, 0.0, 0.0,
                               max_eval, req_abs_error, req_rel_error, INTG_MHEM);
    if (err != 0) printf("integration error\n");
    err = impedances_images(potzli, potzti, electrodes, images, ne,
                            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, max_eval,
                            req_abs_error, req_rel_error, INTG_MHEM);
    if (err != 0) printf("integration error\n");

    // BLAS and LAPACK helper variables
    char trans = 'N';
    _Complex double one = 1.0;
    _Complex double zero = 0.0;
    double begin = omp_get_wtime();  // to estimate time until completion
    #pragma omp parallel private(err)
    {
        #pragma omp single
        {
            printf("avaible threads: %i\n", omp_get_num_threads());
        }
        _Complex double* zl = malloc(ne2 * sizeof(_Complex double));
        _Complex double* zt = malloc(ne2 * sizeof(_Complex double));
        _Complex double* yn = malloc(nn2 * sizeof(_Complex double));
        _Complex double* yla = malloc((ne * nn) * sizeof(_Complex double));
        _Complex double* ytb = malloc((ne * nn) * sizeof(_Complex double));
        _Complex double* ie = malloc(NRHS * nn * sizeof(_Complex double));
        _Complex double* il = malloc(NRHS * ne * sizeof(_Complex double));
        _Complex double* it = malloc(NRHS * ne * sizeof(_Complex double));
        if (!zl || !zt || !ie || !yn || !yla || !ytb || !ie || !il || !it) {
            printf("Can't allocate memory\n");
            exit(200);
        }
        double field_point[] = {0.0, 0.0, 0.0};
        double field_point2[] = {0.0, 0.0, 0.0};
        _Complex double sigma, epsr, kappa, gamma;  // soil parameters
        _Complex double ref_l, ref_t;  // reflection coefficients (images)
        _Complex double iwu_4pi, one_4pik, exp_gr;
        double rbar, r0, r1, r2;
        //_Complex double gpd_var[NRHS];
        _Complex double efvector[3];  // helper variable of a 3d vector
        int index;

        // use raw BLAS and LAPACK to solve system with multiple RHS
        int nn1 = (int) nn;
        int ldie = nn1;  // leading dimension of ie
        int* ipiv = malloc(nn1 * sizeof(int));  // pivot indices
        int info;
        int nrhs = NRHS;
        char uplo = 'L';  // matrices are symmetric, only Lower Half is set
        int lwork = -1;  // signal to Query the optimal workspace
        _Complex double wkopt;
        zsysv_(&uplo, &nn1, &nrhs, yn, &nn1, ipiv, ie, &ldie, &wkopt, &lwork, &info);
        lwork = creal(wkopt);
        _Complex double* work = malloc(lwork * sizeof(_Complex double));
        #pragma omp for
        for (size_t i = 0; i < ns; i++) {
            //printf("i = %li from thread %d\n", i, omp_get_thread_num());
            sigma = sigma0;
            epsr = 4.0;
            kappa = (sigma + s[i] * epsr * EPS0);  // soil complex conductivity
            gamma = csqrt(s[i] * MU0 * kappa);  // soil propagation constant
            iwu_4pi = s[i] * mur * MU0 / (FOUR_PI);
            one_4pik = 1.0 / (FOUR_PI * kappa);
            // reflection coefficient, soil to air
            ref_t = (kappa - s[i] * EPS0) / (kappa + s[i] * EPS0);
            ref_l = 1.0;
            // modified HEM (mHEM):
            for (size_t m = 0; m < ne; m++) {
                for (size_t k = m; k < ne; k++) {
                    rbar = vector_length(electrodes[k].middle_point,
                                         electrodes[m].middle_point);
                    exp_gr = cexp(-gamma * rbar);
                    zl[m * ne + k] = exp_gr * potzl[m * ne + k];
                    zt[m * ne + k] = exp_gr * potzt[m * ne + k];
                    rbar = vector_length(electrodes[k].middle_point,
                                         images[m].middle_point);
                    exp_gr = cexp(-gamma * rbar);
                    zl[m * ne + k] += ref_l * exp_gr * potzli[m * ne + k];
                    zt[m * ne + k] += ref_t * exp_gr * potzti[m * ne + k];
                    zl[m * ne + k] *= iwu_4pi;
                    zt[m * ne + k] *= one_4pik;
                }
            }
            // Traditional HEM (highly discouraged):
            /*calculate_impedances(zl, zt, electrodes, ne, gamma, s[i], 1.0, kappa,
                                 max_eval, req_abs_error, req_rel_error, INTG_DOUBLE);
            impedances_images(zl, zt, electrodes, images, ne, gamma, s[i], mur,
                              kappa, ref_l, ref_t, max_eval, req_abs_error,
                              req_rel_error, INTG_DOUBLE);*/

            calculate_yla_ytb(yla, ytb, zl, zt, a, b, ne, nn);
            // YN = A^T * inv(ZL) * A + B^T * inv(ZT) * B
            fill_impedance_adm2(yn, yla, ytb, a, b, ne, nn);
            for (size_t m = 0; m < (NRHS * nn); m++) {
                ie[m] = 0.0;
            }
            for (size_t m = 0; m < (NRHS); m++) {
                ie[inj_node + ldie * m] = 1.0;
            }
            // solve
            zsysv_(&uplo, &nn1, &nrhs, yn, &nn1, ipiv, ie, &ldie, work, &lwork, &info);
            // Check for the exact singularity
            if (info > 0) {
                printf("The diagonal element of the triangular factor of YN,\n");
                printf("U(%i,%i) is zero, so that YN is singular;\n", info, info);
                printf("the solution could not be computed.\n");
                printf("  in the %li-th frequency\n", i);
                exit(info);
            } else if (info < 0) {
                printf("the %i-th parameter to zsysv had an illegal value.\n", info);
            }
            err = 1;  // used in zgemv_ to specify incx, incy
            for (unsigned k = 0; k < NRHS; k++) {
                gpr_s[i + k * ns] = ie[inj_node + k * ldie];
                // IT = inv(ZT) * B * IE + 0 * IT
                zgemv_(&trans, &ne, &nn, &one, yla, &ne, (ie + nn1 * k), &err,
                       &zero, (il + ne * k), &err);
                // IL = inv(ZL) * A * IE + 0 * IL
                zgemv_(&trans, &ne, &nn, &one, ytb, &ne, (ie + nn1 * k), &err,
                       &zero, (it + ne * k), &err);
                // if using intel MKL, replace the above by:
                /*cblas_zgemv(CblasColMajor, CblasNoTrans, ne, nn, &one, ytb, ne,
                            (ie + nn1 * k), 1, &zero, (it + ne * k), 1);
                cblas_zgemv(CblasColMajor, CblasNoTrans, ne, nn, &one, yla, ne,
                            (ie + nn1 * k), 1, &zero, (il + ne * k), 1);*/
            }
            // images' effect as (1 + ref) because we're calculating on ground level
            for (unsigned k = 0; k < NRHS; k++) {
                for (size_t m = 0; m < ne; m++) {
                    it[m + k * ne] *= (1.0 + ref_t);
                    il[m + k * ne] *= (1.0 + ref_l);
                }
            }
            // Ground Potential Distribution (GPD) =============================
            for (size_t p = 0; p < num_points; p++) {
                for (unsigned k = 0; k < NRHS; k++) {
                    index = (i * NRHS * num_points) + (p * NRHS) + k;
                    ground_pot_s[index] = 0.0;
                    // Traditional HEM:
                    /*ground_pot_s[index] = electric_potential((points + p*3),
                                                    electrodes, ne, (it + ne * k),
                                                    gamma, kappa, max_eval,
                                                    req_abs_error, req_rel_error);*/

                }
                for (size_t m = 0; m < ne; m++) {
                // calculate using a simplification to the numerical integrals (mHEM):
                    rbar = vector_length((points + p*3), electrodes[m].middle_point);
                    r1 = vector_length((points + p*3), electrodes[m].start_point);
                    r2 = vector_length((points + p*3), electrodes[m].end_point);
                    r0 = (r1 + r2 + electrodes[m].length) / (r1 + r2 - electrodes[m].length);
                    exp_gr = cexp(-gamma * rbar) * log(r0) / electrodes[m].length;
                    for (unsigned k = 0; k < NRHS; k++) {
                        index = (i * NRHS * num_points) + (p * NRHS) + k;
                        ground_pot_s[index] += one_4pik * it[m + ne * k] * exp_gr;
                    }
                }
            }
            // Electric Field ==================================================
            for (size_t p = 0; p < np_field; p++) {
                field_point[0] = points[3 * p];
                for (size_t k = 0; k < NRHS; k++) {
                    // conservative field
                    for (size_t m = 0; m < 3; m++) {
                        efvector[m] = 0.0;
                    }
                    // pass "s = 0.0" so that the non-conservative part is ignored
                    electric_field(field_point, electrodes, ne,
                                   (il + ne * k), (it + ne * k),
                                   gamma, 0.0, mur, kappa, max_eval,
                                   req_abs_error, req_rel_error, efvector);
                    for (size_t m = 0; m < 3; m++) {
                        index = m + (k * 6) + (p * 6 * NRHS) + (i * 6 * NRHS * np_field);
                        efield_s[index] = efvector[m];
                    }
                    // non-conservative field
                    for (size_t m = 0; m < 3; m++) {
                        efvector[m] = 0.0;
                    }
                    magnetic_potential(field_point, electrodes, ne,
                                      (il + ne * k), gamma, mur, max_eval,
                                      req_abs_error, req_rel_error, efvector);
                    for (size_t m = 3; m < 6; m++) {
                        index = m + (k * 6) + (p * 6 * NRHS) + (i * 6 * NRHS * np_field);
                        efield_s[index] = efvector[m - 3] * (-s[i]);
                    }
                    // voltage =================================================
                    index = k + (p * NRHS) + (i * NRHS * np_field);
                    field_point2[0] = field_point[0] + 1.0;
                    voltage_s[index] = voltage(field_point, field_point2, electrodes, ne,
                                               (il + ne * k), (it + ne * k),
                                               gamma, s[i], mur, kappa,
                                               max_eval, req_abs_error, req_rel_error);
                }
            }
            if (i == 0) {
                printf("Expected more time until completion of the frequency loop: %.2f min.\n",
                       (omp_get_wtime() - begin) * ns / 60.0 / omp_get_num_threads());
            }
        }
        free(zl);
        free(zt);
        free(yn);
        free(yla);
        free(ytb);
        free(ie);
        free(il);
        free(it);
        free(ipiv);
        free(work);
    } // end parallel
    free(potzl);
    free(potzt);
    free(potzli);
    free(potzti);
    free(a);
    free(b);
    free(electrodes);
    free(images);
    free(nodes);
    printf("Frequency loop ended. Saving results.\n");
    _Complex double var;
    // GPR  ==============================================================
    char gpr_file_name[60];
    sprintf(gpr_file_name, "alipio83powsys_gpr_%.0f.csv", 1 / sigma0);
    FILE* gpr_file = fopen(gpr_file_name, "w");
    fprintf(gpr_file, "f");
    for (unsigned k = 0; k < NRHS; k++) fprintf(gpr_file, ",i%d", k+1);
    fprintf(gpr_file, "\n");
    for (size_t i = 0; i < ns; i++) {
        fprintf(gpr_file, "%e", freq[i]);
        for (size_t k = 0; k < NRHS; k++) {
            var = gpr_s[i + ns * k];
            fprintf(gpr_file, ",%e%+e%s", creal(var), cimag(var), "im");
        }
        fprintf(gpr_file, "\n");
    }
    fclose(gpr_file);
    free(gpr_s);

    // GPD  ==============================================================
    size_t index;
    char gpd_file_name[60];
    sprintf(gpd_file_name, "alipio83powsys_gpd_%.0f.csv", 1 / sigma0);
    FILE* gpd_file = fopen(gpd_file_name, "w");
    fprintf(gpd_file, "f,x");
    for (unsigned k = 0; k < NRHS; k++) fprintf(gpd_file, ",i%d", k+1);
    fprintf(gpd_file, "\n");
    for (size_t i = 0; i < ns; i++) {
        for (size_t p = 0; p < num_points; p++) {
            fprintf(gpd_file, "%e,%f", freq[i], points[p*3]);
            for (size_t k = 0; k < NRHS; k++) {
                index = (i * NRHS * num_points) + (p * NRHS) + k;
                var = ground_pot_s[index];
                fprintf(gpd_file, ",%e%+e%s", creal(var), cimag(var), "im");
            }
            fprintf(gpd_file, "\n");
        }
    }
    fclose(gpd_file);
    free(points);
    free(ground_pot_s);

    // Electric Fields  ===================================================
    char efield_file_name[60];
    sprintf(efield_file_name, "alipio83powsys_efield_%.0f.csv", 1 / sigma0);
    FILE* efield_file = fopen(efield_file_name, "w");
    fprintf(efield_file, "f,x");
    for (unsigned k = 0; k < NRHS; k++) {
        fprintf(efield_file, ",ecx%d,ecy%d,encx%d,ency%d", k+1, k+1, k+1, k+1);
    }
    fprintf(efield_file, "\n");
    for (size_t i = 0; i < ns; i++) {
        for (size_t p = 0; p < np_field; p++) {
            fprintf(efield_file, "%e,%f", freq[i], points[3 * p]);
            for (size_t k = 0; k < NRHS; k++) {
                for (size_t m = 0; m < 6; m++) {
                    if (m != 2 && m != 5) {  // exclude EZ
                        index = m + (k * 6) + (p * 6 * NRHS) + (i * 6 * NRHS * np_field);
                        var = efield_s[index];
                        fprintf(efield_file, ",%e%+e%s", creal(var), cimag(var), "im");
                    }
                }
            }
            fprintf(efield_file, "\n");
        }
    }
    fclose(efield_file);
    free(efield_s);

    // Voltage  ===================================================
    char voltage_file_name[60];
    sprintf(voltage_file_name, "alipio83powsys_voltage_%.0f.csv", 1 / sigma0);
    FILE* voltage_file = fopen(voltage_file_name, "w");
    fprintf(voltage_file, "f,x");
    for (unsigned k = 0; k < NRHS; k++) {
        fprintf(voltage_file, ",v%i", k+1);
    }
    fprintf(voltage_file, "\n");
    for (size_t i = 0; i < ns; i++) {
        for (size_t p = 0; p < np_field; p++) {
            fprintf(voltage_file, "%e,%f", freq[i], points[3 * p]);
            for (size_t k = 0; k < NRHS; k++) {
                index = k + (p * NRHS) + (i * NRHS * np_field);
                var = voltage_s[index];
                fprintf(voltage_file, ",%e%+e%s", creal(var), cimag(var), "im");
            }
            fprintf(voltage_file, "\n");
        }
    }
    fclose(voltage_file);
    free(voltage_s);
    return 0;
}

int
main (int argc, char *argv[])
{
    double start_time = omp_get_wtime();
    double lmax = 0.01;
    printf("============================\nCase 1\n");
    run_case(lmax, 1 / 100.0);
    printf("============================\nCase 2\n");
    run_case(lmax, 1 / 1000.0);
    double end_time = omp_get_wtime();
    printf("Elapsed time: %.2f minutes\n", (end_time - start_time) / 60.0);
    return 0;
}
