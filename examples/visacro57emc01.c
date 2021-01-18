/*
Reproduce Figure 4 in [1] of the experimental results of a
20 x 16 grounding grid for a slow and a fast current injected at the corner
and at the center of the grid. That is, there are four simulations in total.
In addition to calculating the GPR (Ground Potential Rise), the ground potential
distribution (GPD) and the electric field (both conservative and nonconservative)
along two horizontal lines are all calculated for each injection.
mHEM formulation is used.

To get the maximum performance, instead of using some of the high level routines
provided by HP-HEM, I opted to call the low level BLAS, LAPACK and FFTW3 routines
directly. This increased complexity in the code results in a very noticeable
gain in performance. The four simulations are done simultaneously, that is,
the same ZL and ZT matrices and the same factorized YN matrix are used to
solve the system for multiple Right Hand Sides (injected current vectors).
The code is parallelized with OpenMP. All arrays are allocated as being one-
dimensional. When talking about second or third dimension to them, that's an
abstraction.

[1] S. Visacro, R. Alipio, C. Pereira, M. Guimarães and M. A. O. Schroeder,
"Lightning Response of Grounding Grids: Simulated and Experimental Results,"
in IEEE Transactions on Electromagnetic Compatibility,
vol. 57, no. 1, pp. 121-127, Feb. 2015, doi: 10.1109/TEMC.2014.2362091
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include <fftw3.h>  // (fftw_complex) is (_Complex double) if fftw3.h is included after complex.h
#include <omp.h>
#include "auxiliary.h"
#include "electrode.h"
#include "linalg.h"
#include "cubature.h"
#include "grid.h"

/*
tmax : maximum simulation time [s]
nt : number of time steps
Lmax : segments maximum length [m]
inj_t : injected current vectors in an 1d array
NRHS : number of injections (simulations)
*/
int
run_case (double tmax, int nt, double Lmax, double* inj_t, const unsigned int NRHS)
{
    // numerical integration parameters
    size_t max_eval = 0;
    double req_abs_error = 1e-6;
    double req_rel_error = 1e-6;

    // soil model (Alipio) parameters
    double mur = 1.0;  // soil rel. magnetic permeability
    double sigma0 = 1.0 / 2000.0;  // soil conductivity in low frequency
    // parameters that I fitted:
    double h_soil = 0.925 * pow(sigma0 * 1e3, -0.73);
    double g_soil = 0.325;
    double eps_ratio = 1.850;  // soil rel. permitivitty ratio

    // Laplace transform of the injection currents =============================
    // do it explicitly for performance reasons, see FFTW3 manual: http://fftw.org/fftw3_doc/
    // its very unfortunate that FFTW3 relies on global variables created by itself
    // It is best at handling sizes of the form 2^a 3^b 5^c 7^d 11^e 13^f, where e+f = 0 or 1.
    // Transforms whose sizes are powers of 2 are especially fast.
    // see: http://www.fftw.org/fftw3_doc/Complex-DFTs.html#Complex-DFTs
    int ns = nt / 2 + 1;
    int rank = 1;  // we are computing 1d transforms
    int n[] = {nt};  // transforms of length N
    int howmany = NRHS;
    int idist = nt;  // the distance in memory  between the first element of the
                     // first array and the first element of the second array
    int odist = ns;
    int istride = 1;  // distance between two elements in the same array
    int ostride = 1;
    int* inembed = n;
    int* onembed = n;
    _Complex double* inj_s = malloc(howmany * ns * sizeof(_Complex double));
    _Complex double* s = malloc(ns * sizeof(_Complex double));
    double* nlt_input = malloc(howmany * nt * sizeof(double));
    _Complex double* nlt_output = malloc(howmany * ns * sizeof(_Complex double));
    // You must create the plan before initializing the input
    // Arrays n, inembed, and onembed are not used after this function returns.
    // You can safely free or reuse them.
    fftw_plan inj_plan = fftw_plan_many_dft_r2c(rank, n, howmany,
                                                nlt_input,  inembed, istride, idist,
                                                nlt_output, onembed, ostride, odist,
                                                FFTW_ESTIMATE);
    double c = log(pow(nt, 2.0)) / tmax;
    double dt = tmax / (nt - 1);
    double dw = TWO_PI / tmax;
    double var;
    for (int k = 0; k < ns; k++) {
        s[k] = c + I * dw * k;
    }
    for (int k = 0; k < nt; k++) {
        var = dt * exp(-c * k * dt);
        for (int i = 0; i < NRHS; i++) {
            nlt_input[k + nt * i] = var * inj_t[k + nt * i];
        }
    }
    fftw_execute(inj_plan);
    for (int k = 0; k < NRHS * ns; k++) {
        inj_s[k] = nlt_output[k];
    }

    // electrodes definition ===================================================
    double Lx = 20.0;
    double Ly = 16.0;
    const double h = -0.5;  // burial depth
    double length = 4.0;
    int div = length / Lmax;
    Grid grid = {6, 5, Lx, Ly, div, div, 5e-3, h};
    int ne = number_segments(grid);
    int nn = number_nodes(grid);
    printf("Num. segments = %i\n", ne);
    printf("Num. nodes    = %i\n", nn);
    Electrode* electrodes = malloc(sizeof(Electrode) * ne);
    double* nodes = malloc(nn * 3 * sizeof(double));
    electrode_grid(grid, electrodes, nodes);
    Electrode* images = malloc(sizeof(Electrode) * ne);
    for (size_t m = 0; m < ne; m++) {
        populate_electrode(images + m, electrodes[m].start_point,
                           electrodes[m].end_point, electrodes[m].radius);
        images[m].start_point[2] = -images[m].start_point[2];
        images[m].end_point[2] = -images[m].end_point[2];
        images[m].middle_point[2] = -images[m].middle_point[2];
    }
    // identify nodes ==========================================================
    double corner_point[] = {0.0, 0.0, h};
    size_t corner_node = 0;
    double tol = 1e-6;
    for (int i = 0; i <= nn; i++) {
        if (i == nn) {
            printf("Could not find corner node [0, 0]\n");
            exit(-1);
        }
        if (equal_points_tol(corner_point, nodes + 3*i, tol)) {
            corner_node = i;
            printf("corner node = %li\n", corner_node);
            break;
        }
    }
    double center_point[] = {8.0, 8.0, h};
    size_t center_node = 0;
    for (int i = 0; i <= nn; i++) {
        if (i == nn) {
            printf("Could not find center node [8, 8]\n");
            exit(-1);
        }
        if (equal_points_tol(center_point, nodes + 3*i, tol)) {
            center_node = i;
            printf("center node = %li\n", center_node);
            break;
        }
    }

    // define an array of points where to calculate ground scalar electric potential (GPD)
    // They form a grid of (1 x 1) [m^2] meshes
    double offset = 6.0;  // distance from groundig grid where to begin calculating GPD and fields
    Lx += offset * 2;
    Ly += offset * 2;
    Grid point_grid = {Lx + 1, Ly + 1, (int) Lx, (int) Ly, 1, 1, 1.0, 0.0};
    size_t num_points = number_nodes(point_grid);
    printf("Num. points to calculate GPD = %li\n", num_points);
    double* points = malloc(num_points * 3 * sizeof(double));
    Electrode* temp_elec = malloc(number_segments(point_grid) * sizeof(Electrode));
    electrode_grid(point_grid, temp_elec, points);
    free(temp_elec);
    // offset points to start at (-offset, -offset)
    for (int i = 0; i < num_points; i++) {
        points[3*i + 0] -= offset;
        points[3*i + 1] -= offset;
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

    // GPR =====================================================================
    // GPR[i,j] => GPR[i + j * ns]
    // dimensions: i-frequency, j-injection
    _Complex double* gpr_s = malloc(NRHS * ns * sizeof(_Complex double));
    double* gpr_t = malloc(NRHS * nt * sizeof(double));
    rank = 1;  // we are computing 1d transforms
    n[0] = nt;  // transforms of length N
    howmany = NRHS;
    idist = ns;  // the distance in memory  between the first element of the
                 // first array and the first element of the second array
    odist = nt;
    istride = 1;  // distance between two elements in the same vector
    ostride = 1;
    inembed = n;
    onembed = n;
    fftw_plan gpr_plan = fftw_plan_many_dft_c2r(rank, n, howmany,
                                                gpr_s, inembed, istride, idist,
                                                gpr_t, onembed, ostride, odist,
                                                FFTW_ESTIMATE);
    // GPD (Ground Potential Distribution) =====================================
    // GPD array: GPD[i,j,k] => GPD[(i * NRHS * num_points) + (p * NRHS) + k]
    // the outermost dimension iterates first
    // dimensions: i-frequency, j-injection, k-point
    // maybe the following helps with the visualization of the GPD array:
    //     potential at point 1 for injection 1 at frequency 1
    //     potential at point 2 for injection 1 at frequency 1
    //     ...
    //     potential at point 1 for injection 2 at frequency 1
    //     ...
    //     potential at point j for injection k at frequency i
    _Complex double* ground_pot_s = calloc(NRHS * num_points * ns, sizeof(_Complex double));
    double* ground_pot_t = malloc(NRHS * num_points * nt * sizeof(double));
    rank = 1;  // we are computing 1d transforms
    n[0] = nt;  // transforms of length N
    howmany = num_points * NRHS;
    idist = 1;  // the distance in memory between the first element of the
                // first array and the first element of the second array
    odist = 1;
    istride = NRHS * num_points;  // distance between two elements in the same vector
    ostride = NRHS * num_points;
    inembed = n;
    onembed = n;
    fftw_plan ground_pot_plan = fftw_plan_many_dft_c2r(rank, n, howmany,
                                            ground_pot_s, inembed, istride, idist,
                                            ground_pot_t, onembed, ostride, odist,
                                            FFTW_ESTIMATE);
    // Electric fields =========================================================
    // Along two lines: y = 0 and y = 8
    double efield_dr = 0.1;  // spatial step
    size_t np_field = (ceil(Lx / efield_dr) + 1) * 2;  // Lx was already added by 2*offset
    printf("Num. points to electric fields = %li\n", np_field);
    _Complex double* efield_s = malloc(NRHS * 6 * np_field * ns * sizeof(_Complex double));
    double* efield_t = malloc(NRHS * 6 * np_field * nt * sizeof(double));
    // efield[i,j,k,m] => efield[m + (k * 6) + (j * 6 * NRHS) + (i * 6 * NRHS * np_field)]
    // dimensions: i-frequency, j-point, k-injection, m-field
    rank = 1;  // we are computing 1d transforms
    n[0] = nt;  // transforms of length N
    howmany = np_field * NRHS * 6;
    idist = 1;  // the distance in memory  between the first element of the
                // first array and the first element of the second array
    odist = 1;
    istride = np_field * NRHS * 6;  // distance between two elements in the same vector
    ostride = np_field * NRHS * 6;
    inembed = n;
    onembed = n;
    fftw_plan efield_plan = fftw_plan_many_dft_c2r(rank, n, howmany,
                                              efield_s, inembed, istride, idist,
                                              efield_t, onembed, ostride, odist,
                                              FFTW_ESTIMATE);

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
        _Complex double* il = calloc(NRHS * ne, sizeof(_Complex double));
        _Complex double* it = calloc(NRHS * ne, sizeof(_Complex double));
        _Complex double* ie = malloc(NRHS * nn * sizeof(_Complex double));
        if (!zl || !zt || !ie || !yn || !yla || !ytb || !il || !it || !ie) {
            printf("Can't allocate memory\n");
            exit(200);
        }
        double field_point[3] = {0.0, 0.0, 0.0};
        _Complex double sigma, epsr, kappa, gamma;  // soil parameters
        _Complex double ref_l, ref_t;  // reflection coefficients (images)
        _Complex double iwu_4pi, one_4pik, exp_gr;
        double rbar, r0, r1, r2;
        //_Complex double gpd_var[NRHS];
        _Complex double efvector[3];  // helper variable of a 3d vector
        size_t index;

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
            alipio_soil(&sigma, &epsr, sigma0, s[i], h_soil, g_soil, eps_ratio);
            kappa = (sigma + s[i] * epsr * EPS0);  // soil complex conductivity
            gamma = csqrt(s[i] * MU0 * kappa);  // soil propagation constant
            iwu_4pi = s[i] * mur * MU0 / (FOUR_PI);
            one_4pik = 1.0 / (FOUR_PI * kappa);
            // reflection coefficient, soil / air
            ref_t = (kappa - s[i] * EPS0) / (kappa + s[i] * EPS0);
            ref_l = 0.0;  // longitudinal current does not have image
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
            ie[corner_node]            = inj_s[i];           // first column of IE
            ie[center_node + ldie]     = inj_s[i + ns];      // second column of IE
            ie[corner_node + ldie * 2] = inj_s[i + ns * 2];  // third column of IE
            ie[center_node + ldie * 3] = inj_s[i + ns * 3];  // fourth column of IE
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
            gpr_s[i] = ie[corner_node];  // first column of IE
            if (NRHS > 1) gpr_s[i + ns] = ie[center_node + ldie];  // second column of IE
            if (NRHS > 2) gpr_s[i + ns * 2] = ie[corner_node + ldie * 2];  // third column of IE
            if (NRHS > 3) gpr_s[i + ns * 3] = ie[center_node + ldie * 3];  // fourth column of IE

            err = 1;  // used in zgemv_ to specify incx, incy
            for (int k = 0; k < NRHS; k++) {
                // IT = inv(ZT) * B * IE + 0 * IT
                zgemv_(&trans, &ne, &nn, &one, yla, &ne, (ie + nn1 * k), &err,
                       &zero, (il + ne * k), &err);
                // IL = inv(ZL) * A * IE + 0 * IL
                zgemv_(&trans, &ne, &nn, &one, ytb, &ne, (ie + nn1 * k), &err,
                       &zero, (it + ne * k), &err);
                // if using intel MKL, replace the above by:
                /*cblas_zgemv(CblasColMajor, CblasNoTrans, ne, nn, &one, yla, ne,
                            (ie + nn1 * k), 1, &zero, (il + ne * k), 1);
                cblas_zgemv(CblasColMajor, CblasNoTrans, ne, nn, &one, ytb, ne,
                            (ie + nn1 * k), 1, &zero, (it + ne * k), 1);*/
            }
            // images' effect as (1 + ref) because we're calculating on ground level
            for (int k = 0; k < NRHS; k++) {
                for (size_t m = 0; m < ne; m++) {
                    it[m + k * ne] *= (1.0 + ref_t);
                    il[m + k * ne] *= (1.0 + ref_l);
                }
            }
            // Ground Potential Distribution (GPD) =============================
            for (size_t p = 0; p < num_points; p++) {
                // multiply by (1 + ref_t) to add the images' effects because
                // we are calculating on ground level
                // Traditional HEM:
                /*gpd_var[k] = electric_potential((points + p*3),
                                                electrodes, ne, (it + ne * k),
                                                gamma, kappa, max_eval,
                                                req_abs_error,
                                                req_rel_error);*/
                for (int k = 0; k < NRHS; k++) {
                    index = (i * NRHS * num_points) + (p * NRHS) + k;
                    ground_pot_s[index] = 0.0;
                }
                for (size_t m = 0; m < ne; m++) {
                // calculate using a simplification to the numerical integrals (mHEM):
                    rbar = vector_length((points + p*3), electrodes[m].middle_point);
                    r1 = vector_length((points + p*3), electrodes[m].start_point);
                    r2 = vector_length((points + p*3), electrodes[m].end_point);
                    r0 = (r1 + r2 + electrodes[m].length) / (r1 + r2 - electrodes[m].length);
                    exp_gr = cexp(-gamma * rbar) * log(r0) / electrodes[m].length;
                    for (int k = 0; k < NRHS; k++) {
                        index = (i * NRHS * num_points) + (p * NRHS) + k;
                        ground_pot_s[index] += one_4pik * it[m + ne * k] * exp_gr;
                    }
                }
            }
            // Electric Field ==================================================
            for (size_t p = 0; p < np_field; p++) {
                if (p < np_field / 2) {
                    field_point[0] = efield_dr * p - offset;
                    field_point[1] = 0.0;
                } else {
                    field_point[0] = efield_dr * (p - np_field / 2) - offset;
                    field_point[1] = 8.0;
                }
                for (size_t k = 0; k < NRHS; k++) {
                    // conservative field
                    // pass "s = 0.0" so that the non-conservative part is ignored
                    for (size_t m = 0; m < 3; m++) efvector[m] = 0.0;
                    electric_field(field_point, electrodes, ne,
                                   (il + ne * k), (it + ne * k),
                                   gamma, 0.0, mur, kappa, max_eval,
                                   req_abs_error, req_rel_error, efvector);
                    for (size_t m = 0; m < 3; m++) {
                        index = m + (k * 6) + (p * 6 * NRHS) + (i * 6 * NRHS * np_field);
                        efield_s[index] = efvector[m];
                    }
                    // non-conservative field
                    for (size_t m = 0; m < 3; m++) efvector[m] = 0.0;
                    magnetic_potential(field_point, electrodes, ne,
                                      (il + ne * k), gamma, mur, max_eval,
                                      req_abs_error, req_rel_error, efvector);
                    for (size_t m = 3; m < 6; m++) {
                        index = m + (k * 6) + (p * 6 * NRHS) + (i * 6 * NRHS * np_field);
                        efield_s[index] = efvector[m - 3] * (-s[i]);
                    }
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
    printf("Frequency loop ended. Saving results.\n");
    // GPR  ==============================================================
    double inlt_scale;
    fftw_execute(gpr_plan);
    char gpr_file_name[60];
    sprintf(gpr_file_name, "visacro57emc01_gpr.csv");
    FILE* gpr_file = fopen(gpr_file_name, "w");
    fprintf(gpr_file, "t");
    for (unsigned k = 0; k < NRHS; k++) fprintf(gpr_file, ",i%d", k+1);
    fprintf(gpr_file, "\n");
    for (size_t i = 0; i < nt; i++) {
        inlt_scale = exp(c * i * dt) / (nt * dt);
        fprintf(gpr_file, "%e", dt * i);
        for (size_t k = 0; k < NRHS; k++) {
            fprintf(gpr_file, ",%e", gpr_t[i + nt * k] * inlt_scale);
        }
        fprintf(gpr_file, "\n");
    }
    fclose(gpr_file);
    // GPD  ==============================================================
    size_t index;
    begin = omp_get_wtime();
    fftw_execute(ground_pot_plan);
    printf("Time to do INLT of the GPD: %f s\n", (omp_get_wtime() - begin));
    char gpd_file_name[60];
    sprintf(gpd_file_name, "visacro57emc01_gpd.csv");
    FILE* gpd_file = fopen(gpd_file_name, "w");
    fprintf(gpd_file, "t,x,y");
    for (unsigned k = 0; k < NRHS; k++) fprintf(gpd_file, ",i%d", k+1);
    fprintf(gpd_file, "\n");
    begin = omp_get_wtime();
    for (size_t i = 0; i < nt; i++) {
        inlt_scale = exp(c * i * dt) / (nt * dt);
        for (size_t p = 0; p < num_points; p++) {
            fprintf(gpd_file, "%e,%f,%f", dt * i, points[p*3], points[p*3 + 1]);
            for (size_t k = 0; k < NRHS; k++) {
                index = (i * NRHS * num_points) + (p * NRHS) + k;
                var = ground_pot_t[index] * inlt_scale;
                fprintf(gpd_file, ",%e", var);
            }
            fprintf(gpd_file, "\n");
        }
    }
    printf("Time to save GPD: %f s\n", (omp_get_wtime() - begin));
    fclose(gpd_file);
    // Electric Fields  ===================================================
    fftw_execute(efield_plan);
    fflush(stdout);
    char efield_file_name[60];
    sprintf(efield_file_name, "visacro57emc01_efield.csv");
    FILE* efield_file = fopen(efield_file_name, "w");
    fprintf(efield_file, "t,x,y");
    for (unsigned k = 0; k < NRHS; k++) {
        fprintf(efield_file, ",ecx%d,encx%d", k+1, k+1);
    }
    fprintf(efield_file, "\n");
    double x, y;
    for (size_t i = 0; i < nt; i++) {
        inlt_scale = exp(c * i * dt) / (nt * dt);
        for (size_t p = 0; p < np_field; p++) {
            if (p < np_field / 2) {
                x = efield_dr * p - offset;
                y = 0.0;
            } else {
                x = efield_dr * (p - np_field / 2) - offset;
                y = 8.0;
            }
            fprintf(efield_file, "%e,%f,%f", dt * i, x, y);
            for (size_t k = 0; k < NRHS; k++) {
                for (size_t m = 0; m < 6; m++) {
                    if (m == 0 || m == 3) {  // only Ex
                        index = m + (k * 6) + (p * 6 * NRHS) + (i * 6 * NRHS * np_field);
                        var = efield_t[index] * inlt_scale;
                        fprintf(efield_file, ",%e", var);
                    }
                }
            }
            fprintf(efield_file, "\n");
        }
    }
    fclose(efield_file);
    free(potzl);
    free(potzt);
    free(potzli);
    free(potzti);
    free(a);
    free(b);
    free(electrodes);
    free(images);
    free(nodes);
    free(nlt_input);
    free(nlt_output);
    free(gpr_s);
    free(gpr_t);
    free(ground_pot_s);
    free(ground_pot_t);
    free(efield_s);
    free(efield_t);
    free(points);
    fftw_destroy_plan(inj_plan);
    fftw_destroy_plan(gpr_plan);
    fftw_destroy_plan(ground_pot_plan);
    fftw_destroy_plan(efield_plan);
    return 0;
}

int
main (int argc, char *argv[])
{
    const unsigned int NRHS = 4;  // number of injections (simulations)
    if (argc != 2 && argc != 3) {
        printf("Wrong number of arguments, the following are needed:\n");
        printf("  L_max : segments maximum length\n");
        printf("  Nt (optional) : number of time steps\n");
        printf("    Nt defaults to 3001 and can be as high as 9000\n");
        printf("    the simulation is done with a time step of 20 ns.\n");
        exit(argc);
    }
    double dt = 20e-9;
    char *p;
    double Lmax = strtod(argv[1], &p);
    int nt = 3001;
    if (argc == 3) {
        nt = strtol(argv[2], &p, 10);
    }
    double tmax = dt * (nt - 1);
    printf("  final time [s] = %g\n", tmax);
    printf("  number of time steps = %i\n", nt);
    printf("  seg. max. length [m] = %g\n", Lmax);
    double start_time = omp_get_wtime();
    // open injected current files
    char fname_corner_fast[] = "examples/visacro57emc01_auxfiles/it_corner_fast.csv";
    FILE* file_corner_fast = fopen(fname_corner_fast, "r");
    if (file_corner_fast == NULL) {
        printf("Cannot open file %s\n", fname_corner_fast);
        return -11;
    }
    char fname_center_fast[] = "examples/visacro57emc01_auxfiles/it_center_fast.csv";
    FILE* file_center_fast = fopen(fname_center_fast, "r");
    if (file_center_fast == NULL) {
        printf("Cannot open file %s\n", fname_center_fast);
        return -12;
    }
    char fname_corner_slow[] = "examples/visacro57emc01_auxfiles/it_corner_slow.csv";
    FILE* file_corner_slow = fopen(fname_corner_slow, "r");
    if (file_corner_slow == NULL) {
        printf("Cannot open file %s\n", fname_corner_slow);
        return -13;
    }
    char fname_center_slow[] = "examples/visacro57emc01_auxfiles/it_center_slow.csv";
    FILE* file_center_slow = fopen(fname_center_slow, "r");
    if (file_center_slow == NULL) {
        printf("Cannot open file %s\n", fname_center_slow);
        return -14;
    }

    double* inj_t = malloc(NRHS * nt * sizeof(double));
    for (int i = 0; i < nt; i++) {
        fscanf(file_corner_fast, "%lf", inj_t + i);
        fscanf(file_center_fast, "%lf", inj_t + i + nt);
        fscanf(file_corner_slow, "%lf", inj_t + i + nt * 2);
        fscanf(file_center_slow, "%lf", inj_t + i + nt * 3);
    }
    fclose(file_corner_fast);
    fclose(file_center_fast);
    fclose(file_corner_slow);
    fclose(file_center_slow);
    run_case(tmax, nt, Lmax, inj_t, NRHS);
    free(inj_t);
    double end_time = omp_get_wtime();
    printf("Elapsed time: %.2f minutes\n", (end_time - start_time) / 60.0);
    return 0;
}
