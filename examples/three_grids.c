/*
Simulates 3 grounding grids from communication towers in close proximity to
each other. They are considered either connected and isolated from each other
and two distinct soil conductivity.
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
connect_grids : connect the grids?
*/
int
time_domain (double tmax, int nt, double Lmax, double* inj_t, double sigma0, bool connect_grids)
{
    const unsigned int NRHS = 1;
    // numerical integration parameters
    size_t max_eval = 0;
    double req_abs_error = 1e-6;
    double req_rel_error = 1e-6;

    // soil model (Alipio) parameters
    double mur = 1.0;  // soil rel. magnetic permeability
    double h_soil = 1.26 * pow(sigma0 * 1e3, -0.73);
    double g_soil = 0.54;
    double eps_ratio = 12;  // soil rel. permitivitty ratio

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
    // Three 12x12 square grids. Make one, then copy and offset the x values.
    // each grid is 6 m from the other
    const int ngrids = 3;
    double grid_offset = 6.0;
    double Lx = 12.0;
    double Ly = 12.0;
    const double h = -0.5;  // burial depth;
    const double radius = 5e-3;
    int div = ceil(Lx / Lmax);
    if (div % 2 != 0) div += 1; // make sure div is even
    Grid grid = {2, 2, Lx, Ly, div, div, radius, h};
    int nseg_1grid = number_segments(grid);
    int nnode_1grid = number_nodes(grid);
    int ne = nseg_1grid * ngrids;
    int nn = nnode_1grid * ngrids;
    Electrode* electrodes = malloc(sizeof(Electrode) * ne);
    double* nodes = malloc(sizeof(double) * nn * 3);
    electrode_grid(grid, electrodes, nodes);  // get grid electrodes and nodes
    // copy and offset other grids
    double start_point[3], end_point[3];
    double dx;
    size_t index;
    for (int k = 1; k < ngrids; k++) {
        for (int m = 0; m < nseg_1grid; m++) {
            for (int i = 0; i < 3; i++) {
                if (i == 0) {
                    dx = (Lx + grid_offset) * k;
                } else {
                    dx = 0.0;
                }
                start_point[i] = electrodes[m].start_point[i] + dx;
                end_point[i] = electrodes[m].end_point[i] + dx;
            }
            index = m + nseg_1grid * k;
            populate_electrode(electrodes + index, start_point,
                               end_point, electrodes[m].radius);
        }
    }
    // copy and offset remaining nodes
    for (int k = 1; k < ngrids; k++) {
        for (int m = 0; m < nnode_1grid; m++) {
            index = (k * nnode_1grid + m) * 3;
            dx = (Lx + grid_offset) * k;
            nodes[index + 0] = nodes[m*3 + 0] + dx;  // x
            nodes[index + 1] = nodes[m*3 + 1];  // y
            nodes[index + 2] = nodes[m*3 + 2];  // z
        }
    }
    // identify nodes ==========================================================
    double inj_point[] = {0.0, 0.0, h};
    size_t inj_node = 0;
    double tol = 1e-6;
    for (int i = 0; i <= nn; i++) {
        if (i == nn) {
            printf("Could not find injection node [0, 0]\n");
            exit(-1);
        }
        if (equal_points_tol(inj_point, nodes + 3*i, tol)) {
            inj_node = i;
            printf("injection node = %li\n", inj_node);
            break;
        }
    }
    // connect grids ===========================================================
    if (connect_grids) {
        double connection_point[3];
        connection_point[1] = Ly / 2;
        connection_point[2] = h;
        size_t connection_nodes[(ngrids - 1) * 2];
        for (int k = 1; k < ngrids; k++) {
            for (int m = 0; m < 2; m++) {
                connection_point[0] = Lx * k + grid_offset * (k - 1 + m);
                for (int i = 0; i <= nn; i++) {
                    if (i == nn) {
                        printf("Could not find connection node\n");
                        exit(-1);
                    }
                    if (equal_points_tol(connection_point, nodes + 3*i, tol)) {
                        connection_nodes[(k - 1) * 2 + m] = i;
                        break;
                    }
                }
            }
        }
        int nseg_hor = ceil(grid_offset / Lmax);
        const int ne_grids = ne;
        ne += nseg_hor * (ngrids - 1);
        electrodes = realloc(electrodes, sizeof(Electrode) * ne);
        // don't repeat the connection nodes
        const int nn_grids = nn;
        nn += (nseg_hor - 1) * (ngrids - 1);
        nodes = realloc(nodes, sizeof(double) * 3 * nn);
        if (electrodes == NULL || nodes == NULL) {
            puts("Error (re)allocating memory");
            exit(1);
        }
        double* temp_nodes = malloc(3 * (nseg_hor + 1) * sizeof(double));
        int node_index;
        for (int k = 1; k < ngrids; k++) {
            for (int m = 0; m < 3; m++) {
                node_index = connection_nodes[(k - 1) * 2 + 0];
                start_point[m] = nodes[node_index * 3 + m];
                node_index = connection_nodes[(k - 1) * 2 + 1];
                end_point[m] = nodes[node_index * 3 + m];
            }
            segment_electrode(electrodes + ne_grids + nseg_hor * (k - 1),
                              temp_nodes,
                              nseg_hor, start_point, end_point, radius);
            for (int m = 0; m < 3 * (nseg_hor - 1); m++) {
                node_index = nn_grids + (nseg_hor - 1) * (k - 1);
                nodes[node_index * 3 + m] = temp_nodes[m + 3];
            }
        }
        free(temp_nodes);
    }
    printf("Num. segments = %i\n", ne);
    printf("Num. nodes    = %i\n", nn);
    // create images
    Electrode* images = malloc(sizeof(Electrode) * ne);
    for (size_t m = 0; m < ne; m++) {
        populate_electrode(images + m, electrodes[m].start_point,
                           electrodes[m].end_point, electrodes[m].radius);
        images[m].start_point[2] = -images[m].start_point[2];
        images[m].end_point[2] = -images[m].end_point[2];
        images[m].middle_point[2] = -images[m].middle_point[2];
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
    // GPR[i,j] => GPR[i + j * nn]
    // dimensions: i-node, j-frequency
    _Complex double* gpr_s = malloc(NRHS * nn * ns * sizeof(_Complex double));
    double* gpr_t = malloc(NRHS * nn * nt * sizeof(double));
    rank = 1;  // we are computing 1d transforms
    n[0] = nt;  // transforms of length N
    howmany = NRHS * nn;
    idist = 1;  // the distance in memory  between the first element of the
                // first array and the first element of the second array
    odist = 1;
    istride = nn;  // distance between two elements in the same vector
    ostride = nn;
    inembed = n;
    onembed = n;
    fftw_plan gpr_plan = fftw_plan_many_dft_c2r(rank, n, howmany,
                                                gpr_s, inembed, istride, idist,
                                                gpr_t, onembed, ostride, odist,
                                                FFTW_ESTIMATE);
    // Electric fields =========================================================
    // along horizontal line
    const double efield_y = Ly / 2;
    double offset = grid_offset;  // distance from groundig grid where to begin calculating fields
    double efield_dr = 0.1;  // spatial step
    double efield_line_length = ngrids * Lx + 2 * offset;
    size_t np_field = (ceil((efield_line_length) / efield_dr) + 1);
    printf("Num. points to calc. electric fields = %li\n", np_field);
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
    #pragma omp parallel private(err, index)
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
        _Complex double* ie =  malloc(nn * sizeof(_Complex double));
        if (!zl || !zt || !yn || !yla || !ytb || !il || !it || !ie) {
            printf("Can't allocate memory\n");
            exit(200);
        }
        double field_point[3] = {0.0, 0.0, 0.0};
        _Complex double sigma, epsr, kappa, gamma;  // soil parameters
        _Complex double ref_l, ref_t;  // reflection coefficients (images)
        _Complex double iwu_4pi, one_4pik, exp_gr;
        double rbar;
        _Complex double efvector[3];  // helper variable of a 3d vector

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
            ref_l = 1.0;  // longitudinal current has +1 image
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
            ie[inj_node] = inj_s[i];
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
            for (size_t m = 0; m < (NRHS * nn); m++) {
                gpr_s[m + i * nn] = ie[m];
            }
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
            // Electric Field ==================================================
            for (size_t p = 0; p < np_field; p++) {
                field_point[0] = efield_dr * (double)p - offset;
                field_point[1] = efield_y;
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
        free(il);
        free(it);
        free(ie);
        free(ipiv);
        free(work);
    }  // end parallel
    printf("Frequency loop ended. Saving results.\n");
    const char* connect_string = (connect_grids) ? "_connected" : "_isolated";
    // GPR  ==============================================================
    double inlt_scale, x, y;
    fftw_execute(gpr_plan);
    char gpr_file_name[60];
    sprintf(gpr_file_name, "three_grids_time_gpr_%.0f%s.csv", 1 / sigma0, connect_string);
    FILE* gpr_file = fopen(gpr_file_name, "w");
    fprintf(gpr_file, "t,x,y");
    for (unsigned k = 0; k < NRHS; k++) fprintf(gpr_file, ",i%d", k+1);
    fprintf(gpr_file, "\n");
    for (size_t i = 0; i < nt; i++) {
        inlt_scale = exp(c * i * dt) / (nt * dt);
        for (size_t k = 0; k < nn; k++) {
            x = nodes[k * 3 + 0];
            y = nodes[k * 3 + 1];
            fprintf(gpr_file, "%e", dt * i);
            fprintf(gpr_file, ",%f,%f", x, y);
            fprintf(gpr_file, ",%e", gpr_t[k + i * nn] * inlt_scale);
            fprintf(gpr_file, "\n");
        }
    }
    fclose(gpr_file);
    // Electric Fields  ===================================================
    fftw_execute(efield_plan);
    fflush(stdout);
    char efield_file_name[60];
    sprintf(efield_file_name, "three_grids_time_efield_%.0f%s.csv", 1 / sigma0, connect_string);
    FILE* efield_file = fopen(efield_file_name, "w");
    fprintf(efield_file, "t,x,y");
    for (unsigned k = 0; k < NRHS; k++) {
        fprintf(efield_file, ",ecx%d,encx%d", k+1, k+1);
    }
    fprintf(efield_file, "\n");
    for (int i = 0; i < nt; i++) {
        inlt_scale = exp(c * i * dt) / (nt * dt);
        for (int p = 0; p < np_field; p++) {
            x = efield_dr * (double)p - offset;
            y = efield_y;
            fprintf(efield_file, "%e,%f,%f", dt * i, x, y);
            for (int k = 0; k < NRHS; k++) {
                for (int m = 0; m < 6; m++) {
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
    free(efield_s);
    free(efield_t);
    fftw_destroy_plan(inj_plan);
    fftw_destroy_plan(gpr_plan);
    fftw_destroy_plan(efield_plan);
    return 0;
}


/**
nf : number of frequencies of interest (logspaced from 10^2 to 10^7)
Lmax : maximum conductor's length
sigma0 : low frequency soil conductivity
connect_grids : connect the grids?
*/
int
harmonic_analysis (int nf, double Lmax, double sigma0, bool connect_grids)
{
    const unsigned int NRHS = 1;
    // numerical integration parameters
    size_t max_eval = 0;
    double req_abs_error = 1e-6;
    double req_rel_error = 1e-6;

    // soil model (Alipio) parameters
    double mur = 1.0;  // soil rel. magnetic permeability
    double h_soil = 1.26 * pow(sigma0 * 1e3, -0.73);
    double g_soil = 0.54;
    double eps_ratio = 12;  // soil rel. permitivitty ratio

    // frequencies
    unsigned ns = nf;
    double freq[ns];
    logspace(2.0, 7.0, ns, freq);
    _Complex double s[ns];
    for (unsigned i = 0; i < ns; i++) {
        s[i] = I * TWO_PI * freq[i];
    }
    // electrodes definition ===================================================
    // Three 12x12 square grids. Make one, then copy and offset the x values.
    // each grid is 6 m from the other
    const int ngrids = 3;
    double grid_offset = 6.0;
    double Lx = 12.0;
    double Ly = 12.0;
    const double h = -0.5;  // burial depth;
    const double radius = 5e-3;
    int div = ceil(Lx / Lmax);
    if (div % 2 != 0) div += 1; // make sure div is even
    Grid grid = {2, 2, Lx, Ly, div, div, radius, h};
    int nseg_1grid = number_segments(grid);
    int nnode_1grid = number_nodes(grid);
    int ne = nseg_1grid * ngrids;
    int nn = nnode_1grid * ngrids;
    Electrode* electrodes = malloc(sizeof(Electrode) * ne);
    double* nodes = malloc(sizeof(double) * nn * 3);
    electrode_grid(grid, electrodes, nodes);  // get grid electrodes and nodes
    // copy and offset other grids
    double start_point[3], end_point[3];
    double dx;
    size_t index;
    for (int k = 1; k < ngrids; k++) {
        for (int m = 0; m < nseg_1grid; m++) {
            for (int i = 0; i < 3; i++) {
                if (i == 0) {
                    dx = (Lx + grid_offset) * k;
                } else {
                    dx = 0.0;
                }
                start_point[i] = electrodes[m].start_point[i] + dx;
                end_point[i] = electrodes[m].end_point[i] + dx;
            }
            index = m + nseg_1grid * k;
            populate_electrode(electrodes + index, start_point,
                               end_point, electrodes[m].radius);
        }
    }
    // copy and offset remaining nodes
    for (int k = 1; k < ngrids; k++) {
        for (int m = 0; m < nnode_1grid; m++) {
            index = (k * nnode_1grid + m) * 3;
            dx = (Lx + grid_offset) * k;
            nodes[index + 0] = nodes[m*3 + 0] + dx;  // x
            nodes[index + 1] = nodes[m*3 + 1];  // y
            nodes[index + 2] = nodes[m*3 + 2];  // z
        }
    }
    // identify nodes ==========================================================
    double inj_point[] = {0.0, 0.0, h};
    size_t inj_node = 0;
    double tol = 1e-6;
    for (int i = 0; i <= nn; i++) {
        if (i == nn) {
            printf("Could not find injection node [0, 0]\n");
            exit(-1);
        }
        if (equal_points_tol(inj_point, nodes + 3*i, tol)) {
            inj_node = i;
            printf("injection node = %li\n", inj_node);
            break;
        }
    }
    // connect grids ===========================================================
    if (connect_grids) {
        double connection_point[3];
        connection_point[1] = Ly / 2;
        connection_point[2] = h;
        size_t connection_nodes[(ngrids - 1) * 2];
        for (int k = 1; k < ngrids; k++) {
            for (int m = 0; m < 2; m++) {
                connection_point[0] = Lx * k + grid_offset * (k - 1 + m);
                for (int i = 0; i <= nn; i++) {
                    if (i == nn) {
                        printf("Could not find connection node\n");
                        exit(-1);
                    }
                    if (equal_points_tol(connection_point, nodes + 3*i, tol)) {
                        connection_nodes[(k - 1) * 2 + m] = i;
                        break;
                    }
                }
            }
        }
        int nseg_hor = ceil(grid_offset / Lmax);
        const int ne_grids = ne;
        ne += nseg_hor * (ngrids - 1);
        electrodes = realloc(electrodes, sizeof(Electrode) * ne);
        // don't repeat the connection nodes
        const int nn_grids = nn;
        nn += (nseg_hor - 1) * (ngrids - 1);
        nodes = realloc(nodes, sizeof(double) * 3 * nn);
        if (electrodes == NULL || nodes == NULL) {
            puts("Error (re)allocating memory");
            exit(1);
        }
        double* temp_nodes = malloc(3 * (nseg_hor + 1) * sizeof(double));
        int node_index;
        for (int k = 1; k < ngrids; k++) {
            for (int m = 0; m < 3; m++) {
                node_index = connection_nodes[(k - 1) * 2 + 0];
                start_point[m] = nodes[node_index * 3 + m];
                node_index = connection_nodes[(k - 1) * 2 + 1];
                end_point[m] = nodes[node_index * 3 + m];
            }
            segment_electrode(electrodes + ne_grids + nseg_hor * (k - 1),
                              temp_nodes,
                              nseg_hor, start_point, end_point, radius);
            for (int m = 0; m < 3 * (nseg_hor - 1); m++) {
                node_index = nn_grids + (nseg_hor - 1) * (k - 1);
                nodes[node_index * 3 + m] = temp_nodes[m + 3];
            }
        }
        free(temp_nodes);
    }
    printf("Num. segments = %i\n", ne);
    printf("Num. nodes    = %i\n", nn);
    // create images ===========================================================
    Electrode* images = malloc(sizeof(Electrode) * ne);
    for (size_t m = 0; m < ne; m++) {
        populate_electrode(images + m, electrodes[m].start_point,
                           electrodes[m].end_point, electrodes[m].radius);
        images[m].start_point[2] = -images[m].start_point[2];
        images[m].end_point[2] = -images[m].end_point[2];
        images[m].middle_point[2] = -images[m].middle_point[2];
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
    _Complex double* gpr = calloc(NRHS * ns * nn, sizeof(_Complex double));
    // Electric fields =========================================================
    // along horizontal line
    const double efield_y = Ly / 2;
    double offset = grid_offset;  // distance from groundig grid where to begin calculating fields
    double efield_dr = 0.1;  // spatial step
    double efield_line_length = ngrids * Lx + 2 * offset;
    size_t np_field = (ceil((efield_line_length) / efield_dr) + 1);
    printf("Num. points to calc. electric fields = %li\n", np_field);
    _Complex double* efield_s = malloc(NRHS * 6 * np_field * ns * sizeof(_Complex double));
    // efield[i,j,k,m] => efield[m + (k * 6) + (j * 6 * NRHS) + (i * 6 * NRHS * np_field)]
    // dimensions: i-frequency, j-point, k-injection, m-field
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
    #pragma omp parallel private(err, index)
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
        _Complex double* ie = malloc((ne * nn) * sizeof(_Complex double));
        if (!zl || !zt ||  !yn || !yla || !ytb || !il || !it) {
            printf("Can't allocate memory\n");
            exit(200);
        }
        double field_point[3] = {0.0, 0.0, 0.0};
        _Complex double sigma, epsr, kappa, gamma;  // soil parameters
        _Complex double ref_l, ref_t;  // reflection coefficients (images)
        _Complex double iwu_4pi, one_4pik, exp_gr;
        double rbar;
        //_Complex double gpd_var[NRHS];
        _Complex double efvector[3];  // helper variable of a 3d vector

        // use raw BLAS and LAPACK to solve system with multiple RHS
        int nn1 = (int) nn;
        int ldie = nn1;  // leading dimension of ie
        int* ipiv = malloc(nn1 * sizeof(int));  // pivot indices
        int info;
        int nrhs = NRHS;
        char uplo = 'L';  // matrices are symmetric, only Lower Half is set
        int lwork = -1;  // signal to Query the optimal workspace
        _Complex double wkopt;
        zsysv_(&uplo, &nn1, &nrhs, yn, &nn1, ipiv, gpr, &ldie, &wkopt, &lwork, &info);
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
            ref_l = 1.0;  // longitudinal current has +1 image
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
            ie[inj_node] = 1.0;
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
            for (size_t m = 0; m < (NRHS * nn); m++) {
                gpr[m + i * nn] = ie[m];
            }
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
            // Electric Field ==================================================
            for (size_t p = 0; p < np_field; p++) {
                field_point[0] = efield_dr * p - offset;
                field_point[1] = efield_y;
                field_point[2] = 0.0;
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
        free(il);
        free(it);
        free(ie);
        free(ipiv);
        free(work);
    }  // end parallel
    printf("Frequency loop ended. Saving results.\n");
    // GPR  ==============================================================
    double x, y;
    _Complex double var;
    const char* connect_string = (connect_grids) ? "_connected" : "_isolated";
    char gpr_file_name[60];
    sprintf(gpr_file_name, "three_grids_freq_gpr_%.0f%s.csv", 1 / sigma0, connect_string);
    FILE* gpr_file = fopen(gpr_file_name, "w");
    fprintf(gpr_file, "f,x,y");
    for (unsigned k = 0; k < NRHS; k++) fprintf(gpr_file, ",i%d", k+1);
    fprintf(gpr_file, "\n");
    for (size_t i = 0; i < ns; i++) {
        for (size_t k = 0; k < nn; k++) {
            x = nodes[k * 3 + 0];
            y = nodes[k * 3 + 1];
            var = gpr[i * nn + k];
            fprintf(gpr_file, "%e", freq[i]);
            fprintf(gpr_file, ",%f,%f", x, y);
            fprintf(gpr_file, ",%e%+e%s", creal(var), cimag(var), "im");
            fprintf(gpr_file, "\n");
        }
    }
    fclose(gpr_file);
    // Electric Fields  ===================================================
    fflush(stdout);
    char efield_file_name[60];
    sprintf(efield_file_name, "three_grids_freq_efield_%.0f%s.csv", 1 / sigma0, connect_string);
    FILE* efield_file = fopen(efield_file_name, "w");
    fprintf(efield_file, "f,x,y");
    for (unsigned k = 0; k < NRHS; k++) {
        fprintf(efield_file, ",ecx%d,encx%d", k+1, k+1);
    }
    fprintf(efield_file, "\n");
    for (size_t i = 0; i < ns; i++) {
        for (size_t p = 0; p < np_field; p++) {
            x = efield_dr * p - offset;
            y = efield_y;
            fprintf(efield_file, "%e,%f,%f", freq[i], x, y);
            for (size_t k = 0; k < NRHS; k++) {
                for (size_t m = 0; m < 6; m++) {
                    if (m == 0 || m == 3) {  // only Ex
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
    free(potzl);
    free(potzt);
    free(potzli);
    free(potzti);
    free(a);
    free(b);
    free(electrodes);
    free(images);
    free(nodes);
    free(gpr);
    free(efield_s);
    return 0;
}

int
main (int argc, char *argv[])
{
    if (argc != 5) {
        printf("Wrong number of arguments, the following are needed:\n");
        printf("  L_max : segments maximum length\n");
        printf("  Nf : Number of frequency steps for harmonic analysis\n");
        printf("  Nt : number of time steps\n");
        printf("  tmax : final simulation time\n");
        exit(argc);
    }
    char *p;
    double Lmax = strtod(argv[1], &p);
    int nf = strtol(argv[2], &p, 10);
    int nt = strtol(argv[3], &p, 10);
    double tmax = strtod(argv[4], &p);
    double dt = tmax / (nt - 1);
    printf("Time step dt = %e\n", dt);
    double start_time = omp_get_wtime();
    double* inj_t = malloc(nt * sizeof(double));
    double alpha = 1.0 / 50e-6;
    double beta = 1.0 / 0.1e-6;
    for (int i = 0; i < nt; i++) {
        inj_t[i] = exp(-alpha * i * dt) - exp(-beta * i * dt);
    }
    printf("Harmonic Analysis =========================\n");
    printf("Disconnected grids ====\n");
    printf("  rho0 = 100\n");
    harmonic_analysis(nf, Lmax,  1 / 100.0, false);
    printf("  rho0 = 1000\n");
    harmonic_analysis(nf, Lmax,  1 / 1000.0, false);
    printf("Connected grids ===\n");
    printf("  rho0 = 100\n");
    harmonic_analysis(nf, Lmax,  1 / 100.0, true);
    printf("  rho0 = 1000\n");
    harmonic_analysis(nf, Lmax,  1 / 1000.0, true);
    printf("Time domain ===============================\n");
    printf("Disconnected grids ====\n");
    printf("  rho0 = 100\n");
    time_domain(tmax, nt, Lmax, inj_t, 1 / 100.0, false);
    printf("  rho0 = 1000\n");
    time_domain(tmax, nt, Lmax, inj_t, 1 / 1000.0, false);
    printf("Connected grids ===\n");
    printf("  rho0 = 100\n");
    time_domain(tmax, nt, Lmax, inj_t, 1 / 100.0, true);
    printf("  rho0 = 1000\n");
    time_domain(tmax, nt, Lmax, inj_t, 1 / 1000.0, true);
    free(inj_t);
    double end_time = omp_get_wtime();
    printf("Elapsed time: %.2f minutes\n", (end_time - start_time) / 60.0);
    return 0;
}
