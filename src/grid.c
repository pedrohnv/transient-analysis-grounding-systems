#include <stdbool.h>
#include <float.h>
#include <complex.h>
#include <math.h>
#include "auxiliary.h"
#include "electrode.h"
#include "grid.h"

int
number_segments (const Grid grid)
{
    const unsigned int N = grid.edge_segments_x;
    const unsigned int vx = grid.vertices_x;
    const unsigned int M = grid.edge_segments_y;
    const unsigned int vy = grid.vertices_y;
    return ( N * vy * (vx - 1) + M * vx * (vy - 1) );
}

int
number_nodes (const Grid grid)
{
    const unsigned int N = grid.edge_segments_x;
    const unsigned int vx = grid.vertices_x;
    const unsigned int M = grid.edge_segments_y;
    const unsigned int vy = grid.vertices_y;
    return ( vx * vy + vx * (vy - 1) * (M - 1) + vy * (vx - 1) * (N - 1) );
}

int
electrode_grid (const Grid grid, Electrode* electrodes, double* nodes)
{
    const size_t N = grid.edge_segments_x;
    const double Lx = grid.length_x;
    const size_t vx = grid.vertices_x;
    const double lx = Lx / (N * (vx - 1));
    const size_t M = grid.edge_segments_y;
    const double Ly = grid.length_y;
    const size_t vy = grid.vertices_y;
    const double ly = Ly / (M * (vy - 1));
    double start_point[] = {0.0, 0.0, grid.depth};
    double end_point[] = {0.0, 0.0, grid.depth};
    //size_t Ns_hor = N*vy*(vx - 1);
    size_t ed = 0;
    size_t nd = 0;
    // Make horizontal electrodes
    for (size_t h = 0; h < vy; h++) {
        for (size_t n = 0; n < (vx - 1); n++) {
            for (size_t k = 0; k < N; k++) {
                //ed = (h*(vx - 1) + n) * N + k;
                start_point[0] = lx * (n * N + k);
                start_point[1] = ly * M * h;
                end_point[0] = start_point[0] + lx;
                end_point[1] = start_point[1];
                populate_electrode(&(electrodes[ed]), start_point, end_point,
                                   grid.radius);
                ed++;
                if (n == 0 && k == 0) {
                    nodes[3*nd + 0] = start_point[0];
                    nodes[3*nd + 1] = start_point[1];
                    nodes[3*nd + 2] = start_point[2];
                    nd++;
                }
                nodes[3*nd + 0] = end_point[0];
                nodes[3*nd + 1] = end_point[1];
                nodes[3*nd + 2] = end_point[2];
                nd++;
            }
        }
    }
    // Make vertical electrodes
    for (size_t g = 0; g < vx; g++) {
        for (size_t m = 0; m < (vy - 1); m++) {
            for (size_t k = 0; k < M; k++) {
                //ed = (g*(vy - 1) + m) * M + k + Ns_hor;
                start_point[0] = lx * N * g;
                start_point[1] = ly * (m * M + k);
                end_point[0] = start_point[0];
                end_point[1] = start_point[1] + ly;
                populate_electrode(&(electrodes[ed]), start_point, end_point,
                                   grid.radius);
                ed++;
                if (k < M - 1) {
                    nodes[3*nd + 0] = end_point[0];
                    nodes[3*nd + 1] = end_point[1];
                    nodes[3*nd + 2] = end_point[2];
                    nd++;
                }
            }
        }
    }
    return 0;
}


int
impedances_grid (_Complex double* zl, _Complex double* zt, const Grid grid,
                 _Complex double gamma, _Complex double s, double mur,
                 _Complex double kappa, size_t max_eval, double req_abs_error,
                 double req_rel_error, int integration_type, bool images)
{
    _Complex double iwu_4pi, one_4pik;
    if (integration_type == INTG_MHEM || integration_type == INTG_NONE) {
        iwu_4pi = 1.0;
        one_4pik = 1.0;
        gamma = 1.0;
        s = 1.0;
        mur = 1.0;
        kappa = 1.0;
    } else {
        iwu_4pi = s * mur * MU0 / (FOUR_PI);
        one_4pik = 1.0 / (FOUR_PI * kappa);
    }
    const size_t N = grid.edge_segments_x;
    const double Lx = grid.length_x;
    const unsigned int vx = grid.vertices_x;
    const double lx = Lx / (N * (vx - 1));
    const size_t M = grid.edge_segments_y;
    const double Ly = grid.length_y;
    const unsigned int vy = grid.vertices_y;
    const double ly = Ly / (M * (vy - 1));
    int square = ((fabs(lx - ly) < DBL_EPSILON) && (N == M) && (vx == vy));
    const double depth1 = grid.depth;
    const double depth2 = (images) ? -depth1 : depth1;
    const size_t num_hor = N * vy * (vx - 1);
    const size_t num_vert = M * vx * (vy - 1);
    const size_t num_seg = num_hor + num_vert;
    Electrode* sender = malloc(sizeof(Electrode));
    Electrode* receiver = malloc(sizeof(Electrode));
    double start_point[] = {0.0, 0.0, depth1};
    double end_point[] = {lx, 0.0, depth1};
    populate_electrode(sender, start_point, end_point, grid.radius);
    start_point[2] = depth2;
    end_point[2] = depth2;
    double x0, y0;
    double result[2], error[2];  // of the integration
    int failure = 0;
    _Complex double intg;
    size_t id1, id2, idx;
    // first SEGMENT to all horizontal others: Z(X[1,1,1]; X[h2,n2,k2])
    for (size_t h2 = 0; h2 < vy; h2++) {
        y0 = ly * M * h2;
        start_point[1] = y0;
        end_point[1] = y0;
        for (size_t n2 = 0; n2 < (vx - 1); n2++) {
            for (size_t k2 = 0; k2 < N; k2++) {
                x0 = lx * (N * n2 + k2);
                start_point[0] = x0;
                end_point[0] = x0 + lx;
                end_point[1] = y0;
                populate_electrode(receiver, start_point, end_point, grid.radius);
                id2 = (h2 * (vx - 1) + n2) * N + k2;
                if (k2 == 0 && n2 == 0 && h2 == 0 && !images) {
                    intg = self_integral(sender);
                } else {
                    failure = integral(sender, receiver, gamma, max_eval,
                                       req_abs_error, req_rel_error,
                                       integration_type, result, error);
                    if (failure) return failure;
                    intg = result[0] + I * result[1];
                }
                zt[id2] = intg;
                zt[num_seg * id2] = zt[id2];
            } // for k2
        } // for n2
    } // for h2
    // first edge to itself: Z(X[1,1,k1]; X[1,1,k2])
    for (size_t k1 = 1; k1 < N; k1++) {
        for (size_t k2 = k1; k2 < N; k2++) {
            zt[k2 + num_seg * k1] = zt[(k2 - 1) + num_seg * (k1 - 1)];
            zt[k1 + num_seg * k2] = zt[k2 + num_seg * k1];
        }
    }
    // first EDGE to all horizontal others: Z(X[1,1,k1]; X[h2,n2,k2])
    for (size_t k1 = 1; k1 < N; k1++) {
        for (size_t h2 = 0; h2 < vy; h2++) {
            for (size_t n2 = 0; n2 < (vx - 1); n2++) {
                for (size_t k2 = 0; k2 < N; k2++) {
                    id2 = (h2 * (vx - 1) + n2) * N + k2;
                    if (n2 == 0 && k2 == 0) {
                        id1 = (h2 * (vx - 1) + n2) * N + k1;
                        zt[id2 + num_seg * k1] = zt[id1];
                    } else {
                        zt[id2 + num_seg * k1] = zt[(id2 - 1) + num_seg * (k1 - 1)];
                    }
                    zt[k1 + num_seg * id2] = zt[id2 + num_seg * k1];
                } // for k2
            } // for n2
        } // for h2
    } // for k1
    // other horizontal to horizontal edges: Z(X[h1,n1,k1]; X[h2,n2,k2])
    for (size_t h1 = 0; h1 < vy; h1++) {
        for (size_t n1 = 0; n1 < (vx - 1); n1++) {
            id1 = (h1 * (vx - 1) + n1) * N;
            if (h1 > 0 || n1 > 0) {  // skip first edge
                for (size_t h2 = 0; h2 < vy; h2++) {
                    for (size_t n2 = 0; n2 < (vx - 1); n2++) {
                        id2 = (h2 * (vx - 1) + n2) * N;
                        if (h1 <= h2 && n1 <= n2) {
                            idx = ((h2 - h1) * (vx - 1) + n2 - n1) * N;
                            matrix_copy(zt + (idx),
                                        zt + (id2 + num_seg * id1),
                                        num_seg, num_seg, N, N);
                        } else if (h1 <= h2 && n1 > n2) {
                            idx = ((h2 - h1) * (vx - 1) + n1 - n2) * N;
                            transpose_copy(zt + (idx),
                                           zt + (id2 + num_seg * id1),
                                           num_seg, num_seg, N, N);
                        } else {
                            transpose_copy(zt + (id1 + num_seg * id2),
                                           zt + (id2 + num_seg * id1),
                                           num_seg, num_seg, N, N);
                        }
                    } // for n2
                } // for h2
            } // if
        } // for n1
    } // for h1
    // first EDGE to vertical ones: Z(X[1,1,k1]; Y[m1,g1,k2])
    bool c1, c2, c3, c4;
    for (size_t k1 = 0; k1 < N; k1++) {
        x0 = lx * k1;
        start_point[0] = x0;
        start_point[1] = 0.0;
        start_point[2] = depth1;
        end_point[0] = x0 + lx;
        end_point[1] = 0.0;
        end_point[2] = depth1;
        populate_electrode(sender, start_point, end_point, grid.radius);
        start_point[2] = depth2;
        end_point[2] = depth2;
        for (int g1 = 0; g1 < vx; g1++) {
            x0 = lx * N * g1;
            start_point[0] = x0;
            end_point[0] = x0;
            for (int m1 = 0; m1 < (vy - 1); m1++) {
                for (int k2 = 0; k2 < M; k2++) {
                    id2 = (g1 * (vy - 1) + m1) * M + k2 + num_hor;
                    c1 = (k1 > N / 2);
                    c2 = (k1 > N / 2 - 1) && (N % 2 == 0);
                    c3 = (g1 == 0);
                    c4 = (g1 == 1);
                    if (c3 && (c1 || c2)) {
                        idx = (vy - 1 + m1) * M + k2 + num_hor;
                        zt[id2 + k1 * num_seg] = zt[idx + num_seg * (N - 1 - k1)];
                    } else if (c4 && (c1 || c2)) {
                        idx = m1 * M + k2 + num_hor;
                        zt[id2 + k1 * num_seg] = zt[idx + num_seg * (N - 1 - k1)];
                    } else {
                        y0 = ly * (M * m1 + k2);
                        start_point[1] = y0;
                        end_point[1] = y0 + ly;
                        populate_electrode(receiver, start_point, end_point, grid.radius);
                        failure = integral(sender, receiver, gamma, max_eval,
                                           req_abs_error, req_rel_error,
                                           integration_type, result, error);
                        if (failure) return failure;
                        intg = result[0] + I * result[1];
                        zt[id2 + num_seg * k1] = intg;
                    } // if
                    zt[k1 + num_seg * id2] = zt[id2 + num_seg * k1];
                } // for k2
            } // for m1
        } // for g1
    } // for k1
    // other horizontal to vertical edges: Z(X[h1,n1,k1]; Y[m1,g1,k2])
    for (size_t h1 = 0; h1 < vy; h1++) {
        for (size_t n1 = 0; n1 < (vx - 1); n1++) {
            id1 = (h1 * (vx - 1) + n1) * N;
            if (h1 > 0 || n1 > 0) {  // skip first horizontal edge
                for (size_t g1 = 0; g1 < vx; g1++) {
                    for (size_t m1 = 0; m1 < (vy - 1); m1++) {
                        id2 = (g1 * (vy - 1) + m1) * M + num_hor;
                        if (h1 <= m1 && n1 <= g1) {
                            idx = ((g1 - n1) * (vy - 1) + (m1 - h1)) * M + num_hor;
                            matrix_copy(zt + (idx),
                                        zt + (id2 + num_seg * id1),
                                        num_seg, num_seg, M, N);
                        } else if (h1 <= m1 && n1 > g1) {
                            idx = ((n1 - g1 + 1)*(vy - 1) + (m1 - h1)) * M + num_hor;
                            pl_copy(zt + (idx),
                                    zt + (id2 + num_seg * id1),
                                    num_seg, num_seg, M, N);
                        } else if (h1 > m1 && n1 <= g1) {
                            idx = ((g1 - n1)*(vy - 1) + (h1 - m1 - 1)) * M + num_hor;
                            pc_copy(zt + (idx),
                                    zt + (id2 + num_seg * id1),
                                    num_seg, num_seg, M, N);
                        } else {
                            idx = ((n1 - g1 + 1)*(vy - 1) + (h1 - m1 - 1)) * M + num_hor;
                            pcl_copy(zt + (idx),
                                     zt + (id2 + num_seg * id1),
                                     num_seg, num_seg, M, N);
                        }
                        transpose_copy(zt + (id2 + num_seg * id1),
                                       zt + (id1 + num_seg * id2),
                                       num_seg, num_seg, M, N);
                    } // for m1
                } // for g1
            } // if
        } // for n1
    } // for h1
    if (square) {  // lx == ly && N == M && vx == vy
        matrix_copy(zt, zt + ((1 + num_seg) * num_hor),
                    num_seg, num_seg, num_hor, num_hor);
    } else {
        // first vertical SEGMENT to all vertical others: Z(Y[1,1,1]; Y[m2,g2,k2])
        start_point[0] = 0.0;
        start_point[1] = 0.0;
        start_point[2] = depth1;
        end_point[0] = 0.0;
        end_point[1] = ly;
        end_point[2] = depth1;
        populate_electrode(sender, start_point, end_point, grid.radius);
        start_point[2] = depth2;
        end_point[2] = depth2;
        id1 = num_hor;
        for (size_t g2 = 0; g2 < vx; g2++) {
            x0 = lx * N * g2;
            start_point[0] = x0;
            end_point[0] = x0;
            for (size_t m2 = 0; m2 < (vy - 1); m2++) {
                for (size_t k2 = 0; k2 < M; k2 ++) {
                    y0 = ly * (M * m2 + k2);
                    start_point[1] = y0;
                    end_point[1] = y0 + ly;
                    populate_electrode(receiver, start_point, end_point, grid.radius);
                    id2 = (g2 * (vy - 1) + m2) * M + k2 + num_hor;
                    failure = integral(sender, receiver, gamma, max_eval,
                                       req_abs_error, req_rel_error,
                                       integration_type, result, error);
                    if (failure) return failure;
                    intg = result[0] + I * result[1];
                    zt[id2 + num_seg * id1] = intg;
                    zt[id1 + num_seg * id2] = zt[id2 + num_seg * id1];
                } // for k2
            } // for m2
        } // for g2
        // first vertical edge to itself: Z(Y[1,1,k1]; Y[1,1,k2])
        for (size_t k1 = 1; k1 < M; k1++) {
            id1 = num_hor + k1;
            for (size_t k2 = k1; k2 < M; k2++) {
                id2 = num_hor + k2;
                zt[id2 + num_seg * id1] = zt[(id2 - 1) + num_seg * (id1 - 1)];
                zt[id1 + num_seg * id2] = zt[(id2 - 1) + num_seg * (id1 - 1)];
            }
        }
        // first vertical EDGE to all vertical others: Z(Y[1,1,k1]; Y[m2,g2,k2])
        for (size_t k1 = 1; k1 < M; k1++) {
            id1 = num_hor + k1;
            for (size_t g2 = 0; g2 < vx; g2++) {
                for (size_t m2 = 0; m2 < (vy - 1); m2++) {
                    for (size_t k2 = 0; k2 < M; k2++) {
                        id2 = (g2 * (vy - 1) + m2) * M + k2 + num_hor;
                        if (m2 == 0 && k2 == 0) {
                            idx = (g2 * (vy - 1) + m2) * M + k1 + num_hor;
                            zt[id2 + num_seg * id1] = zt[idx +  num_seg * num_hor];
                        } else {
                            zt[id2 + num_seg * id1] = zt[(id2 - 1) + num_seg * (id1 - 1)];
                        }
                        zt[id1 + num_seg * id2] = zt[id2 + num_seg * id1];
                    } // for k2
                } // for m2
            } // for g2
        } // for k1
        // other vertical to vertical edges: Z(Y[m1,g1,k1]; Y[m2,g2,k2])
        for (size_t g1 = 0; g1 < vx; g1++) {
            for (size_t m1 = 0; m1 < (vy - 1); m1++) {
                id1 = (g1 * (vy - 1) + m1) * M + num_hor;
                if (m1 > 0 || g1 > 0) {  // skip first vertical edge
                    for (size_t g2 = 0; g2 < vx; g2++) {
                        for (size_t m2 = 0; m2 < (vy - 1); m2++) {
                            id2 = (g2 * (vy - 1) + m2) * M + num_hor;
                            if (g1 <= g2 && m1 <= m2) {
                                idx = ((g2 - g1) * (vy - 1) + m2 - m1) * M + num_hor;
                                matrix_copy(zt + (idx + num_seg * num_hor),
                                            zt + (id2 + num_seg * id1),
                                            num_seg, num_seg, M, M);
                            } else if (g1 <= g2 && m1 > m2) {
                                idx = ((g2 - g1) * (vy - 1) + m1 - m2) * M + num_hor;
                                transpose_copy(zt + (idx + num_seg * num_hor),
                                               zt + (id2 + num_seg * id1),
                                               num_seg, num_seg, M, M);
                            } else {
                                transpose_copy(zt + (id1 + num_seg * id2),
                                               zt + (id2 + num_seg * id1),
                                               num_seg, num_seg, M, M);
                            }
                            transpose_copy(zt + (id2 + num_seg * id1),
                                           zt + (id1 + num_seg * id2),
                                           num_seg, num_seg, M, M);
                        } // for m2
                    } // for g2
                } // if
            } // for m1
        } // for g1
    } // if square *************************************
    free(sender);
    free(receiver);
    const _Complex double one_4pik_lx2 = one_4pik / (lx * lx);
    const _Complex double one_4pik_ly2 = one_4pik / (ly * ly);
    const _Complex double one_4pik_lxly = one_4pik / (lx * ly);
    for (size_t k = 0; k < num_hor; k++) {
        for (size_t i = 0; i < num_hor; i++) {
            zl[i + num_seg * k] = iwu_4pi * zt[k + num_seg * i];
            zt[i + num_seg * k] *= one_4pik_lx2;
        }
        for (size_t i = num_hor; i < num_seg; i++) {
            zl[i + num_seg * k] = 0.0;
            zt[i + num_seg * k] *= one_4pik_lxly;
        }
    }
    for (size_t k = num_hor; k < num_seg; k++) {
        for (size_t i = 0; i < num_hor; i++) {
            zl[i + num_seg * k] = 0.0;
            zt[i + num_seg * k] *= one_4pik_lxly;
        }
        for (size_t i = num_hor; i < num_seg; i++) {
            zl[i + num_seg * k] = iwu_4pi * zt[i + num_seg * k];
            zt[i + num_seg * k] *= one_4pik_ly2;
        }
    }
    return 0;
}
