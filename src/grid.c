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
