/* High Performance implementation of the Hybrid Electromagnetic Model
Released under the General Public License 3 (GPLv3).
All parameters' units are in the SI base if omitted.

Routines to manipulate the Grid (of Electrodes) structure and do calculations,
i.e., numerical integration and impedances calculation for this specialized case.
*/
#ifndef GRID_H_
#define GRID_H_

#include <complex.h>
#include <stdlib.h>
#include <stdbool.h>
#include "electrode.h"

/** Strutcture to represent an uniform rectangular grid to be used in specialized
routines. This grid has dimensions \f$ L_x \, L_y \f$, a total of (before segmentation)
    \f$ nv = (v_x \, v_y) \f$
vertices and
    \f$ ne = v_y\, (v_x - 1) + v_x\, (v_y - 1) \f$
edges. Each edge is divided into \f$ N \f$ segments so that the total number of nodes
after segmentation is
    \f$ n = v_x\,v_y + v_x\,(v_y - 1)(N_y - 1) + v_y\,(v_x - 1)(N_x - 1) \f$
and the total number of segments is
    \f$ m = N_x\,v_x\,(v_y - 1) + N_y\,v_y\,(v_x - 1) \f$
\n \verbatim
     1          vx
 -  o---o---o---o 1
 |  |   |   |   |
Ly  o---o---o---o
 |  |   |   |   |
 -  o---o---o---o vy
    |---- Lx ---| \endverbatim
*/
typedef struct {
    /** number of vertices \f$v_x\f$ in the \f$\vec x\f$ direction */
    unsigned int vertices_x;
    /** number of vertices \f$v_y\f$ in the \f$\vec y\f$ direction */
    unsigned int vertices_y;
    /** total grid length \f$L_x\f$ in the \f$\vec x\f$ direction; */
    double length_x;
    /** total grid length \f$L_y\f$ in the \f$\vec y\f$ direction; */
    double length_y;
    /** number of segments \f$N_x\f$ that each edge has in the \f$\vec x\f$ direction */
    unsigned int edge_segments_x;
    /** number of segments \f$N_y\f$ that each edge has in the \f$\vec y\f$ direction */
    unsigned int edge_segments_y;
    /** Grid conductors' radius (the same for all segments) */
    double radius;
    /** Grid coordinate in the \f$\vec z\f$ direction (the same for all segments) */
    double depth;
} Grid;

/** Returns how many segments the Grid has. */
int
number_segments (const Grid grid);

/** Returns how many nodes the grid has. */
int
number_nodes (const Grid grid);

/** Populates electrodes and nodes arrays given a Grid structure.
@param grid structure
@param electrodes array to be filled
@param nodes array to be filled
@returns 0 on success
@see number_segments for the minimum size of the array electrodes
@see number_nodes for the minimum size of the array nodes
*/
int
electrode_grid (const Grid grid, Electrode* electrodes, double* nodes);

/** Calculates the impedance matrices for a regular grid, based on [1].
To add the effect of images, call this function separately with 'images'
parameter as true. It will be considered that the medium interface is at \f$ z = 0 \f$.

[1] VIEIRA, Pedro H. N. et al. Symmetry exploitation to reduce impedance
evaluations in grounding grids. International Journal of Electrical Power &
Energy Systems, v. 123, p. 106268, 2020.

@param grid structure
@param zl longitudinal impedance matrix \f$ Z_L \f$ as a flat array of size \f$ m^2 \f$
@param zt transversal impedance matrix \f$ Z_T \f$ as a flat array of size \f$ m^2 \f$
@param gamma medium propagation constant
    \f$ \gamma = \sqrt{j\omega\mu(\sigma + j\omega\varepsilon)} \f$
@param s complex angular frequency \f$ c + j\omega \f$ in [rad/s]
@param mur relative magnetic permeability of the medium \f$ \mu_r \f$
@param kappa medium complex conductivity \f$ \sigma + j\omega\varepsilon \f$ in [S/m]
@param max_eval specifies a maximum number of function evaluations (0 for no limit)
@param req_abs_error the absolute error requested (0 to ignore)
@param req_rel_error the relative error requested (0 to ignore)
@param integration_type type of integration to be done.
@param images calculating impedances with images?
@return 0 on success
@see integral
@see Integration_type
@see https://github.com/stevengj/cubature
*/
int
impedances_grid (_Complex double* zl, _Complex double* zt, const Grid grid,
                 _Complex double gamma, _Complex double s, double mur,
                 _Complex double kappa, size_t max_eval, double req_abs_error,
                 double req_rel_error, int integration_type, bool images);

#endif /* GRID_H_ */
