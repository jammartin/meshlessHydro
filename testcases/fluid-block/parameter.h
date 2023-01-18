//
// Created by Johannes Martin on 21.09.22.
//

#ifndef DEMONSTRATOR_PARAMETER_H
#define DEMONSTRATOR_PARAMETER_H

/// possible values 2 or 3 for 2D or 3D simulations
#define DIM 3

/// define if periodic boundaries should be employed
#define PERIODIC_BOUNDARIES 0

/// maximum number of interactions for each particle
#define MAX_NUM_INTERACTIONS 400
/** maximum interactions with ghost particles
 *  ignored when `PERIODIC_BOUNDARIES` is not set
**/
#define MAX_NUM_GHOST_INTERACTIONS 300

/// flag for slope limiting, 0: no slope limiting
#define SLOPE_LIMITING 1

/// slope limiting parameter, ignored when `SLOPE_LIMITING` is false
#define BETA 4.

/// enforcing flux symmetry by only calculating on side
#define ENFORCE_FLUX_SYM 1

/// define if particles should move, otherwise a fixed grid is used
#define MOVE_PARTICLES 1

/** define debug level to enable additional output:
 * 0: no debug additions
 * 1: additional checks
 * 2: dump NNL and ghosts to files
**/
#define DEBUG_LVL 1

/// use first order quadrature point for Riemann problems
#define FIRST_ORDER_QUAD_POINT 1

/// define verbosity for each VERBOSITY_PARTICLES particles
// TODO: use this flag when debug level 1
#define VERBOSITY_PARTICLES 10

/// define how much tolerance of flux antisymmetry is allowed in checkFluxSymmetry
#define FLUX_SYM_TOL 1e-20

/// define how much difference is tolerated for I as rotation matrix in 3D for aligned vectors
//#define ROT_3D_ALIGN_TOL 1e-20

#endif //DEMONSTRATOR_PARAMETER_H
