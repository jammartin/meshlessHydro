//
// Created by Johannes Martin on 21.09.22.
//

#ifndef DEMONSTRATOR_PARAMETER_H
#define DEMONSTRATOR_PARAMETER_H

/// possible values 2 or 3 for 2D or 3D simulations
#define DIM 2

/// define if periodic boundaries should be employed
#define PERIODIC_BOUNDARIES 1

/// maximum number of interactions for each particle
#define MAX_NUM_INTERACTIONS 1000
/** maximum interactions with ghost particles
 *  ignored when `PERIODIC_BOUNDARIES` is not set
**/
#define MAX_NUM_GHOST_INTERACTIONS 750

#endif //DEMONSTRATOR_PARAMETER_H
