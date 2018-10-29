#ifndef C_H
#define C_H

#include "math.h"

inline float pair_force(float r, float K);
inline float pair_energ(float r, float energy_cut, float K);

float forces(float *x, long int* pairs, long int npairs,
              float K, float rcut, float *force);

#endif
