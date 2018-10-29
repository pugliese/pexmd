#ifndef NN_H
#define NN_H

#include "math.h"

inline float pair_force(float r, float mu, float D);
inline float pair_energ(float r, float energy_cut, float mu, float D);

float forces(float *x, long int* pairs, long int npairs,
              float D, float mu, float rcut, float *force);

#endif
