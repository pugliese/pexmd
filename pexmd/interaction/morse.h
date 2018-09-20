#ifndef MF_H
#define MF_H

#include "math.h"

inline float pair_force(float r, float req, float D, float alpha);
inline float pair_energ(float r, float energy_cut, float req, float D, float alpha);

float forces(float *x, long int* pairs, long int npairs,
              float alpha,float D, float req, float rcut, float *force);

#endif
