#ifndef Q_H
#define Q_H

#include "math.h"

inline float pair_force(float r, float Vo, float r1, float p1,
                  float r2, float p2, float d, float a);
inline float pair_energ(float r, float energy_cut, float Vo, float r1, float p1,
                  float r2, float p2, float d, float a);

float forces(float *x, long int* pairs, long int npairs, float Vo, float r1,
      float p1, float r2, float p2, float d, float a, float rcut, float *force);

#endif
