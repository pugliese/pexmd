#ifndef P_H
#define P_H

#include "math.h"

inline float pair_force(float s2, float D, float qo2);
inline float pair_gorce(float s2, float D, float po2);
inline float pair_energ(float s2, float energy_cut, float D);

float forces(float *x, float *p, long int* pairs, long int npairs,
          float D, float qo, float po, float scut, float *force);

float gorces(float *x, float *p, long int* pairs, long int npairs,
          float D, float qo, float po, float scut, float *gorce);

float fgorces(float *x, float *p, long int* pairs, long int npairs,
          float D, float qo, float po, float scut, float *force, float *gorce);
#endif
