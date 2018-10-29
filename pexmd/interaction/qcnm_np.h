#ifndef NP_H
#define NP_H

#include "math.h"

float pair_force(float r, float mur, float Dr,
              float mua, float Da);
float pair_energ(float r, float energy_cut, float mur,
              float Dr, float mua, float Da);

float forces(float *x, long int* pairs, long int npairs, float Dr,
          float mur, float Da, float mua, float rcut, float *force);

#endif
