#ifndef P_H
#define P_H

#include "math.h"
float pressure_lj_PBC(float *x, long int* pairs, long int npairs, float eps,
             float sigma, float rcut, float *force, float L);

float pair_force(float s2, float D, float qo2);
float pair_gorce(float s2, float D, float po2);
float pressure_pauli_PBC(float *x, float *p, long int* pairs, long int npairs,
             float D, float qo, float po, float scut, float *force, float *gorce, float L);
 #endif
