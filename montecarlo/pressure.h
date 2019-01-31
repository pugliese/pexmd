#ifndef P_H
#define P_H

#include "math.h"
float pressure_lj_PBC(float *x, long int* pairs, long int npairs, float eps,
             float sigma, float rcut, float *force, float L);

float pressure_qcnm_layers(float *x, long int* pairs, long int npairs, float Vo, float r1,
     float p1, float r2, float p2, float d, float a, float rcut, float *force, float L, int ls);
float pair_force_qcnm(float r, float Vo, float r1, float p1, float r2, float p2, float d, float a);

float pair_force(float s2, float D, float qo2);
float pair_gorce(float s2, float D, float po2);
float pressure_pauli_PBC(float *x, float *p, long int* pairs, long int npairs,
             float D, float qo, float po, float scut, float *force, float *gorce, float L);
float delta_fases(float *x, float *p, long int* pairs, long int npairs, float *dq, float *dp, float L);
float delta_fases_sin_PBC(float *x, float *p, long int* pairs, long int npairs, float *dq, float *dp);
 #endif
