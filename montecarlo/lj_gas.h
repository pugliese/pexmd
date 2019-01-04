#ifndef LJG_H
#define LJG_H

#include "math.h"

struct Particles {
  int n;
  float* q;
  float* p;
  float mass;
  float energy;
  float kinetic;
};

struct LJ {
  float rcut2;
  float shift;
};

struct Externos {
  float L;
  float T;
  float delta_q;
  float delta_p;
};

float sub_vector(float* V, int i, float* vec);
int elegir_part(int N);
float uniform();
float boltzmann(float sigma);
float distancia(float* q1, float* q2, struct LJ *lj, float L);
float interaction(float *q1, float *q2, struct LJ *lj, float L);
int energia(struct Particles *parts, struct LJ *lj, float L);
float delta_energia_kin(struct Particles *parts, float *new_p, int i);
float delta_energia_pot(struct Particles *parts, struct LJ *lj, float L, float *new_q, float *new_p, int i);
float set_box(struct Particles *parts, float L);
float set_p(struct Particles *parts, float T);
int step(struct Particles *parts, struct LJ *lj, struct Externos *params);
int N_steps(struct Particles *parts, struct LJ *lj, struct Externos *params, int Nsteps);
int muestrear_impulsos(char *filename, struct Particles *parts, struct LJ *lj, struct Externos *params, int Nsamp, int factor, int factor_term);
int save_checkpoint(char *filename, struct Particles *parts, struct LJ *lj, struct Externos *params);
int load_checkpoint(char *filename, struct Particles *parts, struct LJ *lj, struct Externos *params);

#endif
