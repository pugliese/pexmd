#ifndef NM_H
#define NM_H

#include "math.h"

struct Particles {
  int n;
  float* q;
  float* p;
  float mass;
  float energy;
  float kinetic;
};

struct Pauli {
  float qo;
  float po;
  float D;
  float scut2;
  float shift;
};

struct Nuclear {
  float Vo;
  float r1;
  float r2;
  float p1;
  float p2;
  float d;
  float a;
  float rcut;
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
float distancia_q(float* q1, float* q2, float L);
float distancia_p(float* p1, float* p2);
float interaction_pauli(float r2, float *p1, float *p2, struct Pauli *pauli);
float interaction_nuc(float r, struct Nuclear *nuc);
int energia(struct Particles *parts, struct Pauli *pauli, struct Nuclear *nuc, float L);
float delta_energia_kin(struct Particles *parts, float *new_p, int i);
float delta_energia_pot(struct Particles *parts, struct Pauli *pauli, struct Nuclear *nuc, float L, float *new_q, float *new_p, int i);
float set_box(struct Particles *parts, float L);
float set_p(struct Particles *parts, float T);
int step(struct Particles *parts, struct Pauli *pauli, struct Nuclear *nuc, struct Externos *params);
int N_steps(struct Particles *parts, struct Pauli *pauli, struct Nuclear *nuc, struct Externos *params, int Nsteps);
int muestrear_impulsos(char *filename, struct Particles *parts, struct Pauli *pauli, struct Nuclear *nuc, struct Externos *params, int Nsamp, int factor, int factor_term);
int muestrear_energias(char *filename, struct Particles *parts, struct Pauli *pauli, struct Nuclear *nuc, struct Externos *params, int Nsamp, int factor, int factor_term);

#endif
