#ifndef NSM_H
#define NSM_H

#include "math.h"

struct Particles {
  int n;
  float* q;
  float* p;
  float mass;
  float pot_pauli;
  float pot_nuc;
  float pot_coul;
  float kinetic;
};
/*
En este orden son:
N/4 protones spin up
N/4 protones spin down
N/4 neutrones spin up
N/4 neutrones spin down
*/

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

struct Coulomb {
  float lambda;
  float q2;
  float shift;
};

struct Externos {
  float L;
  float T;
  float delta_q;
  float delta_p;
  int ls;
};

float sub_vector(float* V, int i, float* vec);
int elegir_part(int N);
float uniform();
float boltzmann(float sigma);
float distancia_q(float* q1, float* q2);
float distancia_p(float* p1, float* p2);
float interaction_pauli(float r2, float *p1, float *p2, struct Pauli *pauli);
float interaction_nuc(float r, struct Nuclear *nuc);
float interaction_coul(float r, struct Coulomb *coul);
int energia(struct Particles *parts, struct Pauli *pauli, struct Nuclear *nuc, struct Coulomb *coul, float L, int ls);
float delta_energia_kin(struct Particles *parts, float *new_p, int i);
float delta_energia_pot(struct Particles *parts, struct Pauli *pauli, struct Nuclear *nuc, struct Coulomb *coul, float L, int ls, float *new_q, float *new_p, int i, float *deltas);
float set_box(struct Particles *parts, float L);
float set_p(struct Particles *parts, float T);
int step(struct Particles *parts, struct Pauli *pauli, struct Nuclear *nuc, struct Coulomb *coul, struct Externos *params);
int N_steps(struct Particles *parts, struct Pauli *pauli, struct Nuclear *nuc, struct Coulomb *coul, struct Externos *params, int Nsteps);
int muestrear_impulsos(char *filename, struct Particles *parts, struct Pauli *pauli, struct Nuclear *nuc, struct Coulomb *coul, struct Externos *params, int Nsamp, int factor, int factor_term);
int muestrear_energias(char *filename, struct Particles *parts, struct Pauli *pauli, struct Nuclear *nuc, struct Coulomb *coul, struct Externos *params, int Nsamp, int factor, int factor_term);

#endif
