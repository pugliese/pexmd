#ifndef FG_H
#define FG_H

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

struct Externos {
  float L;
  float T;
  float delta_q;
  float delta_p;
};

float distancia_fases(float* q1, float* q2, float* p1, float* p2, struct Pauli *pauli, float L);
float sub_vector(float* V, int i, float* vec);
float interaction(float *q1, float *q2, float *p1, float *p2, struct Pauli *pauli, float L);
int energia(struct Particles *parts, struct Pauli *pauli, float L);
float delta_energia_kin(struct Particles *parts, float *new_p, int i);
float delta_energia_pot(struct Particles *parts, struct Pauli *pauli, float L, float *new_q, float *new_p, int i);
float uniform();
float boltzmann(float sigma);
int step(struct Particles *parts, struct Pauli *pauli, struct Externos *param);

#endif
