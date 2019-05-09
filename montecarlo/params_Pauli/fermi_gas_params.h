#ifndef FGP_H
#define FGP_H

#include "math.h"

struct Particles {
  int n;
  float* q;
  float* p;
  float mass;
  float energy;
  float kinetic;
  int* siguiente;
  int* anterior;
  int* primero;
  int* celda;
	int M;
  float l;
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

int armar_lista(struct Particles *parts, float rcut, float L);
int actualizar_lista(struct Particles *parts, int i, int prev_idx);
float distancia_fases(float *q1, float *q2, float *p1, float *p2, int *delta_idx, float l, struct Pauli *pauli);
float interaction(float *q1, float *q2, float *p1, float *p2, int *delta_idx, float l, struct Pauli *pauli);
int energia(struct Particles *parts, struct Pauli *pauli);
float delta_energia_kin(struct Particles *parts, float *new_p, int i);
float delta_energia_pot(struct Particles *parts, struct Pauli *pauli, float *new_q, float *new_p, int i, int *new_ms, int *ms);
float uniform();
float boltzmann(float sigma);
int step(struct Particles *parts, struct Pauli *pauli, struct Externos *param);

#endif
