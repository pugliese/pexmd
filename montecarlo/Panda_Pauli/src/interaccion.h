#ifndef I_H
#define I_H

#include "math.h"

struct Interaction{
  struct Pauli* pauli;
  struct Panda_nn* panda_nn;
  struct Panda_np* panda_np;
  float rcut;
};

struct Pauli {
  float qo;
  float po;
  float D;
  float scut2;
  float shift;
};

struct Panda_np {
  float mu_r;
  float mu_a;
  float V_r;
  float V_a;
  float rcut;
  float shift;
};

struct Panda_nn {
  float mu_o;
  float V_o;
  float rcut;
  float shift;
};

float distancia(float *q1, float *q2, int *delta_idx, float l);
float distancia_p(float *p1, float *p2);
float interaction_panda_nn(float r, struct Panda_nn *panda_nn);
float interaction_panda_np(float r, struct Panda_np *panda_np);
float interaction_pauli(float rsq, float *p1, float *p2, struct Pauli *pauli);
float interaction(int t1, int t2, float rsq, float *p1, float *p2, struct Interaction *pot_tot, float *pot_pauli);

#endif
