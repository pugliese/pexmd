#include "interaccion.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

inline float distancia(float *q1, float *q2, int *delta_idx, float l){
  double dq2 = 0;
  for(int k = 0; k < 3; k++){
    dq2 += (q1[k] - q2[k] + delta_idx[k]*l)*(q1[k] - q2[k] + delta_idx[k]*l);
  }
  return dq2;
}

inline float distancia_p(float *p1, float *p2){
  double dp2 = 0;
  for(int k = 0; k < 3; k++){
    dp2 = dp2 + (p1[k] - p2[k])*(p1[k] - p2[k]);
  }
  return dp2;
}


inline float interaction_panda_nn(float r, struct Panda_nn *panda_nn){
  float pot = 0;
  if (r <= panda_nn->rcut){
    pot = panda_nn->V_o*exp(-panda_nn->mu_o*r)/r - panda_nn->shift;
  }
  return pot;
}

inline float interaction_panda_np(float r, struct Panda_np *panda_np){
  float pot = 0;
  if (r <= panda_np->rcut){
    pot = panda_np->V_r*exp(-panda_np->mu_r*r)/r - panda_np->V_a*exp(-panda_np->mu_a*r)/r - panda_np->shift;
  }
  return pot;
}

inline float interaction_pauli(float rsq, float *p1, float *p2, struct Pauli *pauli){
  float pot = 0;
  float s2 = rsq/(pauli->qo*pauli->qo) + distancia_p(p1, p2)/(pauli->po*pauli->po);
  if (s2 < pauli->scut2){
    pot = pauli->D*exp(-0.5*s2) - pauli->shift;
  }
  return pot;
}

inline float interaction(int t1, int t2, float rsq, float *p1, float *p2, struct Interaction *pot_tot, float *pot_pauli){
  float pot, pot_pauli_aux;
  float r = sqrt(rsq);
  if (t1/2 == t2/2){
    pot = interaction_panda_nn(r, pot_tot->panda_nn);
    if (t1 == t2){
      pot_pauli_aux = interaction_pauli(rsq, p1, p2, pot_tot->pauli);
      *pot_pauli += pot_pauli_aux;
      pot += pot_pauli_aux;
    }
  }else{
    pot = interaction_panda_np(r, pot_tot->panda_np);
  }
  return pot;
}
