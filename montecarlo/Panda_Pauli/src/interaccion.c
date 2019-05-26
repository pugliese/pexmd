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
    pot = (panda_np->V_r*exp(-panda_np->mu_r*r) - panda_np->V_a*exp(-panda_np->mu_a*r))/r - panda_np->shift;
  }
  return pot;
}

inline float interaction_pauli(float rsq, float *p1, float *p2, struct Pauli *pauli){
  float pot = 0;
  float s2 = rsq/(pauli->qo*pauli->qo) + distancia_p(p1, p2)/(pauli->po*pauli->po);
  if (s2 <= pauli->scut2){
    pot = pauli->D*exp(-0.5*s2) - pauli->shift;
  }
  return pot;
}

inline float interaction_sin_LUT(int t1, int t2, float rsq, float *p1, float *p2, struct Interaction *pot_tot, float *pot_pauli){
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

// ---------------------------- TABLAS -----------------------------------------

inline float interaction(int t1, int t2, float rsq, float *p1, float *p2, struct Interaction *pot_tot, float *pot_pauli){
  float pot, pot_pauli_aux;
  if (t1/2 == t2/2){
    pot = eval_LUT(rsq, pot_tot->panda_nn->LUT, pot_tot->panda_nn->rcut2, pot_tot->panda_nn->dr2);
    if (t1 == t2){
      float s2 = rsq/(pot_tot->pauli->qo*pot_tot->pauli->qo) + distancia_p(p1, p2)/(pot_tot->pauli->po*pot_tot->pauli->po);
      pot_pauli_aux = eval_LUT(s2, pot_tot->pauli->LUT, pot_tot->pauli->scut2, pot_tot->pauli->ds2);
      *pot_pauli += pot_pauli_aux;
      pot += pot_pauli_aux;
    }
  }else{
    pot = eval_LUT(rsq, pot_tot->panda_np->LUT, pot_tot->panda_np->rcut2, pot_tot->panda_np->dr2);
  }
  return pot;
}

float eval_LUT(float x, float *LUT, float xcut, float dx){
  float res = 0;
  if (x < xcut){
    float a = x/dx;
    int i = (int) a;
    res = (LUT[i] - LUT[i-1])*(a-i) + LUT[i-1];
  }
  return res;
}

float build_LUT_np(struct Panda_np *panda_np, int N){
  float dr2 = panda_np->rcut2/N, r2 = 0, r;
  panda_np->shift = 0;
  panda_np->shift = interaction_panda_np(panda_np->rcut, panda_np);
  panda_np->LUT = (float *) malloc(N*sizeof(float));
  for(int i = 0; i < N; i++){
    r2 = (i+1)*dr2;
    r = sqrt(r2);
    panda_np->LUT[i] = interaction_panda_np(r, panda_np);
  }
  panda_np->dr2 = dr2;
  return dr2;
}

float build_LUT_nn(struct Panda_nn *panda_nn, int N){
  double dr2 = panda_nn->rcut2/N, r2 = 0, r;
  panda_nn->shift = 0;
  panda_nn->shift = interaction_panda_nn(panda_nn->rcut, panda_nn);
  panda_nn->LUT = (float *) malloc(N*sizeof(float));
  for(int i = 0; i < N; i++){
    r2 = (i+1)*dr2;
    r = sqrt(r2);
    panda_nn->LUT[i] = interaction_panda_nn(r, panda_nn);
  }
  panda_nn->dr2 = dr2;
  return dr2;
}

float build_LUT_pauli(struct Pauli *pauli, int N){
  float ds2 = pauli->scut2/N, s2 = 0;
  float p[3] = {0, 0, 0};
  pauli->shift = 0;
  pauli->shift = interaction_pauli(pauli->scut2, p, p, pauli);
  pauli->LUT = (float *) malloc(N*sizeof(float));
  for(int i = 0; i < N; i++){
    s2 = (i+1)*ds2;
    pauli->LUT[i] = interaction_pauli(s2, p, p, pauli);
  }
  pauli->ds2 = ds2;
  return ds2;
}

int build_LUTs(struct Interaction *pot_tot, int Nnp, int Nnn, int Np){
  build_LUT_np(pot_tot->panda_np, Nnp);
  build_LUT_nn(pot_tot->panda_nn, Nnn);
  build_LUT_pauli(pot_tot->pauli, Np);
  return 0;
}

int liberar_LUTs(struct Interaction *pot_tot){
  free(pot_tot->pauli->LUT);
  free(pot_tot->panda_nn->LUT);
  free(pot_tot->panda_np->LUT);

  return 0;
}
