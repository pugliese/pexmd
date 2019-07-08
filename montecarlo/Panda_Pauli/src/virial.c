#include "general.h"
#include "celda.h"
#include "interaccion.h"
#include "avanzar.h"
#include "virial.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

inline float force_mod_panda_nn(float r, struct Panda_nn *panda_nn){
  float f_mod = 0;
  if (r <= panda_nn->rcut){
    f_mod = panda_nn->V_o*exp(-panda_nn->mu_o*r)*(1+panda_nn->mu_o*r)/(r*r*r);
  }
  return f_mod;
}

inline float force_mod_panda_np(float r, struct Panda_np *panda_np){
  float f_mod = 0;
  if (r <= panda_np->rcut){
    f_mod = (panda_np->V_r*exp(-panda_np->mu_r*r)*(1+panda_np->mu_r*r) - panda_np->V_a*exp(-panda_np->mu_a*r)*(1+panda_np->mu_a*r))/(r*r*r);
  }
  return f_mod;
}

inline float fgorce_mod_pauli(float rsq, float *p1, float *p2, struct Pauli *pauli, float *gorce_mod){
  float f_mod;
  float r2_red = rsq/(pauli->qo*pauli->qo);
  if (r2_red < pauli->scut2){
    float p2_red = distancia_p(p1, p2)/(pauli->po*pauli->po);
    float exp_q = exp(-0.5*r2_red);
    float exp_p = exp(-0.5*p2_red);
    f_mod = pauli->D*exp_p*exp_q/(pauli->qo*pauli->qo);
    *gorce_mod = (pauli->D*exp_q - pauli->shift)*exp_p/(pauli->po*pauli->po);
  }else{
    f_mod = 0;
    *gorce_mod = 0;
  }
  return f_mod;
}

inline float fgorce_mod(int t1, int t2, float rsq, float *p1, float *p2, struct Interaction *pot_tot, float *gorce_mod){
  float force_mod;
  *gorce_mod = 0;
  float r = sqrt(rsq);
  if (t1/2 == t2/2){
    force_mod = force_mod_panda_nn(r, pot_tot->panda_nn);
    if (t1 == t2){
      force_mod += fgorce_mod_pauli(rsq, p1, p2, pot_tot->pauli, gorce_mod);
    }
  }else{
    force_mod = force_mod_panda_np(r, pot_tot->panda_np);
  }
  return force_mod;
}

float presion_temperatura(struct Particles *parts, struct Interaction *pot_tot, struct Externos *params, float *PyT){
  double qF = 0, qp = 0;
  int M = parts->M, M3 = M*M*M;
  int mx, my, mz, idxs[3*14], i, j, idx;
  for (int k = 0; k < 9; k++){
    idxs[3*k] = -1;
    idxs[3*k+1] = k/3 - 1;
    idxs[3*k+2] = k%3 - 1;
  }
  for (int k = 9; k < 13; k++){
    idxs[3*k] = 0;
    idxs[3*k+1] = (k-9)/2;
    idxs[3*k+2] = (k-9)%2;
  }
  idxs[3*13] = 0;
  idxs[3*13+1] = 1;
  idxs[3*13+2] = -1;
  float delta_q[3], delta_p[3], force_mod, gorce_mod = 0, rsq;
  for (int m = 0; m < M3; m++){
    mx = m/(M*M);
    my = (m/M) % M;
    mz = m % M;
    i = parts->primero[m];
    for (int k = 0; k < 14; k++){
      idx = (((mx + idxs[3*k] + M)%M)*M + (my + idxs[3*k+1] + M)%M)*M + (mz + idxs[3*k+2] + M)%M;
      i = parts->primero[m];
      while (i != -1){
        j = parts->primero[idx];
        while (j !=-1 && j != i){
          rsq = 0;
          for(int l = 0; l < 3; l++){
            delta_q[l] = parts->q[3*i+l] - parts->q[3*j+l] - idxs[3*k+l]*parts->l;
            delta_p[l] = parts->p[3*i+l] - parts->p[3*j+l];
            rsq += delta_q[l]*delta_q[l];
          }
          force_mod = fgorce_mod(parts->type[i], parts->type[j], rsq, parts->p+3*i, parts->p+3*j, pot_tot, &gorce_mod);
          for (int l = 0; l < 3; l++){
            qF += force_mod*delta_q[l]*delta_q[l];
            qp -= gorce_mod*delta_p[l]*delta_p[l];
          }
          j = parts->siguiente[j];
        }
        i = parts->siguiente[i];
      }
    }
  }
  for (int l = 0; l < 3*parts->n; l++){
    qp += parts->p[l]*parts->p[l]/parts->mass;
  }
  PyT[0] = (qF+qp)/(3*params->L*params->L*params->L);
  PyT[1] = qp/(3*parts->n);
  return 0;
}


float fgorces(struct Particles *parts, struct Interaction *pot_tot, struct Externos *params, float *forces, float *gorces){
  int M = parts->M, M3 = M*M*M;
  int mx, my, mz, idxs[3*14], i, j, idx;
  for (int k = 0; k < 9; k++){
    idxs[3*k] = -1;
    idxs[3*k+1] = k/3 - 1;
    idxs[3*k+2] = k%3 - 1;
  }
  for (int k = 9; k < 13; k++){
    idxs[3*k] = 0;
    idxs[3*k+1] = (k-9)/2;
    idxs[3*k+2] = (k-9)%2;
  }
  idxs[3*13] = 0;
  idxs[3*13+1] = 1;
  idxs[3*13+2] = -1;
  float delta_q[3], delta_p[3], force_mod, gorce_mod = 0, rsq;
  for (int l = 0; l < 3*parts->n; l++){
    forces[l] = 0;
    gorces[l] = 0;
  }
  for (int m = 0; m < M3; m++){
    mx = m/(M*M);
    my = (m/M) % M;
    mz = m % M;
    i = parts->primero[m];
    for (int k = 0; k < 14; k++){
      idx = (((mx + idxs[3*k] + M)%M)*M + (my + idxs[3*k+1] + M)%M)*M + (mz + idxs[3*k+2] + M)%M;
      i = parts->primero[m];
      while (i != -1){
        j = parts->primero[idx];
        while (j !=-1 && j != i){
          rsq = 0;
          for(int l = 0; l < 3; l++){
            delta_q[l] = parts->q[3*i+l] - parts->q[3*j+l] - idxs[3*k+l]*parts->l;
            delta_p[l] = parts->p[3*i+l] - parts->p[3*j+l];
            rsq += delta_q[l]*delta_q[l];
          }
          force_mod = fgorce_mod(parts->type[i], parts->type[j], rsq, parts->p+3*i, parts->p+3*j, pot_tot, &gorce_mod);
          for (int l = 0; l < 3; l++){
            forces[3*i+l] += force_mod*delta_q[l];
            forces[3*j+l] -= force_mod*delta_q[l];
            gorces[3*i+l] += gorce_mod*delta_p[l];
            gorces[3*j+l] -= gorce_mod*delta_p[l];
          }
          j = parts->siguiente[j];
        }
        i = parts->siguiente[i];
      }
    }
  }
  return 0;
}


float punto_fijo(struct Particles *parts, struct Interaction *pot_tot, struct Externos *params, float *new_p){
  int M = parts->M, M3 = M*M*M;
  int mx, my, mz, idxs[3*27], i, j, idx;
  for (int k = 0; k < 27; k++){
    idxs[3*k] = k/9 - 1;
    idxs[3*k+1] = (k/3)%3 - 1;
    idxs[3*k+2] = k%3 - 1;
  }
  float delta_q[3], rsq, s2_ij, sum_i = 0, D = parts->mass*pot_tot->pauli->D/(pot_tot->pauli->po*pot_tot->pauli->po);
  float *exp_ij = (float *) malloc(parts->n*sizeof(float));
  for (int l = 0; l < 3*parts->n; l++)  new_p[l] = 0;
  for (int m = 0; m < M3; m++){
    mx = m/(M*M);
    my = (m/M) % M;
    mz = m % M;
    i = parts->primero[m];
    for (int k = 0; k < 14; k++){
      idx = (((mx + idxs[3*k] + M)%M)*M + (my + idxs[3*k+1] + M)%M)*M + (mz + idxs[3*k+2] + M)%M;
      i = parts->primero[m];
      while (i != -1){
        for (int l = 0; l < parts->n; l++)  exp_ij[l] = 0;
        j = parts->primero[idx];
        while (j !=-1){
          if (i==j){
            j = parts->siguiente[j];
            continue;
          }
          for(int l = 0; l < 3; l++){
            delta_q[l] = parts->q[3*i+l] - parts->q[3*j+l] + idxs[3*k+l]*parts->l;
          }
          rsq = norma(delta_q);
          s2_ij = rsq/(pot_tot->pauli->qo*pot_tot->pauli->qo) + distancia_p(parts->p+3*i, parts->p+3*j)/(pot_tot->pauli->po*pot_tot->pauli->po);
          exp_ij[j] = D*exp(-0.5*s2_ij);
          j = parts->siguiente[j];
        }
        sum_i = 0;
        for(int j = 0; j < parts->n; j++){
          sum_i += exp_ij[j];
        }
        for(int j = 0; j < parts->n; j++){
          for(int k = 0; k < 3; k++){
            new_p[3*i+k] += parts->p[3*j+k]*exp_ij[j]/(1+sum_i);
          }
        }
        i = parts->siguiente[i];
      }
    }
  }
  free(exp_ij);
  return 0;
}
