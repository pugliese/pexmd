#include "general.h"
#include "celda.h"
#include "interaccion.h"
#include "avanzar.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

float delta_energia_kin(struct Particles *parts, float *new_p, int i){
  double delta_kin = 0;
  for(int k = 0; k < 3; k++){
    delta_kin = delta_kin + new_p[k]*new_p[k] - parts->p[3*i+k]*parts->p[3*i+k];
  }
  delta_kin = 0.5*delta_kin/parts->mass;
  return delta_kin;
}

float delta_energia_trap(struct Particles *parts, float *new_q, int i, int *new_ms, int *ms, float L, struct Trap *trap){
  if(trap->k != 0){
    float true_new_q[3], true_q[3];
    for(int k = 0; k < 3; k++){
      true_new_q[k] = new_q[k] + new_ms[k]*parts->l - L/2;
      true_q[k] = parts->q[3*i+k] + ms[k]*parts->l - L/2;
    }
    float new_rsq = norma(true_new_q);
    float rsq = norma(true_q);
    return (eval_LUT(new_rsq, trap->LUT, trap->rcut2, trap->dr2) - eval_LUT(rsq, trap->LUT, trap->rcut2, trap->dr2));
  }else{
    return 0;
  }
}

float delta_energia_trap_sin_LUT(struct Particles *parts, float *new_q, int i, int *new_ms, int *ms, float L, struct Trap *trap){
  if(trap->k != 0){
    float true_new_q[3], true_q[3];
    for(int k = 0; k < 3; k++){
      true_new_q[k] = new_q[k] + new_ms[k]*parts->l - L/2;
      true_q[k] = parts->q[3*i+k] + ms[k]*parts->l - L/2;
    }
    float new_rsq = norma(true_new_q);
    float rsq = norma(true_q);
    return (interaction_trap(new_rsq, trap) - interaction_trap(rsq, trap));
  }else{
    return 0;
  }
}

float delta_energia_pot_sin_LUT(struct Particles *parts, struct Interaction *pot_tot, float *new_q, float *new_p, int i, int *new_ms, int *ms){
  float rsq, d = pot_tot->coul->rcut / parts->l, d2 = d*d, dist;
  int idxs[3], temp;
  int new_idx, idx, j, M = parts->M, K_c = (int) ceil(d);
  int neutron_i = parts->type[i]/2, primer_vecino = 1, toca;
  if(neutron_i)  K_c = 1;
  int K = 2*K_c+1, K2 = K*K, K3 = K2*K;
  double delta_Es[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  for(int k = 0; k < K3; k++){
    idxs[0] = k/K2 - K_c;
    idxs[1] = (k/K)%K - K_c;
    idxs[2] = k%K - K_c;
    if(!neutron_i){
      dist = 0;
      for(int l = 0; l < 3; l++){
        temp = max((abs(idxs[l])-1), 0);
        dist += temp*temp;
      }
      if(d2 <= dist) continue;
      primer_vecino = (idxs[0]<=1 && idxs[1]<=1 && idxs[2]<=1 && idxs[0]>=-1 && idxs[1]>=-1 && idxs[2]>=-1);
    }
    reset_energies(pot_tot);
    idx = (((ms[0] + idxs[0] + M) % M)*M + (ms[1] + idxs[1] + M) % M)*M + (ms[2] + idxs[2] + M) % M;
    j = parts->primero[idx];
    while (j != -1){
      toca = primer_vecino || (parts->type[j]/2 == 0);
      if (j != i && toca){
        rsq = distancia(parts->q+3*i, parts->q+3*j, idxs, parts->l);
        interaction_sin_LUT(parts->type[i], parts->type[j], rsq, parts->p+3*i, parts->p+3*j, pot_tot);
      }
      j = parts->siguiente[j];
    }
    delta_Es[0] -= pot_tot->energy_pauli;
    delta_Es[1] -= pot_tot->energy_panda_nn;
    delta_Es[2] -= pot_tot->energy_panda_np;
    delta_Es[3] -= pot_tot->energy_qcnm;
    delta_Es[4] -= pot_tot->energy_coul;
    reset_energies(pot_tot);
    new_idx = (((new_ms[0]+idxs[0]+M) % M)*M + (new_ms[1]+idxs[1]+M) % M)*M + (new_ms[2]+idxs[2]+M) % M;
    j = parts->primero[new_idx];
    while (j != -1){
      toca = primer_vecino || (parts->type[j]/2 == 0);
      if (j != i && toca){
        rsq = distancia(new_q, parts->q+3*j, idxs, parts->l);
        interaction_sin_LUT(parts->type[i], parts->type[j], rsq, new_p, parts->p+3*j, pot_tot);
      }
      j = parts->siguiente[j];
    }
    delta_Es[0] += pot_tot->energy_pauli;
    delta_Es[1] += pot_tot->energy_panda_nn;
    delta_Es[2] += pot_tot->energy_panda_np;
    delta_Es[3] += pot_tot->energy_qcnm;
    delta_Es[4] += pot_tot->energy_coul;
  }
  pot_tot->energy_pauli = delta_Es[0];
  pot_tot->energy_panda_nn = delta_Es[1];
  pot_tot->energy_panda_np = delta_Es[2];
  pot_tot->energy_qcnm = delta_Es[3];
  pot_tot->energy_coul = delta_Es[4];
  return (delta_Es[0] + delta_Es[1] + delta_Es[2] + delta_Es[3] + delta_Es[4]);
}

float delta_energia_pot(struct Particles *parts, struct Interaction *pot_tot, float *new_q, float *new_p, int i, int *new_ms, int *ms){
  float rsq, d = pot_tot->coul->rcut / parts->l, d2 = d*d, dist;
  int idxs[3], temp;
  int new_idx, idx, j, M = parts->M, K_c = (int) ceil(d);
  int neutron_i = parts->type[i]/2, primer_vecino = 1, toca;
  if(neutron_i)  K_c = 1;
  int K = 2*K_c+1, K2 = K*K, K3 = K2*K;
  double delta_Es[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  for(int k = 0; k < K3; k++){
    idxs[0] = k/K2 - K_c;
    idxs[1] = (k/K)%K - K_c;
    idxs[2] = k%K - K_c;
    if(!neutron_i){
      dist = 0;
      for(int l = 0; l < 3; l++){
        temp = max((abs(idxs[l])-1), 0);
        dist += temp*temp;
      }
      if(d2 <= dist) continue;
      primer_vecino = (idxs[0]<=1 && idxs[1]<=1 && idxs[2]<=1 && idxs[0]>=-1 && idxs[1]>=-1 && idxs[2]>=-1);
    }
    reset_energies(pot_tot);
    idx = (((ms[0] + idxs[0] + M) % M)*M + (ms[1] + idxs[1] + M) % M)*M + (ms[2] + idxs[2] + M) % M;
    j = parts->primero[idx];
    while (j != -1){
      toca = primer_vecino || (parts->type[j]/2 == 0);
      if (j != i && toca){
        rsq = distancia(parts->q+3*i, parts->q+3*j, idxs, parts->l);
        interaction(parts->type[i], parts->type[j], rsq, parts->p+3*i, parts->p+3*j, pot_tot);
      }
      j = parts->siguiente[j];
    }
    delta_Es[0] -= pot_tot->energy_pauli;
    delta_Es[1] -= pot_tot->energy_panda_nn;
    delta_Es[2] -= pot_tot->energy_panda_np;
    delta_Es[3] -= pot_tot->energy_qcnm;
    delta_Es[4] -= pot_tot->energy_coul;
    reset_energies(pot_tot);
    new_idx = (((new_ms[0]+idxs[0]+M) % M)*M + (new_ms[1]+idxs[1]+M) % M)*M + (new_ms[2]+idxs[2]+M) % M;
    j = parts->primero[new_idx];
    while (j != -1){
      toca = primer_vecino || (parts->type[j]/2 == 0);
      if (j != i && toca){
        rsq = distancia(new_q, parts->q+3*j, idxs, parts->l);
        interaction(parts->type[i], parts->type[j], rsq, new_p, parts->p+3*j, pot_tot);
      }
      j = parts->siguiente[j];
    }
    delta_Es[0] += pot_tot->energy_pauli;
    delta_Es[1] += pot_tot->energy_panda_nn;
    delta_Es[2] += pot_tot->energy_panda_np;
    delta_Es[3] += pot_tot->energy_qcnm;
    delta_Es[4] += pot_tot->energy_coul;
  }
  pot_tot->energy_pauli = delta_Es[0];
  pot_tot->energy_panda_nn = delta_Es[1];
  pot_tot->energy_panda_np = delta_Es[2];
  pot_tot->energy_qcnm = delta_Es[3];
  pot_tot->energy_coul = delta_Es[4];
  return (delta_Es[0] + delta_Es[1] + delta_Es[2] + delta_Es[3] + delta_Es[4]);
}

int step_sin_LUT(struct Particles *parts, struct Interaction *pot_tot, struct Externos *params){
  int i = elegir_part(parts->n);
  int m = parts->celda[i], M = parts->M;
  int ms[3], new_m, new_ms[3], idx;
  ms[0] = m/(M*M);
  ms[1] = (m/M) % M;
  ms[2] = m % M;
  float new_q[3], new_p[3];
  for(int k = 0; k < 3; k++){
    new_q[k] = parts->q[3*i+k] + (2*uniform()-1)*params->delta_q;
    idx = corregir_celda(new_q+k, parts->l);
    new_p[k] = parts->p[3*i+k] + (2*uniform()-1)*params->delta_p;
    new_ms[k] = (ms[k] + idx + M) % M;
  }
  float delta_kin = delta_energia_kin(parts, new_p, i);
  float delta_trap = delta_energia_trap_sin_LUT(parts, new_q, i, new_ms, ms, params->L, pot_tot->trap);
  float delta_E = delta_kin + delta_trap + delta_energia_pot_sin_LUT(parts, pot_tot, new_q, new_p, i, new_ms, ms);
  int acepto = (delta_E <= 0) || (uniform() < exp(-delta_E/params->T));
  if (acepto){
    new_m = (new_ms[0]*M + new_ms[1])*M + new_ms[2];
    actualizar_lista(parts, i, new_m);
    for(int k = 0; k < 3; k++){
      parts->q[3*i+k] = new_q[k];
      parts->p[3*i+k] = new_p[k];
    }
    parts->kinetic += delta_kin / (parts->mass * parts->n);
    parts->energy_panda_nn += (float) pot_tot->energy_panda_nn / parts->n;
    parts->energy_panda_np += (float) pot_tot->energy_panda_np / parts->n;
    parts->energy_coul += (float) pot_tot->energy_coul / parts->n;
    parts->energy_qcnm += (float) pot_tot->energy_qcnm / parts->n;
    parts->energy_pauli += (float) pot_tot->energy_pauli / parts->n;
  }
  return acepto;
}

int step(struct Particles *parts, struct Interaction *pot_tot, struct Externos *params){
  int i = elegir_part(parts->n);
  int m = parts->celda[i], M = parts->M;
  int ms[3], new_m, new_ms[3], idx;
  ms[0] = m/(M*M);
  ms[1] = (m/M) % M;
  ms[2] = m % M;
  float new_q[3], new_p[3];
  for(int k = 0; k < 3; k++){
    new_q[k] = parts->q[3*i+k] + (2*uniform()-1)*params->delta_q;
    idx = corregir_celda(new_q+k, parts->l);
    new_p[k] = parts->p[3*i+k] + (2*uniform()-1)*params->delta_p;
    new_ms[k] = (ms[k] + idx + M) % M;
  }
  float delta_kin = delta_energia_kin(parts, new_p, i);
  float delta_trap = delta_energia_trap(parts, new_q, i, new_ms, ms, params->L, pot_tot->trap);
  float delta_E = delta_kin + delta_trap + delta_energia_pot(parts, pot_tot, new_q, new_p, i, new_ms, ms);
  int acepto = (delta_E <= 0) || (uniform() < exp(-delta_E/params->T));
  if (acepto){
    new_m = (new_ms[0]*M + new_ms[1])*M + new_ms[2];
    actualizar_lista(parts, i, new_m);
    for(int k = 0; k < 3; k++){
      parts->q[3*i+k] = new_q[k];
      parts->p[3*i+k] = new_p[k];
    }
    parts->kinetic += delta_kin / parts->n;
    parts->energy_panda_nn += (float) (pot_tot->energy_panda_nn / parts->n);
    parts->energy_panda_np += (float) (pot_tot->energy_panda_np / parts->n);
    parts->energy_coul += (float) (pot_tot->energy_coul / parts->n);
    parts->energy_qcnm += (float) (pot_tot->energy_qcnm / parts->n);
    parts->energy_pauli += (float) (pot_tot->energy_pauli / parts->n);
  }
  return acepto;
}

int N_steps(struct Particles *parts, struct Interaction *pot_tot, struct Externos *params, int Nsteps){
  int aceptados = 0;
  for (int k = 0; k < Nsteps; k++) {
    aceptados = aceptados + step(parts, pot_tot, params);
  }
  return aceptados;
}

int N_steps_sin_LUT(struct Particles *parts, struct Interaction *pot_tot, struct Externos *params, int Nsteps){
  int aceptados = 0;
  for (int k = 0; k < Nsteps; k++) {
    aceptados = aceptados + step_sin_LUT(parts, pot_tot, params);
  }
  return aceptados;
}
