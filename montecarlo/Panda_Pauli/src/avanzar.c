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

float delta_energia_pot(struct Particles *parts, struct Interaction *pot_tot, float *new_q, float *new_p, int i, int *new_ms, int *ms, float* delta_pauli){
  double new_pot = 0, old_pot = 0;
  float new_pot_pauli = 0, old_pot_pauli = 0;
  float rsq;
  int idxs[3];
  int new_idx, idx, j, M = parts->M;
  for(int k = 0; k < 27; k++){
    idxs[0] = k/9 - 1;
    idxs[1] = (k/3)%3 - 1;
    idxs[2] = k%3 - 1;
    idx = (((ms[0] + idxs[0]+M) % M)*M + (ms[1] + idxs[1]+M) % M)*M + (ms[2] + idxs[2]+M) % M;
    j = parts->primero[idx];
    while (j != -1){
      if (j != i){
        rsq = distancia(parts->q+3*i, parts->q+3*j, idxs, parts->l);
        old_pot = old_pot + interaction(parts->type[i], parts->type[j], rsq, parts->p+3*i, parts->p+3*j, pot_tot, &old_pot_pauli);
      }
      j = parts->siguiente[j];
    }
    new_idx = (((new_ms[0]+idxs[0]+M) % M)*M + (new_ms[1]+idxs[1]+M) % M)*M + (new_ms[2]+idxs[2]+M) % M;
    j = parts->primero[new_idx];
    while (j != -1){
      if (j != i){
        rsq = distancia(new_q, parts->q+3*j, idxs, parts->l);
        new_pot = new_pot + interaction(parts->type[i], parts->type[j], rsq, new_p, parts->p+3*j, pot_tot, &new_pot_pauli);
      }
      j = parts->siguiente[j];
    }
  }
  *delta_pauli = new_pot_pauli - old_pot_pauli;
  return new_pot - old_pot;
}


int step(struct Particles *parts, struct Interaction *pot_tot, struct Externos *params){
  int i = elegir_part(parts->n);
  int m = parts->celda[i], M = parts->M;
  int ms[3], new_m, new_ms[3], idxs[3];
  float delta_pauli;
  ms[0] = m/(M*M);
  ms[1] = (m/M) % M;
  ms[2] = m % M;
  float new_q[3], new_p[3];
  for(int k = 0; k < 3; k++){
    new_q[k] = parts->q[3*i+k] + (2*uniform()-1)*params->delta_q;
    idxs[k] = (int) floor(new_q[k]/parts->l);
    new_q[k] = new_q[k] - parts->l*idxs[k];
    new_p[k] = parts->p[3*i+k] + (2*uniform()-1)*params->delta_p;
    new_ms[k] = (ms[k]+idxs[k]+M) % M;
  }
  float delta_kin = delta_energia_kin(parts, new_p, i);
  float delta_E = delta_kin + delta_energia_pot(parts, pot_tot, new_q, new_p, i, new_ms, ms, &delta_pauli);
  int acepto = (delta_E <= 0) || (uniform() < exp(-delta_E/params->T));
  if (acepto){
    new_m = (new_ms[0]*M + new_ms[1])*M + new_ms[2];
    actualizar_lista(parts, i, new_m);
    for(int k = 0; k < 3; k++){
      parts->q[3*i+k] = new_q[k];
      parts->p[3*i+k] = new_p[k];
    }
    parts->kinetic = parts->kinetic + delta_kin;
    parts->energy_pauli = parts->energy_pauli + delta_pauli;
    parts->energy_panda = parts->energy_panda + (delta_E - delta_pauli - delta_kin);
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
