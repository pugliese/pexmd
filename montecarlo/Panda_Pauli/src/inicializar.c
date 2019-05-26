#include "general.h"
#include "celda.h"
#include "interaccion.h"
#include "inicializar.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>


int energia(struct Particles *parts, struct Interaction *pot_tot){
  float kin = 0, pot = 0, rsq, pot_pauli = 0;
  int i, j, mx, my, mz, idx, idxs[3*27];
  int M = parts->M, M3 = M*M*M;
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
  for (int m = 0; m < M3; m++){
    mx = m/(M*M);
    my = (m/M) % M;
    mz = m % M;
    for (int k = 0; k < 14; k++){
      idx = (((mx+idxs[3*k]+M) % M)*M + (my+idxs[3*k+1]+M) % M)*M + (mz+idxs[3*k+2]+M) % M;
      i = parts->primero[m];
      while (i != -1){
        j = parts->primero[idx];
        while (j !=-1 && j != i){
          rsq = distancia(parts->q+3*i, parts->q+3*j, idxs+3*k, parts->l);
          pot = pot + interaction_sin_LUT(parts->type[i], parts->type[j], rsq, parts->p+3*i, parts->p+3*j, pot_tot, &pot_pauli);
          j = parts->siguiente[j];
        }
        i = parts->siguiente[i];
      }
    }
  }
  // Energia cinetica
  for(int i = 0; i < parts->n; i++){
    for(int k = 0; k < 3; k++){
      kin = kin + (parts->p[3*i+k])*(parts->p[3*i+k]);
    }
  }
  parts->kinetic = 0.5*kin/parts->mass;
  parts->energy_panda = pot - pot_pauli;
  parts->energy_pauli = pot_pauli;
  return 0;
}

float set_box(struct Particles *parts, float rcut, float L){
  int n_lado = 1, i = 0, t;
  while (n_lado*n_lado*n_lado < parts->n) n_lado++;
  float dL = L/n_lado;
  for(int x = 0; x < n_lado; x++){
    for(int y = 0; y < n_lado; y++){
      for(int z = 0; z < n_lado; z++){
        parts->q[3*i] = dL*(0.5 + x);
        parts->q[3*i+1] = dL*(0.5 + y);
        parts->q[3*i+2] = dL*(0.5 + z);
        t = (x + y + z)%4;
        parts->type[i] = 2*(t%2) + (t>1);
        i++;
      }
    }
  }
  armar_lista(parts, rcut, L);
  return dL;
}

float set_p(struct Particles *parts, float T){
  float sigma = sqrt(T*parts->mass);
  for(int k = 0; k < 3*parts->n; k++){
    parts->p[k] = boltzmann(sigma);
  }
  float res = 0;
  for(int k = 0; k < 3*parts->n; k++){
    res = res + parts->p[k]*parts->p[k]/parts->mass;
  }
  return res/(3*parts->n);
}

int inicializar(struct Particles *parts, int *comps, int n_types, float mass, float L, float T, struct Interaction *pot_tot){
  parts->n = 0;
  for (int l = 0; l < n_types; l++){
    parts->n += comps[l];
  }
  parts->type = (int *) malloc(parts->n*sizeof(int));
  parts->q = (float *) malloc(3*parts->n*sizeof(float));
  parts->p = (float *) malloc(3*parts->n*sizeof(float));
  parts->mass = mass;
  set_box(parts, pot_tot->rcut, L);
  set_p(parts, T);

  energia(parts, pot_tot);

  return 0;
}

int liberar(struct Particles *parts){
  free(parts->q);
  free(parts->p);
  free(parts->type);

  free(parts->siguiente);
  free(parts->anterior);
  free(parts->primero);
  free(parts->celda);

  return 0;
}
