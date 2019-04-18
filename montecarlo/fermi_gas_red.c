#include "fermi_gas_red.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

// ------------------------------- GENERALES -------------------------------- //

int elegir_part(int N){
  int i = rand();
  while (i > (RAND_MAX/N)*N){
    i = rand();
  }
  i = i % N;
  return i;
}

float uniform(){
  float p = ((float) rand())/RAND_MAX;
  return p;
}

float boltzmann(float sigma){
  float sum = 0;
  int n_samp = 24;
  for(int k = 0; k < n_samp; k++){
    sum = sum + uniform();
  }
  sum = (sum - 0.5*n_samp)*sqrt(12.0/n_samp)*sigma;
  return sum;
}

// -------------------------------- VECINOS --------------------------------- //

int armar_lista(struct Particles *parts, float rcut, float L){
  int idx, N = parts->n;
  parts->M = (int) floor(L/rcut);
  parts->l = L/parts->M;
  int M3 = parts->M*parts->M*parts->M;
  parts->lista = (int *) malloc((2*N+M3)*sizeof(int));
  parts->celda = (int *) malloc(2*N*sizeof(int));
  for(int i = 0; i < M3; i++){
    parts->lista[2*N+i] = -1;
  }
  for(int i = 0; i < N; i++){
    idx = 0;
    for (int k = 0; k < 3; k++){
      idx = idx*parts->M;
      idx += (int) floor(parts->q[3*i+k]/parts->l);
      parts->q[3*i+k] -= parts->l*floor(parts->q[3*i+k]/parts->l);
    }
    parts->celda[i] = idx;
    parts->lista[i] = parts->lista[2*N + idx];    // Su siguiente es el primero previo (o -1 si no habia)
    parts->lista[N+i] = -1;              // Es el nuevo primero; no tiene anterior
    if (parts->lista[2*N + idx] != -1){
      parts->lista[N + parts->lista[2*N + idx]] = i;  // Si habia primero previo, ahora tiene como anterior a i
    }
    parts->lista[2*N + idx] = i;         // i es el nuevo primero
  }
  return 0;
}

int actualizar_lista(struct Particles *parts, int i, int idx){
  int N = parts->n;
  if (idx!=parts->celda[i]){
    printf("Aca la cagamos?\n");
    // Remuevo i de su celda anterior
    int sig_i = parts->lista[i];
    int ant_i = parts->lista[N+i];
    if (sig_i != -1){
      parts->lista[N+sig_i] = ant_i;
    }
    if (ant_i != -1){
      parts->lista[ant_i] = sig_i;
    }else{
      parts->lista[idx] = sig_i;
    }
    // Lo pongo en su nueva celda
    parts->celda[i] = idx;
    parts->lista[i] = parts->lista[2*N + idx];  // Su siguiente es el primero previo (o -1 si no habia)
    parts->lista[N+i] = -1;            // Es el nuevo primero; no tiene anterior
    if (parts->lista[2*N + idx] != -1){
      parts->lista[parts->lista[2*N + idx]] = i; // Si habia primero previo, ahora tiene como anterior a i
    }
    parts->lista[2*N + idx] = i;       // i es el nuevo primero
  }
  return idx;
}

// -------------------------------- ENERGIA --------------------------------- //


float distancia_fases(float *q1, float *q2, float *p1, float *p2, int *delta_idx, float l, struct Pauli *pauli){
  double dq2 = 0;
  double dp2 = 0;
  for(int k = 0; k < 3; k++){
    dq2 = dq2 + (p1[k] - p2[k] + delta_idx[k]*l)*(p1[k] - p2[k] + delta_idx[k]*l);
    dp2 = dp2 + (p1[k] - p2[k])*(p1[k] - p2[k]);
  }
  dq2 = dq2/(pauli->qo*pauli->qo);
  dp2 = dp2/(pauli->po*pauli->po);
  float s2 = dq2 + dp2;
  return s2;
}

float interaction(float *q1, float *q2, float *p1, float *p2, int *delta_idx, float l, struct Pauli *pauli){
  float pot = 0;
  float s2 = distancia_fases(q1, q2, p1, p2, delta_idx, l, pauli);
  if (s2 < pauli->scut2){
    pot = pauli->D*exp(-0.5*s2) - pauli->shift;
  }
  return pot;
}

int energia(struct Particles *parts, struct Pauli *pauli){
  float kin = 0;
  float pot = 0;
  int i, j, mx, my, mz, idx, idxs[3*14];
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
  idxs[3*13] = 1;
  idxs[3*13+1] = 1;
  idxs[3*13+2] = 1;
  printf("Arme piola los s_k\n");
  for (int m = 0; m < M3; m++){
    mx = m/(M*M);
    my = (m/M)%M;
    mz = m%M;
    for (int k = 0; k < 14; k++){
      idx = (((mx + idxs[3*k])%M)*M + (my + idxs[3*k+1])%M)*M + (mz + idxs[3*k+2])%M;
      printf("El idx anda bien: %d %d\n", m, k);
      i = parts->lista[m];
      while (i != -1){
        j = parts->lista[idx];
        while (j != -1){
          if (j < i){
            pot = pot + interaction(parts->q+3*i, parts->p+3*i, parts->q+3*j, parts->p+3*j, idxs+3*k, parts->l, pauli);
          }
          j = parts->lista[j];
        }
        i = parts->lista[i];
      }
    }
  }
  // Energia cinetica
  for(int i = 0; i < parts->n; i++){
    for(int k = 0; k < 3; k++){
      kin = kin + (parts->p[3*i+k])*(parts->p[3*i+k]);
    }
  }
  parts->energy = kin + pot;
  parts->kinetic = 0.5*kin/parts->mass;
  return 0;
}

float delta_energia_kin(struct Particles *parts, float *new_p, int i){
  double delta_kin = 0;
  for(int k = 0; k < 3; k++){
    delta_kin = delta_kin + new_p[k]*new_p[k] - parts->p[3*i+k]*parts->p[3*i+k];
  }
  delta_kin = 0.5*delta_kin/parts->mass;
  return delta_kin;
}

float delta_energia_pot(struct Particles *parts, struct Pauli *pauli, float *new_q, float *new_p, int i, int new_m){
  double new_pot = 0, old_pot = 0;
  int idxs[3];
  int m = parts->celda[i];
  int new_idx, idx, j, M = parts->M;
  int mx = m/(M*M);
  int my = (m/M)%M;
  int mz = m%M;
  int new_mx = new_m/(M*M);
  int new_my = (new_m/M)%M;
  int new_mz = new_m%M;
  for(int k = 0; k < 27; k++){
    idxs[0] = k/9 - 1;
    idxs[1] = (k/3)%3 - 1;
    idxs[2] = k%3 - 1;
    idx = (((mx + idxs[0])%M)*M + (my + idxs[1])%M)*M + (mz + idxs[2])%M;
    //printf("%d\n", idx);
    j = parts->lista[idx];
    while (j != -1){
      if (j != i){
        //printf("%d\n", j);
        old_pot = old_pot + interaction(parts->q+3*i, parts->p+3*i, parts->q+3*j, parts->p+3*j, idxs, parts->l, pauli);
      }
      j = parts->lista[j];
    }
    new_idx = (((new_mx + idxs[0])%M)*M + (new_my + idxs[1])%M)*M + (new_mz + idxs[2])%M;
    //printf("%d\n", new_idx);
    j = parts->lista[new_idx];
    while (j != -1){
      if (j != i){
        new_pot = new_pot + interaction(new_q, new_p, parts->q+3*j, parts->p+3*j, idxs, parts->l, pauli);
      }
      j = parts->lista[j];
    }
  }
  return new_pot - old_pot;
}

// ------------------------------- MUESTREO --------------------------------- //

float set_box(struct Particles *parts, float rcut, float L){
  int n_lado = pow(parts->n, 1.0/3.0);
  float dL = L/n_lado;
  for(int x = 0; x < n_lado; x++){
    for(int y = 0; y < n_lado; y++){
      for(int z = 0; z < n_lado; z++){
        int i = x*n_lado*n_lado + y*n_lado + z;
        parts->q[3*i] = dL*(0.5 + x);
        parts->q[3*i+1] = dL*(0.5 + y);
        parts->q[3*i+2] = dL*(0.5 + z);
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

int step(struct Particles *parts, struct Pauli *pauli, struct Externos *params){
  int i = elegir_part(parts->n);
  int m = parts->celda[i], res = 1, M = parts->M;
  int mx = m/(M*M);
  int my = (m/M)%M;
  int mz = m%M;
  int new_m, idxs[3];
  float new_q[3];
  float new_p[3];
  for(int k = 0; k < 3; k++){
    new_q[k] = parts->q[3*i+k] + (2*uniform()-1)*params->delta_q;
    idxs[k] = (new_q[k]<0) - (new_q[k]>parts->l);
    new_q[k] = new_q[k] + parts->l*((new_q[k]<0) - (new_q[k]>parts->l));
    new_p[k] = parts->p[3*i+k] + (2*uniform()-1)*params->delta_p;
  }
  new_m = (((mx + idxs[0])%M)*M + (my + idxs[1])%M)*M + (mz + idxs[2])%M;
  float delta_kin = delta_energia_kin(parts, new_p, i);
  float delta_E = delta_kin + delta_energia_pot(parts, pauli, new_q, new_p, i, new_m);
  if (delta_E <= 0){
    actualizar_lista(parts, i, new_m);
    for(int k = 0; k < 3; k++){
      parts->q[3*i+k] = new_q[k];
      parts->p[3*i+k] = new_p[k];
    }
    parts->energy = parts->energy + delta_E;
    parts->kinetic = parts->kinetic + delta_kin;
  } else {
    float p = exp(-delta_E/params->T);
    if (uniform() < p){
      actualizar_lista(parts, i, new_m);
      for(int k = 0; k < 3; k++){
        parts->q[3*i+k] = new_q[k];
        parts->p[3*i+k] = new_p[k];
      }
      parts->energy = parts->energy + delta_E;
      parts->kinetic = parts->kinetic + delta_kin;
    } else {
      res = 0;
    }
  }
  return res;
}

int N_steps(struct Particles *parts, struct Pauli *pauli, struct Externos *params, int Nsteps){
  int aceptados = 0;
  for (int k = 0; k < Nsteps; k++) {
    aceptados = aceptados + step(parts, pauli, params);
  }
  return aceptados;
}
/*

int muestrear_impulsos(char *filename, struct Particles *parts, struct Pauli *pauli, struct Externos *params, int Nsamp, int factor, int factor_term){
  // Termalizacion
  printf("%s\n", filename);
  N_steps(parts, pauli, params, factor_term*parts->n);
  // Muestreo
  FILE *f = fopen(filename, "w");
  int aceptados = 0;
  for(int k = 0; k < Nsamp; k++){
    aceptados = aceptados + N_steps(parts, pauli, params, factor*parts->n);
    for(int l = 0; l < 3*parts->n; l++){
      fprintf(f, "%f %f\n", parts->q[l], parts->p[l]);
    }
  }
  fclose(f);
  return aceptados;
}

int muestrear_energias(char *filename, struct Particles *parts, struct Pauli *pauli, struct Externos *params, int Nsamp, int factor, int factor_term){
  // Termalizacion
  printf("%s\n", filename);
  N_steps(parts, pauli, params, factor_term*parts->n);
  // Muestreo
  FILE *f = fopen(filename, "w");
  int aceptados = 0;
  for(int k = 0; k < Nsamp; k++){
    aceptados = aceptados + N_steps(parts, pauli, params, factor*parts->n+1);
    fprintf(f, "%f %f %d\n", parts->kinetic, parts->energy, aceptados);
  }
  fclose(f);
  return aceptados;
}

int save_checkpoint(char *filename, struct Particles *parts, struct Pauli *pauli, struct Externos *params){
  FILE *f = fopen(filename, "w");
  fprintf(f, "%d %f %f %f\n", parts->n, parts->mass, parts->energy, parts->kinetic);
  for(int k = 0; k < 3*parts->n-1; k++){
    fprintf(f, "%f ", parts->q[k]);
  }
  fprintf(f, "%f\n", parts->q[3*parts->n-1]);
  for(int k = 0; k < 3*parts->n-1; k++){
    fprintf(f, "%f ", parts->p[k]);
  }
  fprintf(f, "%f\n", parts->p[3*parts->n-1]);
  fprintf(f, "%f %f %f %f %f\n", pauli->qo, pauli->po, pauli->D, pauli->scut2, pauli->shift);
  fprintf(f, "%f %f %f %f\n", params->L, params->T, params->delta_q, params->delta_p);
  fclose(f);
  return 0;
}

int load_checkpoint(char *filename, struct Particles *parts, struct Pauli *pauli, struct Externos *params){
  FILE *f = fopen(filename, "r");
  int i = fscanf(f, "%d", &parts->n);
  i = fscanf(f, "%f", &parts->mass);
  i = fscanf(f, "%f", &parts->energy);
  i = fscanf(f, "%f", &parts->kinetic);
  for(int k = 0; k < 3*parts->n; k++){
    i = fscanf(f, "%f", &parts->q[k]);
  }
  for(int k = 0; k < 3*parts->n; k++){
    i = fscanf(f, "%f", &parts->p[k]);
  }
  i = fscanf(f, "%f", &pauli->qo);
  i = fscanf(f, "%f", &pauli->po);
  i = fscanf(f, "%f", &pauli->D);
  i = fscanf(f, "%f", &pauli->scut2);
  i = fscanf(f, "%f", &pauli->shift);
  i = fscanf(f, "%f", &params->L);
  i = fscanf(f, "%f", &params->T);
  i = fscanf(f, "%f", &params->delta_q);
  i = fscanf(f, "%f", &params->delta_p);
  fclose(f);
  return 0;
}

int load_checkpoint_parts(char *filename, struct Particles *parts){
  FILE *f = fopen(filename, "r");
  for(int k = 0; k < 3*parts->n; k++){
    int i = fscanf(f, "%f %f\n", &parts->q[k], &parts->p[k]);
  }
  fclose(f);
  return 0;
}
*/


// --------------------------------- MAIN ----------------------------------- //


int main(int argc, char *argv[]){

  int rep_init = 1;
  if(argc==2){
    int i = sscanf(argv[1], "%d\n", &rep_init);
  }

// Particulas
  struct Particles parts;
  int N = 4;
  parts.n = N*N*N;
  parts.q = (float *) malloc(3*parts.n*sizeof(float));
  parts.p = (float *) malloc(3*parts.n*sizeof(float));
  for(int j = 0; j < 3*parts.n; j++){
    parts.q[j] = 0;
    parts.p[j] = 0;
  }

// Potencial
  struct Pauli pauli;
  // PARAMETROS REDUCIDOS
  pauli.qo = 1;
  pauli.po = 1;
  pauli.D = 0;
  pauli.scut2 = 2.25;
  pauli.shift = pauli.D*exp(-0.5*pauli.scut2);
  parts.mass = 1;


// Parametros
  struct Externos params;
  params.L = N*pauli.qo; // mayor a 2*qo*scut
  params.T = 3; // MeV
  params.delta_q = pauli.qo/2;
  params.delta_p = pauli.po/2;

  set_box(&parts, sqrt(pauli.scut2)*pauli.qo, params.L);
  set_p(&parts, 2);
  energia(&parts, &pauli);
  printf("La caga energia?\n");
  int M3 = parts.M*parts.M*parts.M;
  int A = 2*parts.n + M3;
/*
  int flag = 1, j;
  for(int i = 0; i < M3; i++){
    j = parts.lista[2*parts.n+i];
    printf("Celda %d:\n%d", i, j);
    while(parts.lista[j] != -1){
      printf(" -> %d", parts.lista[j]);
      j = parts.lista[j];
    }
    printf("\n%d", j);
    while(parts.lista[parts.n + j] != -1){
      printf(" <- %d", parts.lista[parts.n + j]);
      j = parts.lista[parts.n+j];
    }
    printf("\n", j);
  }
*/
  int segs = time(NULL);
  int aceptados = N_steps(&parts, &pauli, &params, parts.n/10);
  segs = time(NULL) - segs;
  printf("%f en %d segundos\n", 10*((float) aceptados)/parts.n, segs);

  free(parts.q);
  free(parts.p);
  free(parts.celda);
  free(parts.lista);
  return 0;
}
