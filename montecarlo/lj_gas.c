#include "lj_gas.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

// ------------------------------ GENERALES -------------------------------- //

float sub_vector(float* V, int i, float* vec){
  for(int k = 0; k < 3; k++){
    vec[k] = V[3*i+k];
  }
  return 0;
}

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


// -------------------------------- ENERGIA --------------------------------- //


float distancia(float* q1, float* q2, struct LJ *lj, float L){
  double dq2 = 0;
  float dq = 0;
  for(int k = 0; k < 3; k++){
    dq = q1[k] - q2[k];
    if (dq < -0.5*L){
      dq = dq + L;
    }else{
      if (dq > 0.5*L){
        dq = dq - L;
      }
    }
    dq2 = dq2 + dq*dq;
  }
  return dq2;
}

float interaction(float *q1, float *q2, struct LJ *lj, float L){
  float pot = 0;
  float r6 = 0;
  float r2 = distancia(q1, q2, lj, L);
  if (r2 <= lj->rcut2){
    r6 = r2*r2*r2;
    pot = 4*(1.0 - 1.0/r6)/r6 - lj->shift;
  }
  return pot;
}

int energia(struct Particles *parts, struct LJ *lj, float L){
  float kin = 0;
  float pot = 0;
  for(int k = 0; k < 3*parts->n; k++){
    kin = kin + 0.5*parts->p[k]*parts->p[k]/parts->mass;
  }
  for(int i = 1; i < parts->n; i++){
    float qi[3];
    sub_vector(parts->q, i, qi);
    for(int j = 0; j < i; j++){
      float qj[3];
      sub_vector(parts->q, j, qj);
      pot = pot + interaction(qi, qj, lj, L);
    }
  }
  parts->energy = kin + pot;
  parts->kinetic = kin;
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

float delta_energia_pot(struct Particles *parts, struct LJ *lj, float L, float *new_q, float *new_p, int i){
  double delta_pot = 0;
  float qi[3];
  sub_vector(parts->q, i, qi);
  for(int j = 0; j < i; j++){
    float qj[3];
    sub_vector(parts->q, j, qj);
    float pot_ij = interaction(qi, qj, lj, L);
    float new_pot_ij = interaction(new_q, qj, lj, L);
    delta_pot = delta_pot + (new_pot_ij - pot_ij);
  }
  for(int j = i+1; j < parts->n; j++){
    float qj[3];
    sub_vector(parts->q, j, qj);
    float pot_ij = interaction(qi, qj, lj, L);
    float new_pot_ij = interaction(new_q, qj, lj, L);
    delta_pot = delta_pot + (new_pot_ij - pot_ij);
  }
  return delta_pot;
}

// ------------------------------- MUESTREO --------------------------------- //

float set_box(struct Particles *parts, float L){
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

int step(struct Particles *parts, struct LJ *lj, struct Externos *params){
  int res = 1;
  int i = elegir_part(parts->n);
  float new_q[3];
  float new_p[3];
  for(int k = 0; k < 3; k++){
    new_q[k] = parts->q[3*i+k] + (2*uniform()-1)*params->delta_q;
    new_q[k] = new_q[k] + params->L*((new_q[k]<0) - (new_q[k]>params->L));
    new_p[k] = parts->p[3*i+k] + (2*uniform()-1)*params->delta_p;
  }
  float delta_kin = delta_energia_kin(parts, new_p, i);
  float delta_E = delta_kin + delta_energia_pot(parts, lj, params->L, new_q, new_p, i);
  if (delta_E <= 0){
    for(int k = 0; k < 3; k++){
      parts->q[3*i+k] = new_q[k];
      parts->p[3*i+k] = new_p[k];
    }
    parts->energy = parts->energy + delta_E;
    parts->kinetic = parts->kinetic + delta_kin;
  } else {
    float p = exp(-delta_E/params->T);
    if (uniform() < p){
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

int N_steps(struct Particles *parts, struct LJ *lj, struct Externos *params, int Nsteps){
  int aceptados = 0;
  for (int k = 0; k < Nsteps; k++) {
    aceptados = aceptados + step(parts, lj, params);
  }
  return aceptados;
}


int muestrear_impulsos(char *filename, struct Particles *parts, struct LJ *lj, struct Externos *params, int Nsamp, int factor, int factor_term){
  // Termalizacion
  printf("%s\n", filename);
  for(int l = 0; l < factor_term*parts->n; l++){
    step(parts, lj, params);
  }
  // Muestreo
  FILE *f = fopen(filename, "w");
  int aceptados = 0;
  for(int k = 0; k < Nsamp; k++){
    aceptados = aceptados + N_steps(parts, lj, params, factor*parts->n);
    for(int l = 0; l < 3*parts->n; l++){
      fprintf(f, "%f %f\n", parts->q[l], parts->p[l]);
    }
  }
  fclose(f);
  return aceptados;
}

int save_checkpoint(char *filename, struct Particles *parts, struct LJ *lj, struct Externos *params){
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
  fprintf(f, "%f %f\n", lj->rcut2, lj->shift);
  fprintf(f, "%f %f %f %f\n", params->L, params->T, params->delta_q, params->delta_p);
  fclose(f);
  return 0;
}

int load_checkpoint(char *filename, struct Particles *parts, struct LJ *lj, struct Externos *params){
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
  i = fscanf(f, "%f", &lj->rcut2);
  i = fscanf(f, "%f", &lj->shift);
  i = fscanf(f, "%f", &params->L);
  i = fscanf(f, "%f", &params->T);
  i = fscanf(f, "%f", &params->delta_q);
  i = fscanf(f, "%f", &params->delta_p);
  fclose(f);
  return 0;
}



// --------------------------------- MAIN ----------------------------------- //


int main(){

// Particulas
  struct Particles parts;
  int N = 10;
  parts.n = N*N*N;
  parts.mass = 1.043916 * 100; // Masa protón, MeV*(10^-22 s/fm)^2
  parts.q = (float *) malloc(3*parts.n*sizeof(float));
  parts.p = (float *) malloc(3*parts.n*sizeof(float));
  for(int j = 0; j < 3*parts.n; j++){
    parts.q[j] = 0;
    parts.p[j] = 0;
  }

// Potencial
  struct LJ lj;
  lj.rcut2 = 10; // Un ~0.7% del máximo
  float rcut6 = lj.rcut2*lj.rcut2*lj.rcut2;
  lj.shift = 4*(1.0 - 1.0/rcut6)/rcut6;

// Parametros
  struct Externos params;
  params.L = 2*N*6.0; // fm ; mayor a 2*qo*rcut
  params.T = 525; // MeV
  params.delta_q = 6.0/2; // fm
  params.delta_p = 2.067/2; // MeV*10^-22 s/fm

  char filename[255];

  // Distribucion de muchos T; muchas realizaciones
  clock_t start, end;
  double time;
  int factor = 2;
  int factor_term = 500;
  int Nsamp = 20;
  int Nrep = 10;

  float Ts[9] = {3, 2.5, 2, 1.5, 1.0, 0.5, 0.1, 0.05, 0.01};

  for (int j = 0; j < Nrep; j++) {
    srand(j);
    printf("Realizacion %d/%d\n", j+1, Nrep);
    params.T = Ts[0];
    set_box(&parts, params.L);
    set_p(&parts, params.T);
    energia(&parts, &lj, params.L);
    printf("Tramo lineal\n");
    for (int k = 0; k < 9; k++) {
      start = clock();
      params.T = Ts[k];
      sprintf(filename, "LJ_fit/rho0/distribucion_10_rep%d_%f.txt", j+1, params.T);
      int aceptados = muestrear_impulsos(filename, &parts, &lj, &params, Nsamp, factor, factor_term);
      end = clock();
      time = ((double) (end - start)) / CLOCKS_PER_SEC;
      printf("T = %f en %f segundos con %d aceptados (%f)\n", params.T, time, aceptados, ((float) aceptados)/(Nsamp*factor*parts.n));
    }
  }
  free(parts.q);
  free(parts.p);
  return 0;
}
