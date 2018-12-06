#include "fermi_gas.h"
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


float distancia_fases(float* q1, float* q2, float* p1, float* p2, struct Pauli *pauli, float L){
  double dq2 = 0;
  double dp2 = 0;
  float dq = 0;
  for(int k = 0; k < 3; k++){
    dp2 = dp2 + (p1[k] - p2[k])*(p1[k] - p2[k]);
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
  dq2 = dq2/(pauli->qo*pauli->qo);
  dp2 = dp2/(pauli->po*pauli->po);
  float s2 = dq2 + dp2;
  return s2;
}

float interaction(float *q1, float *q2, float *p1, float *p2, struct Pauli *pauli, float L){
  float pot = 0;
  float s2 = distancia_fases(q1, q2, p1, p2, pauli, L);
  if (s2 <= pauli->scut2){
    pot = pauli->D*exp(-0.5*s2) - pauli->shift;
  }
  return pot;
}

int energia(struct Particles *parts, struct Pauli *pauli, float L){
  float kin = 0;
  float pot = 0;
  for(int i = 0; i < parts->n; i++){
    float qi[3];
    float pi[3];
    sub_vector(parts->q, i, qi);
    sub_vector(parts->p, i, pi);
    for(int k = 0; k < 3; k++){
      kin = kin + 0.5*pi[k]*pi[k]/parts->mass;
    }
    for(int j = 0; j < i; j++){
      float qj[3];
      float pj[3];
      sub_vector(parts->q, j, qj);
      sub_vector(parts->p, j, pj);
      pot = pot + interaction(qi, qj, pi, pj, pauli, L);
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

float delta_energia_pot(struct Particles *parts, struct Pauli *pauli, float L, float *new_q, float *new_p, int i){
  double delta_pot = 0;
  float qi[3];
  float pi[3];
  sub_vector(parts->q, i, qi);
  sub_vector(parts->p, i, pi);
  for(int j = 0; j < i; j++){
    float qj[3];
    float pj[3];
    sub_vector(parts->q, j, qj);
    sub_vector(parts->p, j, pj);
    float pot_ij = interaction(qi, qj, pi, pj, pauli, L);
    float new_pot_ij = interaction(new_q, qj, new_p, pj, pauli, L);
    delta_pot = delta_pot + (new_pot_ij - pot_ij);
  }
  for(int j = i+1; j < parts->n; j++){
    float qj[3];
    float pj[3];
    sub_vector(parts->q, j, qj);
    sub_vector(parts->p, j, pj);
    float pot_ij = interaction(qi, qj, pi, pj, pauli, L);
    float new_pot_ij = interaction(new_q, qj, new_p, pj, pauli, L);
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

int step(struct Particles *parts, struct Pauli *pauli, struct Externos *params){
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
  float delta_E = delta_kin + delta_energia_pot(parts, pauli, params->L, new_q, new_p, i);
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


int muestrear_impulsos(struct Particles *parts, struct Pauli *pauli, struct Externos *params, int Nsamp, int factor){
  // Termalizacion
  for(int l = 0; l < 10*factor*parts->n; l++){
    step(parts, pauli, params);
  }
  // Muestreo
  char filename[255];
  sprintf(filename, "impulsos/impulsos_mx100_Lx1.25_%f.txt", params->T);
  printf("%s\n", filename);
  FILE *f = fopen(filename, "w");
  for(int k = 0; k < Nsamp; k++){
    for(int l = 0; l < factor*parts->n; l++){
      step(parts, pauli, params);
    }
    for(int l = 0; l < 3*parts->n; l++){
      fprintf(f, "%f ", parts->p[l]);
    }
  }
  fclose(f);
  return 0;
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
  return 0;
}



// --------------------------------- MAIN ----------------------------------- //


int main(){

// Particulas
  struct Particles parts;
  parts.n = 1000;
  parts.mass = 1.043916 * 100; // Masa protón, MeV*(10^-22 s/fm)^2
  parts.q = (float *) malloc(3*parts.n*sizeof(float));
  parts.p = (float *) malloc(3*parts.n*sizeof(float));
  for(int j = 0; j < 3*parts.n; j++){
    parts.q[j] = 0;
    parts.p[j] = 0;
  }

// Potencial
  struct Pauli pauli;
  float h_barra = 6.582119;
  pauli.qo = 2.067; // fm
  pauli.po = 6; // MeV*10^-22 s/fm
  pauli.D = 34.32*pow(h_barra/(pauli.po*pauli.qo), 3); // MeV
  pauli.scut2 = 10; // Un ~0.7% del máximo
  //pauli.scut2 = 6; // Un ~4.98% del máximo
  pauli.shift = pauli.D*exp(-0.5*pauli.scut2);

// Parametros
  struct Externos params;
  //params.L = 40; // fm ; mayor a 2*qo*scut
  //params.L = 30; // fm ; mayor a 2*qo*scut
  params.L = 50; // fm ; mayor a 2*qo*scut
  params.T = 30; // MeV
  params.delta_q = 3; // fm
  params.delta_p = 1; // MeV*10^-22 s/fm

// Inicializacion
  srand(0);
  set_box(&parts, params.L);
  set_p(&parts, params.T);
  energia(&parts, &pauli, params.L);


  for (int k = 0; k < 6; k++) {
    params.T = 30 - 5*k;
    muestrear_impulsos(&parts, &pauli, &params, 100, 1);
  }
  save_checkpoint("checkpoint_mx100_Lx1.25_5.txt", &parts, &pauli, &params);
  for (int k = 0; k < 6; k++) {
    params.T = 3 - 0.5*k;
    muestrear_impulsos(&parts, &pauli, &params, 100, 1);
  }
  save_checkpoint("checkpoint_mx100_Lx1.25_0,5.txt", &parts, &pauli, &params);
  for (int k = 0; k < 6; k++) {
    params.T = 0.3 - 0.05*k;
    muestrear_impulsos(&parts, &pauli, &params, 200, 3);
  }
  save_checkpoint("checkpoint_mx100_Lx1.25_0,05.txt", &parts, &pauli, &params);
  for (int k = 0; k < 6; k++) {
    params.T = 0.03 - 0.005*k;
    muestrear_impulsos(&parts, &pauli, &params, 200, 3);
  }
  save_checkpoint("checkpoint_mx100_Lx1.25_0,005.txt", &parts, &pauli, &params);
  for (int k = 0; k < 6; k++) {
    params.T = 0.003 - 0.0005*k;
    muestrear_impulsos(&parts, &pauli, &params, 400, 10);
  }
  save_checkpoint("checkpoint_mx100_Lx1.25_0,0005.txt", &parts, &pauli, &params);
  for (int k = 0; k < 6; k++) {
    params.T = 0.0003 - 0.00005*k;
    muestrear_impulsos(&parts, &pauli, &params, 400, 10);
  }
  save_checkpoint("checkpoint_mx100_Lx1.25_0,00005.txt", &parts, &pauli, &params);

/*
  load_checkpoint("checkpoint_mx1000_0,00005.txt", &parts, &pauli, &params);
  for (int k = 0; k < 6; k++) {
    params.T = 0.00003 - 0.000005*k;
    muestrear_impulsos(&parts, &pauli, &params, 1000, 10);
  }
  save_checkpoint("checkpoint_mx1000_0,000005.txt", &parts, &pauli, &params);
*/
  free(parts.q);
  free(parts.p);
  return 0;
}
