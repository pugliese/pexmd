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
  if (s2 < pauli->scut2){
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
      kin = kin + pi[k]*pi[k];
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

int N_steps(struct Particles *parts, struct Pauli *pauli, struct Externos *params, int Nsteps){
  int aceptados = 0;
  for (int k = 0; k < Nsteps; k++) {
    aceptados = aceptados + step(parts, pauli, params);
  }
  return aceptados;
}


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
  struct Pauli pauli;
  float h_barra = 6.582119;
  pauli.qo = 6; // fm
  pauli.po = 2.067; // MeV*10^-22 s/fm
  pauli.D = 34.32*pow(h_barra/(pauli.po*pauli.qo), 3); // MeV
  pauli.scut2 = 10; // Un ~0.7% del máximo
  //pauli.scut2 = 6; // Un ~4.98% del máximo
  pauli.shift = pauli.D*exp(-0.5*pauli.scut2);

// Parametros
  struct Externos params;
  //params.L = 40; // fm ; mayor a 2*qo*scut
  //params.L = 30; // fm ; mayor a 2*qo*scut
  //params.L = 50; // fm ; mayor a 2*qo*scut
  params.L = 2*N*pauli.qo/1; // fm ; mayor a 2*qo*scut
  params.T = 525; // MeV
  params.delta_q = pauli.qo/2; // fm
  params.delta_p = pauli.po/2; // MeV*10^-22 s/fm

  char filename[255];
/*
// Inicializacion
  srand(0);
  set_box(&parts, params.L);
  set_p(&parts, params.T);
  energia(&parts, &pauli, params.L);
  printf("%f\n", parts.kinetic/(1.5*parts.n));
  printf("%f\n", parts.energy - parts.kinetic);
  printf("%f\n", parts.kinetic);
*/
/*
  FILE *f = fopen("prueba_boltzmann.txt", "w");
  for(int i = 0; i < 10000; i++){
    float p = boltzmann(sqrt(params.T*parts.mass));
    fprintf(f, "%f ", p);
  }
  fclose(f);
*/
/*
// Distribucion de muchos T (linscale + logscale)
  clock_t start, end;
  double time;
  int factor = 2;
  int Nsamp = 200;
  int N_temp_lin = 4;
  float lin_Tsup = 7;
  float lin_Tinf = 4;
  float delta_lin = (lin_Tsup - lin_Tinf)/(N_temp_lin - 1);
  params.T = lin_Tsup;
  set_p(&parts, params.T);
  energia(&parts, &pauli, params.L);
  printf("Tramo lineal\n");
  for (int k = 0; k < N_temp_lin; k++) {
    start = clock();
    //params.T = lin_Tsup - delta_lin*(k-0.5);
    //N_steps(&parts, &pauli, &params, 100*parts.n);
    params.T = lin_Tsup - k*delta_lin;
    sprintf(filename, "FD_fit/distribucion_rho0_trial3_%f.txt", params.T);
    int aceptados = muestrear_impulsos(filename, &parts, &pauli, &params, Nsamp, factor);
    end = clock();
    time = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("T = %f en %f segundos con %d aceptados (%f)\n", params.T, time, aceptados, ((float) aceptados)/(Nsamp*factor*parts.n));
  }

  printf("Tramo logaritimico\n");
  int N_temp_log = 5;
  float log10_Tsup = 0.5;
  float log10_Tinf = -1.5;
  float delta_log = (log10_Tsup - log10_Tinf)/(N_temp_log - 1);
  for (int k = 0; k < N_temp_log; k++) {
    start = clock();
    params.T = pow(10, log10_Tsup - k*delta_log);
    sprintf(filename, "FD_fit/distribucion_rho0_trial3_10E%f.txt", log10_Tsup - k*delta_log);
    int aceptados = muestrear_impulsos(filename, &parts, &pauli, &params, Nsamp, factor);
    end = clock();
    time = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("T = %f en %f segundos con %d aceptados (%f)\n", params.T, time, aceptados, ((float) aceptados)/(Nsamp*factor*parts.n));
  }
  */
  /*
  // Distribucion de un solo T; muchas realizaciones
  int factor = 2;
  int Nsamp = 200;
  int Nrea = 6;
  int aceptados;
  clock_t start, end;
  double time, time_global;
  for (int i = 0; i < Nrea; i++) {
    time_global = 0;
    params.T = 3;
    srand(i);
    set_box(&parts, params.L);
    set_p(&parts, params.T);
    energia(&parts, &pauli, params.L);
    for (int k = 0; k < 5; k++) {
      params.T = 3 - 0.5*k;
      start = clock();
      aceptados = N_steps(&parts, &pauli, &params, parts.n*5000);
      end = clock();
      time = ((double) (end - start)) / CLOCKS_PER_SEC;
      printf("T = %f (%d/5) en %f segundos\n", params.T, k+1, time);
      time_global = time_global + time;
    }
    params.T = 0.5;
    start = clock();
    sprintf(filename, "FD_fit/distribucion_rho0_rep%d_%f.txt", i, 0.5);
    aceptados = muestrear_impulsos(filename, &parts, &pauli, &params, Nsamp, factor);
    end = clock();
    time = ((double) (end - start)) / CLOCKS_PER_SEC;
    time_global = time_global + time;
    printf("Iteracion %d/%d para T = %f en %f segundos con %d aceptados (%f)\n", i+1, Nrea, params.T, time_global, aceptados, ((float) aceptados)/(Nsamp*factor*parts.n));
  }
  */
  // Distribucion de muchos T; muchas realizaciones
  clock_t start, end;
  double time;
  int factor = 2;
  int factor_term = 2500;
  int Nsamp = 200;
  int Nrep = 10;

  float Ts[17] = {3, 2.75, 2.5, 2.25, 2, 1.75, 1.5, 1.25, 1.0, 0.75, 0.5, 0.3, 0.1, 0.075, 0.05, 0.03, 0.01};

  for (int j = 0; j < Nrep; j++) {
    srand(j);
    printf("Realizacion %d/%d\n", j+1, Nrep);
    params.T = Ts[0];
    set_box(&parts, params.L);
    set_p(&parts, params.T);
    energia(&parts, &pauli, params.L);
    printf("Tramo lineal\n");
    for (int k = 0; k < 17; k++) {
      start = clock();
      params.T = Ts[k];
      sprintf(filename, "FD_fit/rho2/distribucion_10_rep%d_%f.txt", j+1, params.T);
      int aceptados = muestrear_impulsos(filename, &parts, &pauli, &params, Nsamp, factor, factor_term);
      end = clock();
      time = ((double) (end - start)) / CLOCKS_PER_SEC;
      printf("T = %f en %f segundos con %d aceptados (%f)\n", params.T, time, aceptados, ((float) aceptados)/(Nsamp*factor*parts.n));
    }
  }

  // Distribucion de un solo T
  /*
    for (int k = 0; k < 6; k++) {
      params.T = 3 - 0.5*k;
      N_steps(&parts, &pauli, &params, parts.n*10);
    }
    for (int k = 0; k < 6; k++) {
      params.T = 0.3 - 0.05*k;
      N_steps(&parts, &pauli, &params, parts.n*10);
    }
    for (int k = 0; k < 6; k++) {
      params.T = 0.03 - 0.005*k;
      N_steps(&parts, &pauli, &params, parts.n*10);
    }
    for (int k = 0; k < 5; k++) {
      params.T = 0.003 - 0.0005*k;
      N_steps(&parts, &pauli, &params, parts.n*10);
    }
    params.T = 0.0005;
    sprintf(filename, "FD_fit/FD_fit_n=8_%f.txt", params->T);
    muestrear_impulsos(filename, &parts, &pauli, &params, 1000, 10);
  */

  // Distribucion de muchos T
  /*
    for (int k = 0; k < 4; k++) {
      params.T = 1 - 0.25*k;
      sprintf(filename, "FD_fit/distribucion_rho0_%f.txt", params.T);
      muestrear_impulsos(filename, &parts, &pauli, &params, 1000, 2);
    }
    //save_checkpoint("checkpoint_mx100_Lx1.25_0,5.txt", &parts, &pauli, &params);
    for (int k = 0; k < 6; k++) {
      params.T = 0.1 - 0.025*k;
      sprintf(filename, "FD_fit/distribucion_rho0_%f.txt", params.T);
      muestrear_impulsos(filename, &parts, &pauli, &params, 1000, 2);
    }
    //save_checkpoint("checkpoint_mx100_Lx1.25_0,05.txt", &parts, &pauli, &params);
    for (int k = 0; k < 6; k++) {
      params.T = 0.01 - 0.0025*k;
      sprintf(filename, "FD_fit/distribucion_rho0_%f.txt", params.T);
      muestrear_impulsos(filename, &parts, &pauli, &params, 1000, 2);
    }
    //save_checkpoint("checkpoint_mx100_Lx1.25_0,005.txt", &parts, &pauli, &params);
    for (int k = 0; k < 6; k++) {
      params.T = 0.001 - 0.00025*k;
      sprintf(filename, "FD_fit/distribucion_rho0_%f.txt", params.T);
      muestrear_impulsos(filename, &parts, &pauli, &params, 1000, 2);
    }
    //save_checkpoint("checkpoint_mx100_Lx1.25_0,0005.txt", &parts, &pauli, &params);
    for (int k = 0; k < 6; k++) {
      params.T = 0.0001 - 0.00005*k;
      sprintf(filename, "FD_fit/distribucion_rho0_%f.txt", params.T);
      muestrear_impulsos(filename, &parts, &pauli, &params, 1000, 2);
    }
    //save_checkpoint("checkpoint_mx100_Lx1.25_0,00005.txt", &parts, &pauli, &params);
  */

  free(parts.q);
  free(parts.p);
  return 0;
}
