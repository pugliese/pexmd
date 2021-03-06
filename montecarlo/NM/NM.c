#include "NM.h"
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


float distancia_q(float* q1, float* q2, float L){
  float dq2 = 0;
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
    dq2 += dq*dq;
  }
  return dq2;
}

float distancia_p(float* p1, float* p2){
  float dp2 = 0;
  for(int k = 0; k < 3; k++){
    dp2 = dp2 + (p1[k] - p2[k])*(p1[k] - p2[k]);
  }
  return dp2;
}

float interaction_pauli(float r2, float *p1, float *p2, struct Pauli *pauli){
  float pot = 0;
  float s2 = r2/(pauli->qo*pauli->qo) + distancia_p(p1, p2)/(pauli->po*pauli->po);
  if (s2 < pauli->scut2){
    pot = (pauli->D)*exp(-0.5*s2) - pauli->shift;
    //printf("%f %f: %f %f %f\n", s2, pauli->scut2, pot, pauli->D*exp(-0.5*s2), pauli->shift);
  }
  return pot;
}

float interaction_nuc(float r, struct Nuclear *nuc){
  float pot = 0;
  if (r < nuc->rcut){
    float r1_pow = pow(nuc->r1/r, nuc->p1);
    float r2_pow = pow(nuc->r2/r, nuc->p2);
    pot = nuc->Vo*(r1_pow - r2_pow)/(1 + exp((r - nuc->d)/nuc->a)) - nuc->shift;
  }
  return pot;
}

int energia(struct Particles *parts, struct Pauli *pauli, struct Nuclear *nuc, float L){
  double kin = 0;
  double pot_n = 0;
  double pot_p = 0;
  float r, rsq;
  int m = parts->n/4;
  int inf, sup, g;
  for (int i = 1; i < parts->n; i++) {
    // Energia cinetica
    float qi[3];
    float pi[3];
    sub_vector(parts->q, i, qi);
    sub_vector(parts->p, i, pi);
    for (int k = 0; k < 3; k++){
      kin += pi[k]*pi[k];
    }
    // Interaccion
    g = i/m;    // Grupo de la particula i
    inf = g*m;      // Limites
    sup = (g+1)*m;  // del grupo
    for (int j = 0; j < i; j++){
      float qj[3];
      sub_vector(parts->q, j, qj);
      rsq = distancia_q(qi, qj, L);
      r = sqrt(rsq);
      pot_n += interaction_nuc(r, nuc);
      if ((inf <= j) && (j < sup)){
        float pj[3];
        sub_vector(parts->p, j, pj);
        pot_p += ((double) interaction_pauli(rsq, pi, pj, pauli));
        //printf("%f\n", pot_p);
      }
    }
  }
  parts->kinetic = ((float) 0.5*kin/parts->mass);
  parts->pot_nuc = ((float) pot_n);
  parts->pot_pauli = ((float) pot_p);
  return parts->kinetic + parts->pot_nuc + parts->pot_pauli;
}

float delta_energia_kin(struct Particles *parts, float *new_p, int i){
  double delta_kin = 0;
  for(int k = 0; k < 3; k++){
    delta_kin += new_p[k]*new_p[k] - parts->p[3*i+k]*parts->p[3*i+k];
  }
  delta_kin = 0.5*delta_kin/parts->mass;
  return ((float) delta_kin);
}

float delta_energia_pot(struct Particles *parts, struct Pauli *pauli, struct Nuclear *nuc, float L, float *new_q, float *new_p, int i, float *deltas){
  double delta_nuc = 0;
  double delta_pauli = 0;
  float qi[3];
  float pi[3];
  float r, rsq, new_r, new_rsq, pot_ij, new_pot_ij;
  sub_vector(parts->q, i, qi);
  sub_vector(parts->p, i, pi);
  int m = parts->n/4;
  int g = i/m;
  int inf = g*m;
  int sup = (g+1)*m;
  for (int j = 0; j < parts->n; j++) {
    if (j == i) continue;

    float qj[3];
    sub_vector(parts->q, j, qj);

    rsq = distancia_q(qi, qj, L);
    r = sqrt(rsq);
    pot_ij = interaction_nuc(r, nuc);

    new_rsq = distancia_q(new_q, qj, L);
    new_r = sqrt(new_rsq);
    new_pot_ij = interaction_nuc(new_r, nuc);
    delta_nuc += (new_pot_ij - pot_ij);

    if ((inf <= j) && (j < sup)){
      float pj[3];
      sub_vector(parts->p, j, pj);

      pot_ij = interaction_pauli(rsq, pi, pj, pauli);
      new_pot_ij = interaction_pauli(new_rsq, new_p, pj, pauli);
      delta_pauli += (new_pot_ij - pot_ij);
    }
  }
  deltas[0] = (float) delta_nuc;
  deltas[1] = (float) delta_pauli;
  return ((float) (delta_nuc + delta_pauli));
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

int step(struct Particles *parts, struct Pauli *pauli, struct Nuclear *nuc, struct Externos *params){
  int res = 1;
  int i = elegir_part(parts->n);
  float deltas[2];
  float new_q[3];
  float new_p[3];
  for(int k = 0; k < 3; k++){
    new_q[k] = parts->q[3*i+k] + (2*uniform()-1)*params->delta_q;
    new_q[k] = new_q[k] + params->L*((new_q[k]<0) - (new_q[k]>params->L));
    new_p[k] = parts->p[3*i+k] + (2*uniform()-1)*params->delta_p;
  }

  float delta_kin = delta_energia_kin(parts, new_p, i);
  float delta_E = delta_kin + delta_energia_pot(parts, pauli, nuc, params->L, new_q, new_p, i, deltas);

  if (delta_E <= 0){
    for(int k = 0; k < 3; k++){
      parts->q[3*i+k] = new_q[k];
      parts->p[3*i+k] = new_p[k];
    }
    parts->kinetic += delta_kin;
    parts->pot_nuc += deltas[0];
    parts->pot_pauli += deltas[1];
  } else {
    float p = exp(-delta_E/params->T);
    if (uniform() < p){
      for(int k = 0; k < 3; k++){
        parts->q[3*i+k] = new_q[k];
        parts->p[3*i+k] = new_p[k];
      }
      parts->kinetic += delta_kin;
      parts->pot_nuc += deltas[0];
      parts->pot_pauli += deltas[1];
    } else {
      res = 0;
    }
  }
  return res;
}

int N_steps(struct Particles *parts, struct Pauli *pauli, struct Nuclear *nuc, struct Externos *params, int Nsteps){
  int aceptados = 0;
  float E_vieja, kin_vieja;
  for (int k = 0; k < Nsteps; k++) {
    aceptados = aceptados + step(parts, pauli, nuc, params);
  }
  return aceptados;
}


int muestrear_impulsos(char *filename, struct Particles *parts, struct Pauli *pauli, struct Nuclear *nuc, struct Externos *params, int Nsamp, int factor, int factor_term){
  // Termalizacion
  printf("%s\n", filename);
  N_steps(parts, pauli, nuc, params, factor_term);
  // Muestreo
  FILE *f = fopen(filename, "w");
  int aceptados = 0;
  for(int k = 0; k < Nsamp; k++){
    aceptados = aceptados + N_steps(parts, pauli, nuc, params, factor*parts->n+1);
    for(int l = 0; l < 3*parts->n; l++){
      fprintf(f, "%f %f\n", parts->q[l], parts->p[l]);
    }
  }
  fclose(f);
  return aceptados;
}

int muestrear_energias(char *filename, struct Particles *parts, struct Pauli *pauli, struct Nuclear *nuc, struct Externos *params, int Nsamp, int factor, int factor_term){
  // Termalizacion
  printf("%s\n", filename);
  N_steps(parts, pauli, nuc, params, factor_term*parts->n);
  // Muestreo
  FILE *f = fopen(filename, "w");
  int aceptados = 0;
  for(int k = 0; k < Nsamp; k++){
    aceptados = aceptados + N_steps(parts, pauli, nuc, params, factor*parts->n+1);
    fprintf(f, "%f %f %f %d\n", parts->kinetic, parts->pot_nuc, parts->pot_pauli, aceptados);
  }
  fclose(f);
  return aceptados;
}


int save_checkpoint(char *filename, struct Particles *parts, struct Pauli *pauli, struct Nuclear *nuc, struct Externos *params){
  FILE *f = fopen(filename, "w");
  printf("%s\n", filename);
  fprintf(f, "%d %f %f %f %f\n", parts->n, parts->mass, parts->kinetic, parts->pot_nuc, parts->pot_pauli);
  for(int k = 0; k < 3*parts->n-1; k++){
    fprintf(f, "%f ", parts->q[k]);
  }
  fprintf(f, "%f\n", parts->q[3*parts->n-1]);
  for(int k = 0; k < 3*parts->n-1; k++){
    fprintf(f, "%f ", parts->p[k]);
  }
  fprintf(f, "%f\n", parts->p[3*parts->n-1]);
  fprintf(f, "%f %f %f %f %f\n", pauli->qo, pauli->po, pauli->D, pauli->scut2, pauli->shift);
  fprintf(f, "%f %f %f %f %f %f %f %f %f\n", nuc->Vo, nuc->r1, nuc->r2, nuc->p1, nuc->p2, nuc->d, nuc->a, nuc->rcut, nuc->shift);
  fprintf(f, "%f %f %f %f\n", params->L, params->T, params->delta_q, params->delta_p);
  fclose(f);
  return 0;
}

int load_checkpoint(char *filename, struct Particles *parts, struct Pauli *pauli, struct Nuclear *nuc, struct Externos *params){
  FILE *f = fopen(filename, "r");
  int i;
  i = fscanf(f, "%d %f %f %f %f\n", &parts->n, &parts->mass, &parts->kinetic, &parts->pot_nuc, &parts->pot_pauli);
  for(int k = 0; k < 3*parts->n; k++){
    i = fscanf(f, "%f", &parts->q[k]);
  }
  for(int k = 0; k < 3*parts->n; k++){
    i = fscanf(f, "%f", &parts->p[k]);
  }
  i = fscanf(f, "%f %f %f %f %f\n", &pauli->qo, &pauli->po, &pauli->D, &pauli->scut2, &pauli->shift);
  i = fscanf(f, "%f %f %f %f %f %f %f %f %f\n", &nuc->Vo, &nuc->r1, &nuc->r2, &nuc->p1, &nuc->p2, &nuc->d, &nuc->a, &nuc->rcut, &nuc->shift);
  i = fscanf(f, "%f %f %f %f\n", &params->L, &params->T, &params->delta_q, &params->delta_p);
  fclose(f);
  return 0;
}

// --------------------------------- MAIN ----------------------------------- //


int main(){

// Particulas
  struct Particles parts;
  int N = 18;
  parts.n = N*N*N;
  parts.mass = 1.043916; // Masa protón, MeV*(10^-22 s/fm)^2
  parts.q = (float *) malloc(3*parts.n*sizeof(float));
  parts.p = (float *) malloc(3*parts.n*sizeof(float));
  for(int j = 0; j < 3*parts.n; j++){
    parts.q[j] = 0;
    parts.p[j] = 0;
  }

// Potencial de Pauli
  struct Pauli pauli;
  float h_barra = 6.582119;
  pauli.qo = 6; // fm
  pauli.po = 2.067; // MeV*10^-22 s/fm
  pauli.D = 34.32*pow(h_barra/(pauli.po*pauli.qo), 3) * 1; // MeV, mantengo el D* de Pauli
  pauli.scut2 = 6; // Un ~5% del máximo
  pauli.shift = pauli.D*exp(-0.5*pauli.scut2);

// Potencial nuclear
  struct Nuclear nuc;
  nuc.Vo = 25.93; // MeV
  nuc.r1 = 1.757; // fm
  nuc.r2 = 1.771; // fm
  nuc.p1 = 6.2;
  nuc.p2 = 3.0;
  nuc.d = 3.35; // fm
  nuc.a = 5.0/6.0; // fm
  nuc.rcut = 6; // fm; Un ~0.45% del máximo
  nuc.shift = 0;
  nuc.shift = interaction_nuc(nuc.rcut, &nuc);

// Parametros
  struct Externos params;
  float rho = 0.2;
  params.L = N/pow(rho, 1.0/3.0); // fm ; mayor a 2*qo*scut
  params.T = 0.5; // MeV
  params.delta_q = pauli.qo/2; // fm
  params.delta_p = pauli.po/2; // MeV*10^-22 s/fm

  char filename[255];
  // Distribucion de muchos T; muchas realizaciones
  clock_t start, end;
  double time;
  int pasos = 1000*parts.n;
  srand(1);

  float rhos[6] = {0.1, 0.125, 0.15, 0.175, 0.2, 0.225};
  //float rhos[2] = {0.01, 0.05};
/*
  params.L = N/pow(0.2, 1.0/3.0);

  set_box(&parts, params.L);
  set_p(&parts, params.T);

  params.delta_q = pauli.qo*params.L/1500;
  params.delta_p = pauli.po/50;

  energia(&parts, &pauli, &nuc, params.L);
  int aceptados = N_steps(&parts, &pauli, &nuc, &params, 10*parts.n);
  printf("%f %% \n", 10.0*aceptados/(parts.n));
  printf("%f + %f + %f = %f \n", parts.kinetic, parts.pot_nuc, parts.pot_pauli, parts.kinetic+parts.pot_nuc+parts.pot_pauli);
  energia(&parts, &pauli, &nuc, params.L);
  printf("%f + %f + %f = %f \n", parts.kinetic, parts.pot_nuc, parts.pot_pauli, parts.kinetic+parts.pot_nuc+parts.pot_pauli);

*/


  for (int k = 3; k < 6; k++){
    sprintf(filename, "x1/checkpoint_%f_18.txt", rhos[k]);
    load_checkpoint(filename, &parts, &pauli, &nuc, &params);

    pauli.D = 34.32*pow(h_barra/(pauli.po*pauli.qo), 3) * 1; // MeV, mantengo el D* de Pauli

    params.L = N/pow(rhos[k], 1.0/3.0);
    pauli.scut2 = 6; // Un ~5% del máximo
    pauli.shift = pauli.D*exp(-0.5*pauli.scut2);

    energia(&parts, &pauli, &nuc, params.L);

    params.delta_q = pauli.qo*params.L/3000;
    params.delta_p = pauli.po/500;

    start = clock();
    sprintf(filename, "x1/energias_%f.txt", rhos[k]);
    int aceptados = muestrear_energias(filename, &parts, &pauli, &nuc, &params, 1000, 1, 0);
    end = clock();
    time = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Muestreo rho = %f en %f segundos con %d aceptados (%f)\n", rhos[k], time, aceptados, ((float) aceptados)/pasos);
    sprintf(filename, "x1/checkpoint_%f_18.txt", rhos[k]);
    energia(&parts, &pauli, &nuc, params.L);
    save_checkpoint(filename, &parts, &pauli, &nuc, &params);

    printf("%f + %f + %f = %f \n", parts.kinetic, parts.pot_nuc, parts.pot_pauli, parts.kinetic+parts.pot_nuc+parts.pot_pauli);

    sprintf(filename, "x1/distribucion_%f.txt", rhos[k]);
    aceptados = muestrear_impulsos(filename, &parts, &pauli, &nuc, &params, 1, 0, 0);
  }

  free(parts.q);
  free(parts.p);
  return 0;
}
