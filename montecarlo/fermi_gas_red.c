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
  parts->primero = (int *) malloc(M3*sizeof(int));
  parts->siguiente = (int *) malloc(N*sizeof(int));
  parts->anterior = (int *) malloc(N*sizeof(int));
  parts->celda = (int *) malloc(N*sizeof(int));
  for(int i = 0; i < M3; i++){
    parts->primero[i] = -1;
  }
  for(int i = 0; i < N; i++){
    idx = 0;
    for (int k = 0; k < 3; k++){
      idx = idx*parts->M;
      idx += (int) floor(parts->q[3*i+k]/parts->l);
      parts->q[3*i+k] -= parts->l*floor(parts->q[3*i+k]/parts->l);
    }
    parts->celda[i] = idx;
    parts->siguiente[i] = parts->primero[idx];    // Su siguiente es el primero previo (o -1 si no habia)
    parts->anterior[i] = -1;              // Es el nuevo primero; no tiene anterior
    if (parts->primero[idx] != -1){
      parts->anterior[parts->primero[idx]] = i; // Si habia primero previo, ahora tiene como anterior a i
    }
    parts->primero[idx] = i;         // i es el nuevo primero
  }
  return 0;
}

int actualizar_lista(struct Particles *parts, int i, int idx){
  int N = parts->n;
  if (idx!=parts->celda[i]){
    // Remuevo i de su celda anterior
    int sig_i = parts->siguiente[i];
    int ant_i = parts->anterior[i];
    if (sig_i != -1){
      parts->anterior[sig_i] = ant_i;
    }
    if (ant_i != -1){
      parts->siguiente[ant_i] = sig_i;
    }else{
      parts->primero[parts->celda[i]] = sig_i;
    }
    // Lo pongo en su nueva celda
    parts->celda[i] = idx;
    parts->siguiente[i] = parts->primero[idx];  // Su siguiente es el primero previo (o -1 si no habia)
    parts->anterior[i] = -1;            // Es el nuevo primero; no tiene anterior
    if (parts->primero[idx] != -1){
      parts->anterior[parts->primero[idx]] = i; // Si habia primero previo, ahora tiene como anterior a i
    }
    parts->primero[idx] = i;       // i es el nuevo primero
  }
  return idx;
}

int print_lista(struct Particles *parts){
  int j, M3 = parts->M*parts->M*parts->M, flag;
  for(int i = 0; i < M3; i++){
    j = parts->primero[i];
    printf("Celda %d:\n%d", i, j);
    while(parts->siguiente[j] != -1){
      printf(" -> %d", parts->siguiente[j]);
      j = parts->siguiente[j];
    }
    printf("\n%d", j);
    while(parts->anterior[j] != -1){
      printf(" <- %d", parts->anterior[j]);
      j = parts->anterior[j];
    }
    printf("\n", j);
  }
}

// -------------------------------- ENERGIA --------------------------------- //


float distancia_fases(float *q1, float *q2, float *p1, float *p2, int *delta_idx, float l, struct Pauli *pauli){
  double dq2 = 0;
  double dp2 = 0;
  for(int k = 0; k < 3; k++){
    dq2 = dq2 + (q1[k] - q2[k] + delta_idx[k]*l)*(q1[k] - q2[k] + delta_idx[k]*l);
    dp2 = dp2 + (p1[k] - p2[k])*(p1[k] - p2[k]);
  }
  dq2 = dq2/(pauli->qo*pauli->qo);
  dp2 = dp2/(pauli->po*pauli->po);
  float s2 = dq2 + dp2;
  return s2;
}

float presion_temperatura(struct Particles *parts, struct Pauli *pauli, struct Externos *params, float *PyT){
  double qF = 0, qp = 0;
  float qo2 = pauli->qo*pauli->qo;
  float po2 = pauli->po*pauli->po;
  float dq2, dp2, s2;
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
  float delta_q[3], delta_p[3], force_mod, gorce_mod;
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
          dq2 = 0, dp2 = 0;
          for(int l = 0; l < 3; l++){
            delta_q[l] = parts->q[3*i+l] - parts->q[3*j+l] + idxs[3*k+l]*parts->l;
            delta_p[l] = parts->p[3*i+l] - parts->p[3*j+l];
            dq2 = dq2 + delta_q[l]*delta_q[l];
            dp2 = dp2 + delta_p[l]*delta_p[l];
          }
          dq2 = dq2/(pauli->qo*pauli->qo);
          dq2 = dq2*dq2*dq2;
          dp2 = dp2/(pauli->po*pauli->po);
          s2 = dq2 + dp2;
          if (s2 < pauli->scut2) {
            //force_mod = 6*pauli->D*(1/dq2-1)/dq2;
            force_mod = pauli->D*exp(-s2/2)/qo2;
            gorce_mod = force_mod*qo2/po2;
            for (int l = 0; l < 3; l++){
              qF += force_mod*delta_q[l]*delta_q[l];
              qp -= gorce_mod*delta_p[l]*delta_p[l];
            }
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

float interaction(float *q1, float *q2, float *p1, float *p2, int *delta_idx, float l, struct Pauli *pauli){
  float pot = 0;
  float s2 = distancia_fases(q1, q2, p1, p2, delta_idx, l, pauli);
  if (s2 < pauli->scut2){
    pot = pauli->D*exp(-0.5*s2) - pauli->shift;
  }
  return pot;
}

int energia(struct Particles *parts, struct Pauli *pauli){
  float kin = 0, pot = 0, f;
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
    i = parts->primero[m];
    for (int k = 0; k < 14; k++){
    //for (int k = 0; k < 27; k++){
      idx = (((mx+idxs[3*k]+M) % M)*M + (my+idxs[3*k+1]+M) % M)*M + (mz+idxs[3*k+2]+M) % M;
      i = parts->primero[m];
      while (i != -1){
        j = parts->primero[idx];
        while (j !=-1 && j != i){
          pot = pot + interaction(parts->q+3*i, parts->q+3*j, parts->p+3*i, parts->p+3*j, idxs+3*k, parts->l, pauli);
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
  parts->energy = parts->kinetic + pot;
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

float delta_energia_pot(struct Particles *parts, struct Pauli *pauli, float *new_q, float *new_p, int i, int *new_ms, int *ms){
  double new_pot = 0, old_pot = 0;
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
        old_pot = old_pot + interaction(parts->q+3*i, parts->q+3*j, parts->p+3*i, parts->p+3*j, idxs, parts->l, pauli);
      }
      j = parts->siguiente[j];
    }
    new_idx = (((new_ms[0]+idxs[0]+M) % M)*M + (new_ms[1]+idxs[1]+M) % M)*M + (new_ms[2]+idxs[2]+M) % M;
    j = parts->primero[new_idx];
    while (j != -1){
      if (j != i){
        new_pot = new_pot + interaction(new_q, parts->q+3*j, new_p, parts->p+3*j, idxs, parts->l, pauli);
      }
      j = parts->siguiente[j];
    }
  }
  return new_pot - old_pot;
}

// ------------------------------- MUESTREO --------------------------------- //

int step(struct Particles *parts, struct Pauli *pauli, struct Externos *params){
  int i = elegir_part(parts->n);
  int m = parts->celda[i], res = 1, M = parts->M;
  int ms[3], new_m, new_ms[3], idxs[3];
  ms[0] = m/(M*M);
  ms[1] = (m/M) % M;
  ms[2] = m % M;
  float new_q[3], new_p[3];
  for(int k = 0; k < 3; k++){
    new_q[k] = parts->q[3*i+k] + (2*uniform()-1)*params->delta_q;
    idxs[k] = (new_q[k]>parts->l) - (new_q[k]<0);
    new_q[k] = new_q[k] - parts->l*idxs[k];
    new_p[k] = parts->p[3*i+k] + (2*uniform()-1)*params->delta_p;
    new_ms[k] = (ms[k]+idxs[k]+M) % M;
  }
  float delta_kin = delta_energia_kin(parts, new_p, i);
  float delta_E = delta_kin + delta_energia_pot(parts, pauli, new_q, new_p, i, new_ms, ms);
  if (delta_E <= 0){
    new_m = (new_ms[0]*M + new_ms[1])*M + new_ms[2];
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
      new_m = (new_ms[0]*M + new_ms[1])*M + new_ms[2];
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

float set_box(struct Particles *parts, float rcut, float L){
  int n_lado = pow(parts->n, 1.0/3.0), i;
  float dL = L/n_lado;
  for(int x = 0; x < n_lado; x++){
    for(int y = 0; y < n_lado; y++){
      for(int z = 0; z < n_lado; z++){
        i = x*n_lado*n_lado + y*n_lado + z;
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

int N_steps(struct Particles *parts, struct Pauli *pauli, struct Externos *params, int Nsteps){
  int aceptados = 0;
  for (int k = 0; k < Nsteps; k++) {
    //printf("Progreso: %d%%\r", k*100/Nsteps);
    aceptados = aceptados + step(parts, pauli, params);
  }
  return aceptados;
}

int muestrear_termo(char *filename, struct Particles *parts, struct Pauli *pauli, struct Externos *params, int Nsamp, int factor, int factor_term){
  // Termalizacion
  int segs = time(NULL);
  int aceptados = N_steps(parts, pauli, params, factor_term*parts->n);
  // Muestreo
  float ekin = 0, etot = 0, PyT_tot[2], PyT[2];
  PyT_tot[0] = 0, PyT_tot[1] = 0;
  FILE *f = fopen(filename, "a");
  for(int k = 0; k < Nsamp; k++){
    aceptados += N_steps(parts, pauli, params, factor*parts->n);
    ekin += parts->kinetic/parts->n;
    etot += parts->energy/parts->n;
    presion_temperatura(parts, pauli, params, PyT);
    PyT_tot[0] += PyT[0];
    PyT_tot[1] += PyT[1];
  }
  ekin = ekin/Nsamp;
  etot = etot/Nsamp;
  PyT_tot[0] = PyT_tot[0]/Nsamp;
  PyT_tot[1] = PyT_tot[1]/Nsamp;
  segs = time(NULL)-segs;
  fprintf(f, "0 %f %f %f %f %f %f %d\n", params->T, ekin, etot-ekin, etot, PyT_tot[0], PyT_tot[1], segs);
  fclose(f);
  return aceptados;
}

int muestrear_todo(char *filename_termo, char *filename_dist, struct Particles *parts, struct Pauli *pauli, struct Externos *params, int Nsamp, int factor, int factor_term){
  // Termalizacion
  int segs = time(NULL);
  int aceptados = N_steps(parts, pauli, params, factor_term*parts->n);
  // Muestreo
  float ekin_i = 0, ekin = 0, etot = 0, PyT_tot[2], PyT[2];
  PyT_tot[0] = 0, PyT_tot[1] = 0;
  float var_ekin = 0, var_etot = 0;
  FILE *f = fopen(filename_dist, "w");
  for(int k = 0; k < Nsamp; k++){
    aceptados += N_steps(parts, pauli, params, factor*parts->n);
    ekin += parts->kinetic/parts->n;
    etot += parts->energy/parts->n;
    var_ekin += parts->kinetic*parts->kinetic/(parts->n*parts->n);
    var_etot += parts->energy*parts->energy/(parts->n*parts->n);
    presion_temperatura(parts, pauli, params, PyT);
    PyT_tot[0] += PyT[0];
    PyT_tot[1] += PyT[1];
    for(int i = 0; i < parts->n; i++){
      ekin_i = 0;
      for(int l = 0; l < 3; l++){
        ekin_i += parts->p[3*i+l]*parts->p[3*i+l];
      }
      ekin_i = 0.5*ekin_i/parts->mass;
      fprintf(f, "%f ", ekin_i);
    }
  }
  fprintf(f, "\n");
  fclose(f);
  ekin = ekin/Nsamp;
  etot = etot/Nsamp;
  var_ekin = var_ekin/Nsamp - ekin*ekin;
  var_etot = var_etot/Nsamp - etot*etot;
  PyT_tot[0] = PyT_tot[0]/Nsamp;
  PyT_tot[1] = PyT_tot[1]/Nsamp;
  segs = time(NULL)-segs;
  printf("%f %f\n", var_ekin, var_etot);
  f = fopen(filename_termo, "a");
  fprintf(f, "0 %f %f %f %f %f %f %d\n", params->T, ekin, etot-ekin, etot, PyT_tot[0], PyT_tot[1], segs);
  fclose(f);
  return aceptados;
}


// --------------------------------- MAIN ----------------------------------- //


int main(int argc, char *argv[]){

  float rho = 0.2;
  float D = 100;
  if (argc >= 2){
    int i = sscanf(argv[1], "%f\n", &rho);
    if (argc >= 3){
      int i = sscanf(argv[2], "%f\n", &D);
    }
  }

// Particulas
  struct Particles parts;
  int N = 20;
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
  parts.mass = 1;
  pauli.qo = 1;
  pauli.po = 1;
  pauli.D = D;
  pauli.scut2 = 10;
  pauli.shift = pauli.D*exp(-0.5*pauli.scut2);

// Parametros
  struct Externos params;
  params.L = N*pauli.qo/pow(rho, 1.0/3); // mayor a 2*qo*scut
  params.T = 30; // MeV
  float delta_qo = 1+0.5*log(100/D)/log(10);
  float delta_po = 1+0.25*log(100/D)/log(10);
  params.delta_q = delta_qo*pow(params.T/2, 0.8);
  params.delta_p = delta_po*pow(params.T/2, 0.8);


  set_box(&parts, sqrt(pauli.scut2)*pauli.qo, params.L);
  float Tcin = set_p(&parts, params.T);
  float ekin_i, ekin_max = 0;
  energia(&parts, &pauli);
  int M3 = parts.M*parts.M*parts.M;
  //printf("%f  %f  %f\n", Tcin, parts.kinetic, parts.energy);
  printf("Sistema con %dx%dx%d = %d celdas\n", parts.M, parts.M, parts.M, M3);

  float dT = 0.25;
  int segs, aceptados, Nsamp = 100, factor = 10, factor_term = 1000;
  char filename_termo[255], filename_dist[255];
  sprintf(filename_termo, "prueba_pauli/termo_pauli_%1.2f_%.2f.txt", rho, D);
  FILE *f = fopen(filename_termo, "w");
  fprintf(f, "# T_MC\t\tEcin\t\tEpot\t\tEtot\t\tP\t\tT_V\t\tTiempo\n");
  fclose(f);
  N_steps(&parts, &pauli, &params, factor*parts.n);
  for (int k = 0; k < 120; k++){
    params.delta_q = delta_qo*pow(params.T/2, 0.8);
    params.delta_p = delta_po*pow(params.T/2, 0.8);
    segs = time(NULL);
    if (k%10 == 0){
      sprintf(filename_dist, "prueba_pauli/dist_%.2f_%.2f_%.2f.txt", rho, D, params.T);
      printf("Guardando distribucion de energia cinetica: %s\n", filename_dist);
      aceptados = muestrear_todo(filename_termo, filename_dist, &parts, &pauli, &params, Nsamp, factor, factor_term);
    }else{
      aceptados = muestrear_termo(filename_termo, &parts, &pauli, &params, factor, Nsamp, factor_term);
    }
    segs = time(NULL) - segs;
    printf("T = %1.2f con %1.2f%% aceptados en %d segundos\n", params.T, 100*((float) aceptados)/((factor_term+factor*Nsamp)*parts.n), segs);
    params.T -= dT;
  }

  free(parts.q);
  free(parts.p);
  free(parts.celda);
  free(parts.primero);
  free(parts.siguiente);
  free(parts.anterior);

  return 0;
}
