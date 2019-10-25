#include "general.h"
#include "celda.h"
#include "interaccion.h"
#include "inicializar.h"
#include "energia.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

float set_box(struct Particles *parts, float rcut, float L){
  int n_lado = 1, i = 0;
  while (n_lado*n_lado*n_lado < parts->n) n_lado++;
  float dL = L/n_lado;
  for(int x = 0; x < n_lado; x++){
    for(int y = 0; y < n_lado; y++){
      for(int z = 0; z < n_lado; z++){
        parts->q[3*i] = dL*(0.5 + x);
        parts->q[3*i+1] = dL*(0.5 + y);
        parts->q[3*i+2] = dL*(0.5 + z);
        if ((x+y+z)%2==0){
          parts->type[i] = (z%2);
        }else{
          parts->type[i] = 3 - (y%2);
        }
        i++;
      }
    }
  }
  armar_lista(parts, rcut, L);
  return dL;
}

float redondear_SC(struct Particles *parts, float rcut, float L){
  int n_lado = 1;
  while (n_lado*n_lado*n_lado < parts->n) n_lado++;
  int idx = min_vec(parts->q, parts->n);
  float dL = L/n_lado;
  float q0[3];
  for(int k = 0; k < 3; k++)  q0[k] = parts->q[3*idx+k];// - 0.5*dL;
  for(int i = 0; i < parts->n; i++){
    for(int k = 0; k < 3; k++){
      parts->q[3*i+k] = dL*(0.5+((int)round((parts->q[3*i+k]-q0[k])/dL)%n_lado));
    }
  }
  armar_lista(parts, rcut, L);
  return dL;
}

float set_box_fund_pauli(struct Particles *parts, float rcut, float L){
  int n_lado = 1, i = 0;
  while (n_lado*n_lado*n_lado < parts->n) n_lado++;
  float dL = L/n_lado;
  for(int x = 0; x < n_lado; x++){
    for(int y = 0; y < n_lado; y++){
      for(int z = 0; z < n_lado; z++){
        parts->q[3*i] = dL*(0.5 + x);
        parts->q[3*i+1] = dL*(0.5 + y);
        parts->q[3*i+2] = dL*(0.5 + z);
        if ((x+y+z)%2==0){
          parts->type[i] = (z%2);
        }else{
          parts->type[i] = 3-(y%2);
        }
        i++;
        if(i==parts->n){
          x = n_lado;
          y = n_lado;
          z = n_lado;
        }
      }
    }
  }
  armar_lista(parts, rcut, L);
  return dL;
}

float set_box_B2(struct Particles *parts, float rcut, float L){
  int n_lado = 1, i = 0;
  while (n_lado*n_lado*n_lado < parts->n/2) n_lado++;
  float dL = L/n_lado;
  for(int x = 0; x < n_lado; x++){
    for(int y = 0; y < n_lado; y++){
      for(int z = 0; z < n_lado; z++){
        parts->q[6*i+0] = dL*(0.25 + x);
        parts->q[6*i+1] = dL*(0.25 + y);
        parts->q[6*i+2] = dL*(0.25 + z);
        parts->type[2*i+0] = (x+y+z)%2;
        parts->q[6*i+3] = dL*(0.75 + x);
        parts->q[6*i+4] = dL*(0.75 + y);
        parts->q[6*i+5] = dL*(0.75 + z);
        parts->type[2*i+1] = 2 + (x+y+z)%2;
        i++;
      }
    }
  }
  armar_lista(parts, rcut, L);
  return dL/2;
}

float set_box_B3(struct Particles *parts, float rcut, float L){
  int n_lado = 1, i = 0;
  while (8*n_lado*n_lado*n_lado < parts->n) n_lado++;
  float dL = L/n_lado;
  for(int x = 0; x < n_lado; x++){
    for(int y = 0; y < n_lado; y++){
      for(int z = 0; z < n_lado; z++){
        parts->q[24*i+0] = dL*(0.125 + x);
        parts->q[24*i+1] = dL*(0.125 + y);
        parts->q[24*i+2] = dL*(0.125 + z);
        parts->type[8*i+0] = (x+y+z)%2;
        parts->q[24*i+3] = dL*(0.125 + x + 0.5);
        parts->q[24*i+4] = dL*(0.125 + y + 0.5);
        parts->q[24*i+5] = dL*(0.125 + z);
        parts->type[8*i+1] = (x+y+z+1)%2;
        parts->q[24*i+6] = dL*(0.125 + x + 0.5);
        parts->q[24*i+7] = dL*(0.125 + y);
        parts->q[24*i+8] = dL*(0.125 + z + 0.5);
        parts->type[8*i+2] = (x+y+z+1)%2;
        parts->q[24*i+9] = dL*(0.125 + x);
        parts->q[24*i+10] = dL*(0.125 + y + 0.5);
        parts->q[24*i+11] = dL*(0.125 + z + 0.5);
        parts->type[8*i+3] = (x+y+z+1)%2;

        parts->q[24*i+12] = dL*(0.375 + x);
        parts->q[24*i+13] = dL*(0.375 + y);
        parts->q[24*i+14] = dL*(0.375 + z);
        parts->type[8*i+4] = 2 + (x+y+z)%2;
        parts->q[24*i+15] = dL*(0.375 + x + 0.5);
        parts->q[24*i+16] = dL*(0.375 + y + 0.5);
        parts->q[24*i+17] = dL*(0.375 + z);
        parts->type[8*i+5] = 2 + (x+y+z+1)%2;
        parts->q[24*i+18] = dL*(0.375 + x + 0.5);
        parts->q[24*i+19] = dL*(0.375 + y);
        parts->q[24*i+20] = dL*(0.375 + z + 0.5);
        parts->type[8*i+6] = 2 + (x+y+z+1)%2;
        parts->q[24*i+21] = dL*(0.375 + x);
        parts->q[24*i+22] = dL*(0.375 + y + 0.5);
        parts->q[24*i+23] = dL*(0.375 + z + 0.5);
        parts->type[8*i+7] = 2 + (x+y+z+1)%2;
        i++;
      }
    }
  }
  armar_lista(parts, rcut, L);
  return dL/8;
}

int check_lattice(struct Particles *parts, float* as, float tol){
  float a1[3], a2[3], a3[3], rp[3], rn[3];
  int res = 1, ip = 0, in = 0;
  while ((parts->type[in] / 2) == 0){
    ip = in;
    in++;
  }
  for(int k = 0; k < 3; k++){
    a1[k] = as[3*k+0];
    a2[k] = as[3*k+1];
    a3[k] = as[3*k+2];
    rp[k] = parts->q[3*ip+3];
    rn[k] = parts->q[3*in+3];
  }
  float a1_x_a2[3], a2_x_a3[3], a3_x_a1[3], r[3], m1, m2, m3;
  producto_vectorial(a1, a2, a1_x_a2);
  producto_vectorial(a2, a3, a2_x_a3);
  producto_vectorial(a3, a1, a3_x_a1);
  float A = producto_interno(a1, a2_x_a3);
  printf("A = %f\n", A);
  printf("   a1       a2       a3      a1xa2    a2xa3     a3xa1     rp       rn\n");
  for(int k = 0; k < 3; k++){
    printf("%f %f %f %f %f %f %f %f\n", a1[k], a2[k], a3[k], a1_x_a2[k], a2_x_a3[k], a3_x_a1[k], rp[k], rn[k]);
  }
  for(int i = 0; i < parts->n; i++){
    if(parts->type[i]/2 == 0){ // Proton
      for(int k = 0; k < 3; k++){
        r[k] = parts->q[3*i+k] - rp[k];
      }
    }else{
      for(int k = 0; k < 3; k++){
        r[k] = parts->q[3*i+k] - rn[k];
      }
    }
    m1 = producto_interno(r, a2_x_a3)/A;
    res = res && (fabs(round(m1) - m1) < tol);
    m2 = producto_interno(r, a3_x_a1)/A;
    res = res && (fabs(round(m2) - m2) < tol);
    m3 = producto_interno(r, a1_x_a2)/A;
    res = res && (fabs(round(m3) - m3) < tol);
    if (!res){
      printf("%f %f %f\n", m1, m2, m3);
      break;
    }
  }
  return res;
}

float load_and_rep(char *filename, struct Particles *parts, float rcut, float *L, int K){
  struct Particles parts2;
  float L2;
  load_lammpstrj(filename, &parts2, &L2, rcut);
  if (L2 < rcut) printf("Fucking forget it\n");
  *L = L2*K;
  parts->n = K*K*K*parts2.n;
  int idx, c, cx, cy, cz, M = parts2.M;
  parts->type = (int *) malloc(parts->n*sizeof(int));
  parts->q = (float *) malloc(3*parts->n*sizeof(float));
  parts->p = (float *) malloc(3*parts->n*sizeof(float));
  for(int x = 0; x < K; x++){
    for(int y = 0; y < K; y++){
      for(int z = 0; z < K; z++){
        idx = (x+K*(y+K*z))*parts2.n;
        for(int i = 0; i < parts2.n; i++){
          c = parts2.celda[i];
          cx = c/(M*M);
          cy = (c/M) % M;
          cz = c % M;
          parts->q[3*(i+idx)] = parts2.q[3*i]+x*L2 + cx*parts2.l;
          parts->q[3*(i+idx)+1] = parts2.q[3*i+1]+y*L2 + cy*parts2.l;
          parts->q[3*(i+idx)+2] = parts2.q[3*i+2]+z*L2 + cz*parts2.l;
          parts->p[3*(i+idx)] = parts2.p[3*i];
          parts->p[3*(i+idx)+1] = parts2.p[3*i+1];
          parts->p[3*(i+idx)+2] = parts2.p[3*i+2];
          parts->type[i+idx] = parts2.type[i];
        }
      }
    }
  }
  armar_lista(parts, rcut, *L);
  liberar(&parts2);
  return 0;
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

int inicializar(struct Particles *parts, int *comps, int n_types, float mass, float L, float T, struct Interaction *pot_tot, int md){
  parts->n = 0;
  for (int l = 0; l < n_types; l++){
    parts->n += comps[l];
  }
  parts->type = (int *) malloc(parts->n*sizeof(int));
  parts->q = (float *) malloc(3*parts->n*sizeof(float));
  parts->p = (float *) malloc(3*parts->n*sizeof(float));
  if(md){
    parts->forces = (double *) malloc(3*parts->n*sizeof(float));
    parts->gorces = (double *) malloc(3*parts->n*sizeof(float));
  }else{
    parts->forces = (double *) malloc(sizeof(float));
    parts->gorces = (double *) malloc(sizeof(float));
  }
  parts->mass = mass;
  set_box(parts, pot_tot->rcut, L);
  set_p(parts, T);
  energia(parts, pot_tot);
  return 0;
}

int inicializar_nucleo(struct Particles *parts, int *comps, int n_types, float mass, float L, float T, float R, struct Interaction *pot_tot, int md){
  parts->n = 0;
  for (int l = 0; l < n_types; l++){
    parts->n += comps[l];
  }
  parts->type = (int *) malloc(parts->n*sizeof(int));
  parts->q = (float *) malloc(3*parts->n*sizeof(float));
  parts->p = (float *) malloc(3*parts->n*sizeof(float));
  if(md){
    parts->forces = (double *) malloc(3*parts->n*sizeof(float));
    parts->gorces = (double *) malloc(3*parts->n*sizeof(float));
  }else{
    parts->forces = (double *) malloc(sizeof(float));
    parts->gorces = (double *) malloc(sizeof(float));
  }
  parts->mass = mass;
  int j = 0;
  for (int l = 0; l < n_types; l++){
    for(int i = 0; i < comps[l]; i++){
      parts->type[j] = l;
      for(int k = 0; k < 3; k++){
        parts->q[3*j+k] = 0.5*L + R*(2*uniform()-1);
      }
      j++;
    }
  }
  armar_lista(parts, pot_tot->rcut, L);
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
  free(parts->celda);
  free(parts->primero);

  free(parts->forces);
  free(parts->gorces);

  return 0;
}
