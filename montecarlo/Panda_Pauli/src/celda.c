#include "celda.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>


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
      //if (parts->q[3*i+k]<parts->l*floor(parts->q[3*i+k]/parts->l)) printf("%f %f\n", parts->q[3*i+k], parts->l*floor(parts->q[3*i+k]/parts->l));
      parts->q[3*i+k] -= parts->l*floor(parts->q[3*i+k]/parts->l);
      //if (parts->q[3*i+k]>parts->l) printf("%d) Ouch1\n", 3*i+k);
      //if (parts->q[3*i+k]<0) printf("%d) Ouch2\n", k);
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
  int j, M3 = parts->M*parts->M*parts->M;
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
    printf("\n");
  }
  return 0;
}

int save_lammpstrj(char *filename, struct Particles *parts, float L, int append){
  FILE *f;
  if (append) f = fopen(filename, "a");
  else f = fopen(filename, "w");
	int m, mx, my, mz, M = parts->M;
	fprintf(f, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp\n", 0, parts->n);
	for(int l = 0; l < 3; l++){
		fprintf(f, "0 %f\n", L);
	}
	fprintf(f, "ITEM: ATOMS id type x y z vx vy vz \n");
	for(int l = 0; l < parts->n; l++){
		m = parts->celda[l];
		mx = m/(M*M);
		my = (m/M) % M;
		mz = m % M;
		fprintf(f, "%d %d %f %f %f %f %f %f\n", l, parts->type[l], mx*parts->l+parts->q[3*l], my*parts->l+parts->q[3*l+1], mz*parts->l+parts->q[3*l+2], parts->p[3*l], parts->p[3*l+1], parts->p[3*l+2]);
	}
  fclose(f);
  return 0;
}
