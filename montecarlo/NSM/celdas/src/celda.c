#include "celda.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

int armar_lista(struct Particles *parts, float rcut, float L){
  int idx, N = parts->n;
  float c;
  parts->M = (int) floor(L/rcut);
  parts->l = L/parts->M;
  int M3 = parts->M*parts->M*parts->M;
  printf("Sistema de %dx%dx%d = %d celdas de lado %f\n", parts->M, parts->M, parts->M, M3, parts->l);
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
      c = parts->q[3*i+k]/parts->l;
      idx += (int) floor(c);
      parts->q[3*i+k] -= parts->l*floor(c);
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

int corregir_celda(float *q, float l){
  int correccion = (int) floor(q[0]/l);
  q[0] -= l*correccion;
	return correccion;
}

int ajustar_celdas(struct Particles *parts){
	float l = parts->l;
	int M = parts->M, idx, new_m, ms[3], i, m = 0;
	for(ms[0] = 0; ms[0] < M; ms[0]++){
		for(ms[1] = 0; ms[1] < M; ms[1]++){
			for(ms[2] = 0; ms[2] < M; ms[2]++){
				i = parts->primero[m];
				while(i != -1){
					new_m = 0;
					for(int k = 0; k < 3; k++){
						idx = corregir_celda(parts->q+i, l);
						new_m = new_m*M + (ms[k] + idx + M) % M;
					}
					actualizar_lista(parts, i, new_m);
					i = parts->siguiente[i];
				}
				m++;
			}
		}
	}
	return 0;
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

int load_lammpstrj(char *filename, struct Particles *parts, float* L, float rcut){
  FILE *f = fopen(filename, "r");
  char buffer[255];
  int id;
  for(int l = 0; l < 3; l++){
    fgets(buffer, 255, f);
  }
  id = fscanf(f, "%d\n", &parts->n);
  for(int l = 0; l < 2; l++){
    fgets(buffer, 255, f);
  }
  id = fscanf(f, "0 %f\n", L);
  for(int l = 0; l < 2; l++){
    fgets(buffer, 255, f);
  }
  parts->type = (int *) malloc(parts->n*sizeof(int));
  parts->q = (float *) malloc(3*parts->n*sizeof(float));
  parts->p = (float *) malloc(3*parts->n*sizeof(float));
	for(int l = 0; l < parts->n; l++){
		id = fscanf(f, "%d %d %f %f %f %f %f %f\n", &id, parts->type+l, parts->q+3*l, parts->q+3*l+1, parts->q+3*l+2, parts->p+3*l, parts->p+3*l+1, parts->p+3*l+2);
	}
  fclose(f);
  armar_lista(parts, rcut, *L);
  return 0;
}
