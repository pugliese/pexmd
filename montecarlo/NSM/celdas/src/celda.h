#ifndef C_H
#define C_H

#include "math.h"

struct Particles {
  int n;
  int* type;  // 0: prot up | 1: prot down | 2: neut up | 3: neut down
  float* q;
  float* p;
  float mass;

  int* siguiente;
  int* anterior;
  int* primero;
  int* celda;
  int M;
  float l;

  float energy_pauli;
  float energy_panda_nn;
  float energy_panda_np;
  float energy_qcnm;
  float energy_coul;
  float kinetic;

  double *forces;
  double *gorces;
};

int armar_lista(struct Particles *parts, float rcut, float L);
int actualizar_lista(struct Particles *parts, int i, int idx);
int corregir_celda(float *q, float l);
int ajustar_celdas(struct Particles *parts);
int print_lista(struct Particles *parts);
int save_lammpstrj(char *filename, struct Particles *parts, float L, int append);
int load_lammpstrj(char *filename, struct Particles *parts, float* L, float rcut);
int load_lammpstrj_guille(char *filename, struct Particles *parts, float* L, float rcut);

#endif
