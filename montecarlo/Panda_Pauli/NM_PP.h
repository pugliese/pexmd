#ifndef NM_PP_H
#define NM_PP_H

#include "math.h"

struct Celdas {
  int nc;
  float* q;
  float* p;
  float mass;
  int* siguiente;
  int* anterior;
  int* primero;
  int* celda;
  int M;
  float l;
};

struct Particles {
  int n;
  struct Celdas** tipos; //p up, p down, n up, n down
  float energy_pauli;
  float energy_panda_np;
  float energy_panda_nn;
  float kinetic;
};

struct Pauli {
  float qo;
  float po;
  float D;
  float scut2;
  float shift;
};

struct Panda_np {
  float mu_r;
  float mu_a;
  float V_r;
  float V_a;
  float rcut2;
  float shift;
};

struct Panda_nn {
  float mu_o;
  float V_o;
  float rcut2;
  float shift;
};

struct Externos {
  float L;
  float T;
  float delta_q;
  float delta_p;
};




#endif
