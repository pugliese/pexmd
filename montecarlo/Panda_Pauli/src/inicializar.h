#ifndef INIT_H
#define INIT_H

#include "math.h"

int energia(struct Particles *parts, struct Interaction *pot_tot);
int energia_sin_LUT(struct Particles *parts, struct Interaction *pot_tot);
int energia_QCNM(struct Particles *parts, struct Interaction *pot_tot);

float set_box(struct Particles *parts, float rcut, float L);
float set_box_fund_pauli(struct Particles *parts, float rcut, float L);
float set_box_B2(struct Particles *parts, float rcut, float L);
float set_box_B3(struct Particles *parts, float rcut, float L);
int check_lattice(struct Particles *parts, float* as, float tol);
float redondear_SC(struct Particles *parts, float rcut, float L);
float set_p(struct Particles *parts, float T);
int inicializar(struct Particles *parts, int *comps, int n_types, float mass, float L, float T, struct Interaction *pot_tot);
float load_and_rep(char *filename, struct Particles *parts, float rcut, float *L, int K);
int liberar(struct Particles *parts);

#endif
