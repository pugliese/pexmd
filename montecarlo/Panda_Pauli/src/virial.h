#ifndef VIRIAL_H
#define VIRIAL_H

float force_mod_panda_nn(float r, struct Panda_nn *panda_nn);
float force_mod_panda_np(float r, struct Panda_np *panda_np);
float fgorce_mod_pauli(float rsq, float *p1, float *p2, struct Pauli *pauli, float *gorce_mod);
float fgorce_mod(int t1, int t2, float rsq, float *p1, float *p2, struct Interaction *pot_tot, float *gorce_mod);
float presion_temperatura(struct Particles *parts, struct Interaction *pot_tot, struct Externos *params, float *PyT);

#endif
