#ifndef AV_H
#define AV_H

struct Externos {
  float L;
  float T;
  float delta_q;
  float delta_p;
};

float delta_energia_kin(struct Particles *parts, float *new_p, int i);

float delta_energia_pot(struct Particles *parts, struct Interaction *pot_tot, float *new_q, float *new_p, int i, int *new_ms, int *ms, float* delta_pauli);
float delta_energia_pot_QCNM(struct Particles *parts, struct Interaction *pot_tot, float *new_q, float *new_p, int i, int *new_ms, int *ms, float* delta_pauli);
float delta_energia_pot_sin_LUT(struct Particles *parts, struct Interaction *pot_tot, float *new_q, float *new_p, int i, int *new_ms, int *ms, float* delta_pauli);

int step(struct Particles *parts, struct Interaction *pot_tot, struct Externos *params);
int N_steps(struct Particles *parts, struct Interaction *pot_tot, struct Externos *params, int Nsteps);
int step_QCNM(struct Particles *parts, struct Interaction *pot_tot, struct Externos *params);

int step_sin_LUT(struct Particles *parts, struct Interaction *pot_tot, struct Externos *params);
int N_steps_sin_LUT(struct Particles *parts, struct Interaction *pot_tot, struct Externos *params, int Nsteps);
int N_steps_QCNM(struct Particles *parts, struct Interaction *pot_tot, struct Externos *params, int Nsteps);

#endif
