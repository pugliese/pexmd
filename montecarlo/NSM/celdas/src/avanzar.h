#ifndef AV_H
#define AV_H

struct Externos {
  float L;
  float T;
  float delta_q;
  float delta_p;
};

float delta_energia_kin(struct Particles *parts, float *new_p, int i);

float delta_energia_pot(struct Particles *parts, struct Interaction *pot_tot, float *new_q, float *new_p, int i, int *new_ms, int *ms);
float delta_energia_trap(struct Particles *parts, float *new_q, int i, int *new_ms, int *ms, float L, struct Trap *trap);
float delta_energia_pot_sin_LUT(struct Particles *parts, struct Interaction *pot_tot, float *new_q, float *new_p, int i, int *new_ms, int *ms);
float delta_energia_trap_sin_LUT(struct Particles *parts, float *new_q, int i, int *new_ms, int *ms, float L, struct Trap *trap);

int step(struct Particles *parts, struct Interaction *pot_tot, struct Externos *params);
int N_steps(struct Particles *parts, struct Interaction *pot_tot, struct Externos *params, int Nsteps);

int step_sin_LUT(struct Particles *parts, struct Interaction *pot_tot, struct Externos *params);
int N_steps_sin_LUT(struct Particles *parts, struct Interaction *pot_tot, struct Externos *params, int Nsteps);

#endif
