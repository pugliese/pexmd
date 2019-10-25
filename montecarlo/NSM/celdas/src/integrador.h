#ifndef I_H
#define I_H

#include "math.h"

int verlet_first_step(struct Particles *parts, float dt);
int verlet_second_step(struct Particles *parts, float dt);
int velocity_verlet(struct Particles *parts, struct Interaction *pot_tot, float dt);

int RK2_half_step(struct Particles *parts, float dt);
int RK2(struct Particles *parts, struct Interaction *pot_tot, float dt);

int MidPoint_Rule_iteration(struct Particles *parts, float *qo, float *po, float dt);
int MidPoint_Rule(struct Particles *parts, struct Interaction *pot_tot, float dt, int K);

#endif
