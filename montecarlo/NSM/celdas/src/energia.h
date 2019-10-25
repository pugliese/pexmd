#ifndef E_H
#define E_H

#include "math.h"

int construir_lista_celdas(int K_c, int *idxs);

int energia(struct Particles *parts, struct Interaction *pot_tot);
int energia_sin_LUT(struct Particles *parts, struct Interaction *pot_tot);

float fgorces_and_energy(struct Particles *parts, struct Interaction *pot_tot);
float fgorces_and_energy_sin_LUT(struct Particles *parts, struct Interaction *pot_tot);

#endif
