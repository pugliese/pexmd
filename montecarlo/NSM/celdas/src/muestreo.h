#ifndef MUESTREO_H
#define MUESTREO_H

float presion_temperatura(struct Particles *parts, struct Interaction *pot_tot, float *PyT);
float presion_temperatura_sin_LUT(struct Particles *parts, struct Interaction *pot_tot, float *PyT);
float mean_radius(struct Particles *parts);

int muestrear_termo(char *filename, struct Particles *parts, struct Interaction *pot_tot, struct Externos *params, int Nsamp, int factor_desc, int factor_term);
int muestrear_nucleo(char *filename, struct Particles *parts, struct Interaction *pot_tot, struct Externos *params, int Nsamp, int factor_desc, int factor_term);

float ajustar_deltas(struct Particles *parts, struct Interaction *pot_tot, struct Externos *params, float acept_min, float acept_max);
int ajustar_deltas_full(struct Particles *parts, struct Interaction *pot_tot, struct Externos *params, float acept_min, float acept_max);

#endif
