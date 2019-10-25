#ifndef G_H
#define G_H

#include "math.h"

int elegir_part(int N);
float uniform();
float boltzmann(float sigma);
float max(float a, float b);
float min(float a, float b);
float producto_interno(float *u, float* v);
int producto_vectorial(float *u, float* v, float* u_x_v);
float norma(float* v);
float distancia(float *q1, float *q2, int *delta_idx, float l);
float distancia_p(float *p1, float *p2);
int min_vec(float* v, int n);
int max_vec(float* v, int n);
int shuffle_array(int *array, int n);

#endif
