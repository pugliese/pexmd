#ifndef G_H
#define G_H

#include "math.h"

int elegir_part(int N);
float uniform();
float boltzmann(float sigma);
float max(float a, float b);
float min(float a, float b);
float norma(float* v);
int min_vec(float* v, int n);
int max_vec(float* v, int n);
int shuffle_array(int *array, int n);

#endif
