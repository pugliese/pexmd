#include "general.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

int elegir_part(int N){
  int i = rand();
  while (i > (RAND_MAX/N)*N){
    i = rand();
  }
  i = i % N;
  return i;
}

float uniform(){
  float p = ((float) rand())/((float) RAND_MAX);
  return p;
}

double normaldist()
{
  double x,y,r2,c;
  r2 = 0.0;
  while(r2 == 0.0 || r2 > 1.0){
    x = 2.0*((double) rand() / (double)RAND_MAX) - 1.0;
    y = 2.0*((double) rand() / (double)RAND_MAX) - 1.0;
    r2 = x*x + y*y;
  }
  c = sqrt(-2.0*log(r2)/r2);
  return x*c;
}

float boltzmann(float sigma){
  return normaldist()*sigma;
}

/*
float boltzmann(float sigma){
  float sum = 0;
  int n_samp = 24;
  for(int k = 0; k < n_samp; k++){
    sum = sum + uniform();
  }
  sum = (sum - 0.5*n_samp)*sqrt(12.0/n_samp)*sigma;
  return sum;
}
*/

float max(float a, float b){
  if (a>b){
    return a;
  }else{
    return b;
  }
}

float min(float a, float b){
  if (a<b){
    return a;
  }else{
    return b;
  }
}

float norma(float* v){
  float rsq = 0;
  for(int k = 0; k < 3; k++){
    rsq += v[k]*v[k];
  }
  return rsq;
}

int min_vec(float* v, int n){
  float rsq_min = norma(v), rsq;
  int idx_min = 0;
  for(int i = 1; i < n; i++){
    rsq = norma(v+3*i);
    if (rsq < rsq_min){
      idx_min = i;
      rsq_min = rsq;
    }
  }
  return idx_min;
}

int max_vec(float* v, int n){
  float rsq_max = norma(v), rsq;
  int idx_max = 0;
  for(int i = 1; i < n; i++){
    rsq = norma(v+3*i);
    if (rsq > rsq_max){
      idx_max = i;
      rsq_max = rsq;
    }
  }
  return idx_max;
}

int shuffle_array(int *array, int n){
  if (n > 1){
    int i, t;
    for (i = 0; i < n - 1; i++){
      size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
      t = array[j];
      array[j] = array[i];
      array[i] = t;
    }
  }
  return 0;
}
