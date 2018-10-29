#include "qcnm_nn.h"
#include <stdio.h>
#include <time.h>
#include <unistd.h>

float pair_energ(float r, float energy_cut, float mur,
              float Dr, float mua, float Da){

  double exp_r = exp(-mur*r)/r;
  double exp_a = exp(-mua*r)/r;
  return Dr*exp_r - Da*exp_a - energy_cut;

}

float pair_force(float r, float mur, float Dr,
              float mua, float Da){

  double exp_r = exp(-mur*r)/r;
  double exp_a = exp(-mua*r)/r;
  return ( Dr*exp_r*(mur+1.0/r) - Da*exp_a*(mua+1.0/r) )/r;

}

// Desacopladas, devuelve modulo, encapsula potencial y calcula energy_cut 1 vez
float forces(float *x, long int* pairs, long int npairs, float Dr,
          float mur, float Da, float mua, float rcut, float *force){

  float energy = 0;
  float energy_cut = 0;
  energy_cut = pair_energ(rcut, 0, mur, Dr, mua, Da);

  for (int l = 0; l < npairs; l++) {
    float delta_r[3];
    int i = pairs[2*l];
    int j = pairs[2*l+1];

    for (int k = 0; k < 3; k++) {
      delta_r[k] = x[3*i+k]-x[3*j+k];
    }

    float r = 0;
    for (int k = 0; k < 3; k++) {
      r = r + delta_r[k]*delta_r[k];
    }
    r = sqrt(r);

    if (r<rcut) {
      float m_force = pair_force(r, mur, Dr, mua, Da);
      for (int k = 0; k < 3; k++){
        force[3*i+k] += m_force*delta_r[k];
        force[3*j+k] -= m_force*delta_r[k];
      }
      energy += pair_energ(r, energy_cut, mur, Dr, mua, Da);
    }
  }

  return energy;

}
