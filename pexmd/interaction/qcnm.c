#include "qcnm.h"
#include <stdio.h>
#include <time.h>
#include <unistd.h>

float pair_energ(float r, float energy_cut, float Vo, float r1, float p1,
                  float r2, float p2, float d, float a){

  double qexp = exp((r-d)/a);
  double r_red1_pot = pow(r1/r, p1);
  double r_red2_pot = pow(r2/r, p2);
  return Vo*(r_red1_pot - r_red2_pot)/(1+qexp) - energy_cut;

}

float pair_force(float r, float Vo, float r1, float p1,
                  float r2, float p2, float d, float a){

  double qexp = exp((r-d)/a);
  double r_red1_pot = pow(r1/r, p1);
  double r_red2_pot = pow(r2/r, p2);
  return ( (p1*r_red1_pot - p2*r_red2_pot)/r + (r_red1_pot - r_red2_pot)*qexp/(a*(1+qexp)) )*Vo/((1+qexp)*r);

}

// Desacopladas, devuelve modulo, encapsula potencial y calcula energy_cut 1 vez
float forces(float *x, long int* pairs, long int npairs, float Vo, float r1,
      float p1, float r2, float p2, float d, float a, float rcut, float *force){

  float energy = 0;
  float energy_cut = 0;
  energy_cut = pair_energ(rcut, 0, Vo, r1, p1, r2, p2, d, a);

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
      float m_force = pair_force(r, Vo, r1, p1, r2, p2, d, a);
      for (int k = 0; k < 3; k++){
        force[3*i+k] += m_force*delta_r[k];
        force[3*j+k] -= m_force*delta_r[k];
      }
      energy += pair_energ(r, energy_cut, Vo, r1, p1, r2, p2, d, a);
    }
  }

  return energy;

}
