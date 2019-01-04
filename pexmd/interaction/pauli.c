#include "pauli.h"
#include <stdio.h>
#include <time.h>
#include <unistd.h>

float pair_energ(float s2, float energy_cut, float D){

  float pexp = exp(-s2/2);
  return D*pexp - energy_cut;

}

float pair_force(float s2, float D, float qo2){

  float pexp = exp(-s2/2);
  return D*pexp/qo2;

}

float pair_gorce(float s2, float D, float po2){

  float pexp = exp(-s2/2);
  return D*pexp/po2;

}

float fgorces(float *x, float *p, long int* pairs, long int npairs,
  float D, float qo, float po, float scut, float *force, float *gorce){

  float energy = 0;
  float qo2 = qo*qo;
  float po2 = po*po;
  float scut2 = scut*scut;
  float energy_cut = pair_energ(scut2, 0, D);

  for (int l = 0; l < npairs; l++) {
    float delta_q[3];
    float delta_p[3];
    int i = pairs[2*l];
    int j = pairs[2*l+1];

    for (int k = 0; k < 3; k++) {
      delta_q[k] = x[3*i+k] - x[3*j+k];
      delta_p[k] = p[3*i+k] - p[3*j+k];
    }

    float q2 = 0;
    for (int k = 0; k < 3; k++) {
      q2 = q2 + delta_q[k]*delta_q[k];
    }

    float p2 = 0;
    for (int k = 0; k < 3; k++) {
      p2 = p2 + delta_p[k]*delta_p[k];
    }

    float s2 = q2/qo2 + p2/po2;

    if (s2<scut2) {
      float force_mod = pair_force(s2, D, qo2);
      float gorce_mod = force_mod*qo2/po2;
      for (int k = 0; k < 3; k++){
        force[3*i+k] += force_mod*delta_q[k];
        force[3*j+k] -= force_mod*delta_q[k];
        gorce[3*i+k] += gorce_mod*delta_p[k];
        gorce[3*j+k] -= gorce_mod*delta_p[k];
      }
      energy += pair_energ(s2, energy_cut, D);
    }
  }

  return energy;

}

float fgorces_PBC(float *x, float *p, long int* pairs, long int npairs,
  float D, float qo, float po, float scut, float *force, float *gorce, float L){

  float energy = 0;
  float qo2 = qo*qo;
  float po2 = po*po;
  float scut2 = scut*scut;
  float energy_cut = pair_energ(scut2, 0, D);
  for (int l = 0; l < npairs; l++) {
    float delta_q[3];
    float delta_p[3];
    int i = pairs[2*l];
    int j = pairs[2*l+1];

    for (int k = 0; k < 3; k++) {
      delta_q[k] = x[3*i+k] - x[3*j+k];
      delta_p[k] = p[3*i+k] - p[3*j+k];
      delta_q[k] = delta_q[k] + L*(delta_q[k] < -0.5*L) - L*(0.5*L < delta_q[k]);
    }

    float q2 = 0;
    for (int k = 0; k < 3; k++) {
      q2 = q2 + delta_q[k]*delta_q[k];
    }

    float p2 = 0;
    for (int k = 0; k < 3; k++) {
      p2 = p2 + delta_p[k]*delta_p[k];
    }

    float s2 = q2/qo2 + p2/po2;

    if (s2<scut2) {
      float force_mod = pair_force(s2, D, qo2);
      float gorce_mod = force_mod*qo2/po2;
      for (int k = 0; k < 3; k++){
        force[3*i+k] += force_mod*delta_q[k];
        force[3*j+k] -= force_mod*delta_q[k];
        gorce[3*i+k] += gorce_mod*delta_p[k];
        gorce[3*j+k] -= gorce_mod*delta_p[k];
      }
      energy += pair_energ(s2, energy_cut, D);
    }
  }

  return energy;

}

float forces(float *x, float *p, long int* pairs, long int npairs,
  float D, float qo, float po, float scut, float *force){

  float energy = 0;
  float qo2 = qo*qo;
  float po2 = po*po;
  float scut2 = scut*scut;
  float energy_cut = pair_energ(scut2, 0, D);

  for (int l = 0; l < npairs; l++) {
    float delta_q[3];
    float delta_p[3];
    int i = pairs[2*l];
    int j = pairs[2*l+1];

    for (int k = 0; k < 3; k++) {
      delta_q[k] = x[3*i+k] - x[3*j+k];
      delta_p[k] = p[3*i+k] - p[3*j+k];
    }

    float q2 = 0;
    for (int k = 0; k < 3; k++) {
      q2 = q2 + delta_q[k]*delta_q[k];
    }

    float p2 = 0;
    for (int k = 0; k < 3; k++) {
      p2 = p2 + delta_p[k]*delta_p[k];
    }

    float s2 = q2/qo2 + p2/po2;

    if (s2<scut2) {
      float force_mod = pair_force(s2, D, qo2);
      for (int k = 0; k < 3; k++){
        force[3*i+k] += force_mod*delta_q[k];
        force[3*j+k] -= force_mod*delta_q[k];
      }
      energy += pair_energ(s2, energy_cut, D);
    }
  }

  return energy;

}


float gorces(float *x, float *p, long int* pairs, long int npairs,
  float D, float qo, float po, float scut, float *gorce){

  float energy = 0;
  // Siempre uso los cuadrados, cambio parametros por sus cuadrados?
  float qo2 = qo*qo;
  float po2 = po*po;
  float scut2 = scut*scut;
  float energy_cut = pair_energ(scut2, 0, D);

  for (int l = 0; l < npairs; l++) {
    float delta_q[3];
    float delta_p[3];
    int i = pairs[2*l];
    int j = pairs[2*l+1];

    for (int k = 0; k < 3; k++) {
      delta_q[k] = x[3*i+k] - x[3*j+k];
      delta_p[k] = p[3*i+k] - p[3*j+k];
    }

    float q2 = 0;
    for (int k = 0; k < 3; k++) {
      q2 = q2 + delta_q[k]*delta_q[k];
    }

    float p2 = 0;
    for (int k = 0; k < 3; k++) {
      p2 = p2 + delta_p[k]*delta_p[k];
    }

    float s2 = q2/qo2 + p2/po2;

    if (s2<scut2) {
      float gorce_mod = pair_gorce(s2, D, po2);
      for (int k = 0; k < 3; k++){
        gorce[3*i+k] += gorce_mod*delta_p[k];
        gorce[3*j+k] -= gorce_mod*delta_p[k];
      }
      energy += pair_energ(s2, energy_cut, D);
    }
  }

  return energy;

}
