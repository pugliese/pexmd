#include "pressure.h"
#include <stdio.h>
#include <time.h>
#include <unistd.h>


float pair_force_qcnm(float r, float Vo, float r1, float p1,
                  float r2, float p2, float d, float a){

  double qexp = exp((r-d)/a);
  double r_red1_pot = pow(r1/r, p1);
  double r_red2_pot = pow(r2/r, p2);
  return Vo*( (p1*r_red1_pot - p2*r_red2_pot)/r + (r_red1_pot - r_red2_pot)*qexp/(a*(1+qexp)) )/((1+qexp)*r);

}

float pressure_qcnm_layers(float *x, long int* pairs, long int npairs, float Vo, float r1,
      float p1, float r2, float p2, float d, float a, float rcut, float *force, float L, int ls){

  float pres = 0.0;
  for(int e1 = -ls; e1 <= ls; e1++){
    for(int e2 = -ls; e2 <= ls; e2++){
      for(int e3 = -ls; e3 <= ls; e3++){

        float celda[3];
        celda[0] = L*e1;
        celda[1] = L*e2;
        celda[2] = L*e3;

        for (int l = 0; l < npairs; l++) {
          float delta_r[3];
          int i = pairs[2*l];
          int j = pairs[2*l+1];

          for (int k = 0; k < 3; k++) {
            delta_r[k] = x[3*i+k] + celda[0] - x[3*j+k];
          }

          float r = 0;
          for (int k = 0; k < 3; k++) {
            r = r + delta_r[k]*delta_r[k];
          }
          r = sqrt(r);

          if (r<rcut) {
            float m_force = pair_force_qcnm(r, Vo, r1, p1, r2, p2, d, a);
            for (int k = 0; k < 3; k++){
              pres += delta_r[k] * m_force * delta_r[k];
              force[3*i+k] += m_force*delta_r[k];
              force[3*j+k] -= m_force*delta_r[k];
            }
          }
        }
      }
    }
  }
  return pres;
}

float pressure_lj_PBC(float *x, long int* pairs, long int npairs, float eps,
             float sigma, float rcut, float *force, float L) {
  float pres = 0.0;
  float ljf1 = 48;
  float ljf2 = 24;
  float lje1 = 4;
  float lje2 = 4;
  float rcutsq = rcut * rcut;
  for (int ii = 0; ii < npairs; ii++) {
    float delr[3];
    long int i = pairs[2*ii];
    long int j = pairs[2*ii + 1];

    for (int k = 0; k < 3; k++) {
      delr[k] = x[3*i + k] - x[3*j + k];
      delr[k] = delr[k] + L*( (delr[k] < -0.5*L) - (0.5*L < delr[k]) );
    }

    float rsq = 0.0;
    for (int k = 0; k < 3; k++) {
      rsq += delr[k] * delr[k];
    }

    if (rsq < rcutsq) {
      float r2inv = 1.0/rsq;
      float r6inv = r2inv * r2inv * r2inv;
      float forcelj = r2inv * r6inv * (ljf1 * r6inv - ljf2);
      for (int k = 0; k < 3; k++) {
        pres += delr[k] * forcelj * delr[k];
        force[3*i + k] += forcelj * delr[k];
        force[3*j + k] -= forcelj * delr[k];
      }
    }
  }
  return pres;
}



float pair_force(float s2, float D, float qo2){

  float pexp = exp(-s2/2);
  return D*pexp/qo2;

}

float pair_gorce(float s2, float D, float po2){

  float pexp = exp(-s2/2);
  return D*pexp/po2;

}

float pressure_pauli_PBC(float *x, float *p, long int* pairs, long int npairs,
  float D, float qo, float po, float scut, float *force, float *gorce, float L){

  float qF = 0;
  double qF_aux = 0;
  float qo2 = qo*qo;
  float po2 = po*po;
  float scut2 = scut*scut;
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
        qF += force_mod*delta_q[k]*delta_q[k];
      }
    }
  }
  qF_aux = qF;
  return qF_aux;

}

float delta_fases(float *x, float *p, long int* pairs, long int npairs, float *dq, float *dp, float L){

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

    dq[l] = 0;
    dp[l] = 0;
    for (int k = 0; k < 3; k++) {
      dq[l] += delta_q[k]*delta_q[k];
      dp[l] += delta_p[k]*delta_p[k];
    }
    dq[l] = sqrt(dq[l]);
    dp[l] = sqrt(dp[l]);

  }
  return 0;
}

float delta_fases2(float *x, float *p, long int n, float *dq, float *dp, float L){
  int l = 0;
  for (int i = 1; i < n; i++) {
    for(int j = 0; j < i; j++) {
      float delta_q[3];
      float delta_p[3];
      for (int k = 0; k < 3; k++) {
        delta_q[k] = x[3*i+k] - x[3*j+k];
        delta_p[k] = p[3*i+k] - p[3*j+k];
        delta_q[k] = delta_q[k] + L*(delta_q[k] < -0.5*L) - L*(0.5*L < delta_q[k]);
      }
      dq[l] = 0;
      dp[l] = 0;
      for (int k = 0; k < 3; k++) {
        dq[l] += delta_q[k]*delta_q[k];
        dp[l] += delta_p[k]*delta_p[k];
      }
      dq[l] = sqrt(dq[l]);
      dp[l] = sqrt(dp[l]);
      l++;
    }
  }
  return 0;
}


float delta_fases_sin_PBC(float *x, float *p, long int* pairs, long int npairs, float *dq, float *dp){

  for (int l = 0; l < npairs; l++) {
    float delta_q[3];
    float delta_p[3];
    int i = pairs[2*l];
    int j = pairs[2*l+1];

    for (int k = 0; k < 3; k++) {
      delta_q[k] = x[3*i+k] - x[3*j+k];
      delta_p[k] = p[3*i+k] - p[3*j+k];
    }

    dq[l] = 0;
    dp[l] = 0;
    for (int k = 0; k < 3; k++) {
      dq[l] += delta_q[k]*delta_q[k];
      dp[l] += delta_p[k]*delta_p[k];
    }
    dq[l] = sqrt(dq[l]);
    dp[l] = sqrt(dp[l]);

  }
  return 0;
}

float Gr(float* x, int N, float dr, float L, float* g, int ng) {
  float rij;
  int k;
  float delta_q[3];

  for(int i = 1; i < N; i++){
    for(int j = 0; j < i; j++){
      rij = 0;
      for (int k = 0; k < 3; k++) {
        delta_q[k] = x[3*i+k] - x[3*j+k];
        delta_q[k] = delta_q[k] + L*(delta_q[k] < -0.5*L) - L*(0.5*L < delta_q[k]);
        rij += delta_q[k]*delta_q[k];
      }
      rij = sqrt(rij);
      k = floor(rij/dr);
      g[k]++;
    }
  }
  float rho = N/(L*L*L);
  for (int i = 1; i < ng; i++){
    g[i] = g[i]/(2*M_PI*dr*dr*dr*i*i*rho*N);
  }
  return 0;
}
