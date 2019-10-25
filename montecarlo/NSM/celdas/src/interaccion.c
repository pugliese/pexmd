#include "interaccion.h"
#include "general.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

int reset_energies(struct Interaction* pot_tot){
  pot_tot->energy_pauli = 0;
  pot_tot->energy_panda_nn = 0;
  pot_tot->energy_panda_np = 0;
  pot_tot->energy_qcnm = 0;
  pot_tot->energy_coul = 0;
  return 0;
}

// ----------------------------- POTENCIAL -------------------------------------

inline float interaction_QCNM(float r, struct QCNM *qcnm){
  float pot = 0;
  if (r <= qcnm->rcut){
    pot = qcnm->V_o*(pow(qcnm->r_1/r, qcnm->p_1) - pow(qcnm->r_2/r, qcnm->p_2))/(1 + exp((r-qcnm->d)/qcnm->a)) - qcnm->shift;
  }
  return pot;
}

inline float interaction_trap(float rsq, struct Trap *trap){
  float pot = trap->k*pow(rsq, trap->power/2);
  return pot;
}

inline float interaction_Coulomb(float r, struct Coulomb *coul){
  float pot = 0;
  if (r <= coul->rcut){
    pot = coul->q2*exp(-r/coul->lambda)/r - coul->shift;
  }
  return pot;
}

inline float interaction_panda_nn(float r, struct Panda_nn *panda_nn){
  float pot = 0;
  if (r <= panda_nn->rcut){
    pot = panda_nn->V_o*exp(-panda_nn->mu_o*r)/r - panda_nn->shift;
  }
  return pot;
}

inline float interaction_panda_np(float r, struct Panda_np *panda_np){
  float pot = 0;
  if (r <= panda_np->rcut){
    pot = (panda_np->V_r*exp(-panda_np->mu_r*r) - panda_np->V_a*exp(-panda_np->mu_a*r))/r - panda_np->shift;
  }
  return pot;
}

inline float interaction_pauli(float rsq, float *p1, float *p2, struct Pauli *pauli){
  float pot = 0;
  float s2 = rsq/(pauli->qo*pauli->qo) + distancia_p(p1, p2)/(pauli->po*pauli->po);
  if (s2 <= pauli->scut2){
    pot = pauli->D*exp(-0.5*s2) - pauli->shift;
  }
  return pot;
}

inline float interaction_pauli_factorizado(float rsq, float *p1, float *p2, struct Pauli *pauli){
  float pot = 0;
  float r2_norm = rsq/(pauli->qo*pauli->qo);
  if (r2_norm <= pauli->scut2){
    float p2_norm = distancia_p(p1, p2)/(pauli->po*pauli->po);
    pot = (pauli->D*exp(-0.5*r2_norm) - pauli->shift)*exp(-0.5*p2_norm);
  }
  return pot;
}

// ----------------------------- F/GUERZAS -------------------------------------

inline float force_mod_panda_nn(float r, struct Panda_nn *panda_nn){
  float f_mod = 0;
  if (r <= panda_nn->rcut){
    f_mod = panda_nn->V_o*exp(-panda_nn->mu_o*r)*(1+panda_nn->mu_o*r)/(r*r*r);
  }
  return f_mod;
}

inline float force_mod_panda_np(float r, struct Panda_np *panda_np){
  float f_mod = 0;
  if (r <= panda_np->rcut){
    f_mod = (panda_np->V_r*exp(-panda_np->mu_r*r)*(1+panda_np->mu_r*r) - panda_np->V_a*exp(-panda_np->mu_a*r)*(1+panda_np->mu_a*r))/(r*r*r);
  }
  return f_mod;
}

inline float fgorce_mod_pauli(float rsq, float *p1, float *p2, struct Pauli *pauli, float *gorce_mod){
  float f_mod;
  float r2_red = rsq/(pauli->qo*pauli->qo);
  if (r2_red < pauli->scut2){
    float p2_red = distancia_p(p1, p2)/(pauli->po*pauli->po);
    float exp_q = exp(-0.5*r2_red);
    float exp_p = exp(-0.5*p2_red);
    f_mod = pauli->D*exp_p*exp_q/(pauli->qo*pauli->qo);
    *gorce_mod = (pauli->D*exp_q - pauli->shift)*exp_p/(pauli->po*pauli->po);
  }else{
    f_mod = 0;
    *gorce_mod = 0;
  }
  return f_mod;
}

inline float force_mod_QCNM(float r, struct QCNM *qcnm){
  float f_mod = 0;
  if (r <= qcnm->rcut){
    float t1 = pow(qcnm->r_1/r, qcnm->p_1);
    float t2 = pow(qcnm->r_2/r, qcnm->p_2);
    float qexp = exp((r - qcnm->d)/qcnm->a);
    f_mod = (qcnm->p_1*t1 - qcnm->p_2*t1)/r + (t1 - t2)/(qcnm->a*(1+1/qexp));
    f_mod = f_mod*qcnm->V_o/(1 + qexp);
  }
  return f_mod;
}

inline float force_mod_Coulomb(float r, struct Coulomb *coul){
  float f_mod = 0;
  if (r <= coul->rcut){
    f_mod = coul->q2*exp(-r/coul->lambda)*(1 + r/coul->lambda)/(r*r*r);
  }
  return f_mod;
}

inline float force_mod_trap(float rsq, struct Trap *trap){
  float f_mod = trap->k*pow(rsq, trap->power/2-1);
  return f_mod;
}

// ---------------------------- INTERACCIONES ----------------------------------

inline int interaction_sin_LUT(int t1, int t2, float rsq, float *p1, float *p2, struct Interaction *pot_tot){
  float r = sqrt(rsq);
  int res = 10*t1 + t2;
  pot_tot->energy_qcnm += interaction_QCNM(r, pot_tot->qcnm);
  if (t1/2 == t2/2){ // Same isospin
    pot_tot->energy_panda_nn += interaction_panda_nn(r, pot_tot->panda_nn);
    if(t1/2 == 0){  // Both protons
      pot_tot->energy_coul += interaction_Coulomb(r, pot_tot->coul);
    }
    if (t1 == t2){ // Also same spin
      pot_tot->energy_pauli += interaction_pauli_factorizado(rsq, p1, p2, pot_tot->pauli);
    }
  }else{  // Different isospin
    pot_tot->energy_panda_np += interaction_panda_np(r, pot_tot->panda_np);
  }
  return res;
}

float fgorce_mod_sin_LUT(int t1, int t2, float rsq, float *p1, float *p2, struct Interaction *pot_tot, float *gorce_mod){
  float r = sqrt(rsq);
  float force_mod = force_mod_QCNM(r, pot_tot->qcnm);
  *gorce_mod = 0;
  if(pot_tot->trap->k != 0) force_mod += force_mod_trap(rsq, pot_tot->trap);
  if (t1/2 == t2/2){
    force_mod += force_mod_panda_nn(r, pot_tot->panda_nn);
    if (t1/2 == 0){
      force_mod += force_mod_Coulomb(r, pot_tot->coul);
    }
    if (t1 == t2){
      force_mod += fgorce_mod_pauli(rsq, p1, p2, pot_tot->pauli, gorce_mod);
    }
  }else{
    force_mod = force_mod_panda_np(r, pot_tot->panda_np);
  }
  return force_mod;
}

// ---------------------------- TABLAS -----------------------------------------

inline int interaction(int t1, int t2, float rsq, float *p1, float *p2, struct Interaction *pot_tot){
  int res = 10*t1 + t2;
  pot_tot->energy_qcnm += eval_LUT(rsq, pot_tot->qcnm->LUT, pot_tot->qcnm->rcut2, pot_tot->qcnm->dr2);
  if (t1/2 == t2/2){ // Same isospin
    pot_tot->energy_panda_nn += eval_LUT(rsq, pot_tot->panda_nn->LUT, pot_tot->panda_nn->rcut2, pot_tot->panda_nn->dr2);
    if(t1/2 == 0){  // Both protons
      pot_tot->energy_coul += eval_LUT(rsq, pot_tot->coul->LUT, pot_tot->coul->rcut2, pot_tot->coul->dr2);
    }
    if (t1 == t2){ // Also same spin
      pot_tot->energy_pauli += interaction_pauli_factorizado(rsq, p1, p2, pot_tot->pauli);
    }
  }else{  // Different isospin
    pot_tot->energy_panda_np += eval_LUT(rsq, pot_tot->panda_np->LUT, pot_tot->panda_np->rcut2, pot_tot->panda_np->dr2);
  }
  return res;
}

inline float fgorce_mod(int t1, int t2, float rsq, float *p1, float *p2, struct Interaction *pot_tot, float *gorce_mod){
  float force_mod = eval_LUT(rsq, pot_tot->qcnm->LUT_F, pot_tot->qcnm->rcut2, pot_tot->qcnm->dr2);
  if(pot_tot->trap->k != 0) force_mod += eval_LUT(rsq, pot_tot->trap->LUT_F, pot_tot->trap->rcut2, pot_tot->trap->dr2);
  if (t1/2 == t2/2){
    force_mod += eval_LUT(rsq, pot_tot->panda_nn->LUT_F, pot_tot->panda_nn->rcut2, pot_tot->panda_nn->dr2);
    if (t1/2 == 0){
      force_mod += eval_LUT(rsq, pot_tot->coul->LUT_F, pot_tot->coul->rcut2, pot_tot->coul->dr2);
    }
    if (t1 == t2){
      force_mod += fgorce_mod_pauli(rsq, p1, p2, pot_tot->pauli, gorce_mod);
    }else{
      *gorce_mod = 0;
    }
  }else{
    force_mod = eval_LUT(rsq, pot_tot->panda_np->LUT_F, pot_tot->panda_np->rcut2, pot_tot->panda_np->dr2);
  }
  return force_mod;
}

float eval_LUT(float x, float *LUT, float xcut, float dx){
  if(x < dx){
    return LUT[0];
  }
  float res = 0;
  if (x < xcut){
    float a = x/dx;
    int i = (int) a;
    res = (LUT[i] - LUT[i-1])*(a-i) + LUT[i-1];
  }
  return res;
}
