#ifndef I_H
#define I_H

#include "math.h"

struct Interaction{
  struct Pauli* pauli;
  struct Panda_nn* panda_nn;
  struct Panda_np* panda_np;
  struct QCNM* qcnm;
  struct Coulomb* coul;
  struct Trap* trap;
  float rcut;
  double energy_pauli;
  double energy_panda_nn;
  double energy_panda_np;
  double energy_qcnm;
  double energy_coul;
};

struct Pauli {
  float qo;
  float po;
  float D;
  float scut2;
  float shift;
  float *LUT;
  float *LUT_F;
  float ds2;
};

struct Panda_np {
  float mu_r;
  float mu_a;
  float V_r;
  float V_a;
  float rcut;
  float rcut2;
  float shift;
  float *LUT;
  float *LUT_F;
  float dr2;
};

struct Panda_nn {
  float mu_o;
  float V_o;
  float rcut;
  float rcut2;
  float shift;
  float *LUT;
  float *LUT_F;
  float dr2;
};

struct QCNM {
  float V_o;
  float p_1;
  float p_2;
  float r_1;
  float r_2;
  float d;
  float a;
  float rcut;
  float rcut2;
  float shift;
  float *LUT;
  float *LUT_F;
  float dr2;
};

struct Coulomb {
  float q2;
  float lambda;
  float rcut;
  float rcut2;
  float shift;
  float *LUT;
  float *LUT_F;
  float dr2;
};

struct Trap{
  float power;
  float k;
  float rcut2;
  float *LUT;
  float *LUT_F;
  float dr2;
};

int reset_energies(struct Interaction* pot_tot);

float interaction_QCNM(float r, struct QCNM *qcnm);
float interaction_trap(float rsq, struct Trap *trap);
float interaction_panda_nn(float r, struct Panda_nn *panda_nn);
float interaction_panda_np(float r, struct Panda_np *panda_np);
float interaction_pauli(float rsq, float *p1, float *p2, struct Pauli *pauli);
float interaction_pauli_factorizado(float rsq, float *p1, float *p2, struct Pauli *pauli);
float interaction_Coulomb(float r, struct Coulomb *coul);

float force_mod_QCNM(float r, struct QCNM *qcnm);
float force_mod_trap(float rsq, struct Trap *trap);
float force_mod_panda_nn(float r, struct Panda_nn *panda_nn);
float force_mod_panda_np(float r, struct Panda_np *panda_np);
float fgorce_mod_pauli(float rsq, float *p1, float *p2, struct Pauli *pauli, float *gorce_mod);
float force_mod_Coulomb(float r, struct Coulomb *coul);

int interaction(int t1, int t2, float rsq, float *p1, float *p2, struct Interaction *pot_tot);
float fgorce_mod(int t1, int t2, float rsq, float *p1, float *p2, struct Interaction *pot_tot, float *gorce_mod);
float eval_LUT(float x, float *LUT, float xcut, float dx);
int interaction_sin_LUT(int t1, int t2, float rsq, float *p1, float *p2, struct Interaction *pot_tot);
float fgorce_mod_sin_LUT(int t1, int t2, float rsq, float *p1, float *p2, struct Interaction *pot_tot, float *gorce_mod);

#endif
