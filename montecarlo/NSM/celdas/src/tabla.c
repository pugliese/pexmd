#include "interaccion.h"
#include "general.h"
#include "tabla.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

// ------------------------- TABLAS DE POTENCIAL -------------------------------

float build_LUT_np(struct Panda_np *panda_np, int N){
  float dr2 = panda_np->rcut2/N, r2 = 0, r;
  panda_np->shift = 0;
  panda_np->shift = interaction_panda_np(panda_np->rcut, panda_np);
  panda_np->LUT = (float *) malloc(N*sizeof(float));
  for(int i = 0; i < N; i++){
    r2 = (i+1)*dr2;
    r = sqrt(r2);
    panda_np->LUT[i] = interaction_panda_np(r, panda_np);
  }
  panda_np->dr2 = dr2;
  return dr2;
}

float build_LUT_nn(struct Panda_nn *panda_nn, int N){
  double dr2 = panda_nn->rcut2/N, r2 = 0, r;
  panda_nn->shift = 0;
  panda_nn->shift = interaction_panda_nn(panda_nn->rcut, panda_nn);
  panda_nn->LUT = (float *) malloc(N*sizeof(float));
  for(int i = 0; i < N; i++){
    r2 = (i+1)*dr2;
    r = sqrt(r2);
    panda_nn->LUT[i] = interaction_panda_nn(r, panda_nn);
  }
  panda_nn->dr2 = dr2;
  return dr2;
}

float build_LUT_pauli(struct Pauli *pauli, int N){
  float ds2 = pauli->scut2/N, s2 = 0;
  float p[3] = {0, 0, 0};
  pauli->shift = 0;
  pauli->shift = interaction_pauli(pauli->scut2*pauli->qo*pauli->qo, p, p, pauli);
  pauli->LUT = (float *) malloc(N*sizeof(float));
  for(int i = 0; i < N; i++){
    s2 = (i+1)*ds2;
    pauli->LUT[i] = interaction_pauli(s2*pauli->qo*pauli->qo, p, p, pauli);
  }
  pauli->ds2 = ds2;
  return ds2;
}

float build_LUT_QCNM(struct QCNM *qcnm, int N){
  float dr2 = qcnm->rcut2/N, r2 = 0;
  qcnm->shift = 0;
  qcnm->shift = interaction_QCNM(qcnm->rcut, qcnm);
  qcnm->LUT = (float *) malloc(N*sizeof(float));
  for(int i = 0; i < N; i++){
    r2 = (i+1)*dr2;
    qcnm->LUT[i] = interaction_QCNM(sqrt(r2), qcnm);
  }
  qcnm->dr2 = dr2;
  return dr2;
}

float build_LUT_Coulomb(struct Coulomb *coul, int N){
  float dr2 = coul->rcut2/N, r2 = 0;
  coul->shift = 0;
  coul->shift = interaction_Coulomb(coul->rcut, coul);
  coul->LUT = (float *) malloc(N*sizeof(float));
  for(int i = 0; i < N; i++){
    r2 = (i+1)*dr2;
    coul->LUT[i] = interaction_Coulomb(sqrt(r2), coul);
  }
  coul->dr2 = dr2;
  return dr2;
}

float build_LUT_Trap(struct Trap *trap, int N){
  float dr2 = trap->rcut2/N, r2 = 0;
  trap->LUT = (float *) malloc(N*sizeof(float));
  for(int i = 0; i < N; i++){
    r2 = (i+1)*dr2;
    trap->LUT[i] = interaction_trap(r2, trap);
  }
  trap->dr2 = dr2;
  return dr2;
}

// -------------------------- TABLAS DE FUERZAS --------------------------------

float build_LUTF_np(struct Panda_np *panda_np, int N){
  float dr2 = panda_np->dr2, r2 = 0, r;
  panda_np->LUT_F = (float *) malloc(N*sizeof(float));
  for(int i = 0; i < N; i++){
    r2 = (i+1)*dr2;
    r = sqrt(r2);
    panda_np->LUT_F[i] = force_mod_panda_np(r, panda_np);
  }
  return dr2;
}

float build_LUTF_nn(struct Panda_nn *panda_nn, int N){
  double dr2 = panda_nn->dr2, r2 = 0, r;
  panda_nn->LUT_F = (float *) malloc(N*sizeof(float));
  for(int i = 0; i < N; i++){
    r2 = (i+1)*dr2;
    r = sqrt(r2);
    panda_nn->LUT_F[i] = force_mod_panda_nn(r, panda_nn);
  }
  return dr2;
}

float build_LUTF_QCNM(struct QCNM *qcnm, int N){
  float dr2 = qcnm->dr2, r2 = 0;
  qcnm->LUT_F = (float *) malloc(N*sizeof(float));
  for(int i = 0; i < N; i++){
    r2 = (i+1)*dr2;
    qcnm->LUT_F[i] = force_mod_QCNM(sqrt(r2), qcnm);
  }
  return dr2;
}

float build_LUTF_Coulomb(struct Coulomb *coul, int N){
  float dr2 = coul->dr2, r2 = 0;
  coul->LUT_F = (float *) malloc(N*sizeof(float));
  for(int i = 0; i < N; i++){
    r2 = (i+1)*dr2;
    coul->LUT_F[i] = force_mod_Coulomb(sqrt(r2), coul);
  }
  return dr2;
}

float build_LUTF_Trap(struct Trap *trap, int N){
  float dr2 = trap->dr2, r2 = 0;
  trap->LUT_F = (float *) malloc(N*sizeof(float));
  for(int i = 0; i < N; i++){
    r2 = (i+1)*dr2;
    trap->LUT[i] = force_mod_trap(r2, trap);
  }
  return dr2;
}

// --------------------- CONSTRUCTORES/DESTRUCTORES ----------------------------

int build_LUTs(struct Interaction *pot_tot, int Nnp, int Nnn, int Np, int Nq, int Nc, int Nt){
  build_LUT_np(pot_tot->panda_np, Nnp);
  build_LUT_nn(pot_tot->panda_nn, Nnn);
  build_LUT_pauli(pot_tot->pauli, Np);
  build_LUT_QCNM(pot_tot->qcnm, Nq);
  build_LUT_Coulomb(pot_tot->coul, Nc);
  build_LUT_Trap(pot_tot->trap, Nt);
  return 0;
}

int liberar_LUTs(struct Interaction *pot_tot){
  free(pot_tot->pauli->LUT);
  free(pot_tot->panda_nn->LUT);
  free(pot_tot->panda_np->LUT);
  free(pot_tot->qcnm->LUT);
  free(pot_tot->coul->LUT);
  free(pot_tot->trap->LUT);
  return 0;
}

int build_LUTs_F(struct Interaction *pot_tot, int Nnp, int Nnn, int Np, int Nq, int Nc, int Nt){
  build_LUT_np(pot_tot->panda_np, Nnp);
  build_LUT_nn(pot_tot->panda_nn, Nnn);
  build_LUT_pauli(pot_tot->pauli, Np);
  build_LUT_QCNM(pot_tot->qcnm, Nq);
  build_LUT_Coulomb(pot_tot->coul, Nc);
  build_LUT_Trap(pot_tot->trap, Nt);
  return 0;
}

int liberar_LUTs_F(struct Interaction *pot_tot){
  free(pot_tot->panda_nn->LUT_F);
  free(pot_tot->panda_np->LUT_F);
  free(pot_tot->qcnm->LUT_F);
  free(pot_tot->coul->LUT_F);
  free(pot_tot->trap->LUT_F);
  return 0;
}
