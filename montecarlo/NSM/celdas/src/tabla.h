#ifndef T_H
#define T_H

#include "math.h"

float build_LUT_np(struct Panda_np *panda_np, int N);
float build_LUT_nn(struct Panda_nn *panda_nn, int N);
float build_LUT_pauli(struct Pauli *pauli, int N);
float build_LUT_QCNM(struct QCNM *qcnm, int N);
float build_LUT_Coulomb(struct Coulomb *coul, int N);

float build_LUTF_np(struct Panda_np *panda_np, int N);
float build_LUTF_nn(struct Panda_nn *panda_nn, int N);
float build_LUTF_QCNM(struct QCNM *qcnm, int N);
float build_LUTF_Coulomb(struct Coulomb *coul, int N);
float build_LUTF_Trap(struct Trap *trap, int N);

int build_LUTs(struct Interaction *pot_tot, int Nnp, int Nnn, int Np, int Nq, int Nc, int Nt);
int liberar_LUTs(struct Interaction *pot_tot);
int build_LUTs_F(struct Interaction *pot_tot, int Nnp, int Nnn, int Np, int Nq, int Nc, int Nt);
int liberar_LUTs_F(struct Interaction *pot_tot);

#endif
