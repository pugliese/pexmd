#include "general.h"
#include "celda.h"
#include "interaccion.h"
#include "energia.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

int construir_lista_celdas(int K_c, int *idxs){
  int idx = 0;
  for(int k = 0; k <= K_c; k++){
    idxs[3*idx+0] = k;
    idxs[3*idx+1] = 0;
    idxs[3*idx+2] = 0;
    idx++;
  }
  for(int kx = -K_c; kx <= K_c; kx++){
    for(int ky = 1; ky <= K_c; ky++){
      idxs[3*idx+0] = kx;
      idxs[3*idx+1] = ky;
      idxs[3*idx+2] = 0;
      idx++;
    }
  }
  for(int kx = -K_c; kx <= K_c; kx++){
    for(int ky = -K_c; ky <= K_c; ky++){
      for(int kz = 1; kz <= K_c; kz++){
        idxs[3*idx+0] = kx;
        idxs[3*idx+1] = ky;
        idxs[3*idx+2] = kz;
        idx++;
      }
    }
  }
  return 0;
}

//------------------------- ENERGIA SOLAMENTE ----------------------------------

int energia(struct Particles *parts, struct Interaction *pot_tot){
  double kin = 0, rsq;
  reset_energies(pot_tot);
  int K_c = (int) ceil(pot_tot->coul->rcut / parts->l);
  int cantidad = ((2*K_c+1)*(2*K_c+1)*(2*K_c+1) + 1)/2;
  int i, j, mx, my, mz, idx, idxs[3*cantidad];
  int M = parts->M, M3 = M*M*M;
  construir_lista_celdas(K_c, idxs);
  for (int m = 0; m < M3; m++){
    mx = m/(M*M);
    my = (m/M) % M;
    mz = m % M;
    for (int k = 0; k < cantidad; k++){
      idx = (((mx+idxs[3*k]+M) % M)*M + (my+idxs[3*k+1]+M) % M)*M + (mz+idxs[3*k+2]+M) % M;
      i = parts->primero[m];
      if (idxs[3*k]<=1 && idxs[3*k+1]<=1 && idxs[3*k+2]<=1 && idxs[3*k]>=-1 && idxs[3*k+1]>=-1 && idxs[3*k+2]>=-1){
        while (i != -1){
          j = parts->primero[idx];
          while (j !=-1 && j != i){
            rsq = distancia(parts->q+3*i, parts->q+3*j, idxs+3*k, parts->l);
            interaction(parts->type[i], parts->type[j], rsq, parts->p+3*i, parts->p+3*j, pot_tot);
            j = parts->siguiente[j];
          }
          i = parts->siguiente[i];
        }
      }else{
        while (i != -1){
          if(parts->type[i]/2 == 0){
            j = parts->primero[idx];
            while (j !=-1 && j != i){
              if (parts->type[j]/2 == 0){
                rsq = distancia(parts->q+3*i, parts->q+3*j, idxs+3*k, parts->l);
                interaction(parts->type[i], parts->type[j], rsq, parts->p+3*i, parts->p+3*j, pot_tot);
              }
              j = parts->siguiente[j];
            }
          }
          i = parts->siguiente[i];
        }
      }
    }
  }
  // Energia cinetica
  for(int i = 0; i < parts->n; i++){
    for(int k = 0; k < 3; k++){
      kin = kin + (parts->p[3*i+k])*(parts->p[3*i+k]);
    }
  }
  parts->kinetic = 0.5*(float) (kin / (parts->mass * parts->n));
  parts->energy_panda_nn = (float) (pot_tot->energy_panda_nn / parts->n);
  parts->energy_panda_np = (float) (pot_tot->energy_panda_np / parts->n);
  parts->energy_coul = (float) (pot_tot->energy_coul / parts->n);
  parts->energy_qcnm = (float) (pot_tot->energy_qcnm / parts->n);
  parts->energy_pauli = (float) (pot_tot->energy_pauli / parts->n);
  return 0;
}

int energia_sin_LUT(struct Particles *parts, struct Interaction *pot_tot){
	double kin = 0, rsq;
	reset_energies(pot_tot);
	int K_c = (int) ceil(pot_tot->coul->rcut / parts->l);
	int cantidad = ((2*K_c+1)*(2*K_c+1)*(2*K_c+1) + 1)/2;
	int i, j, mx, my, mz, idx, idxs[3*cantidad];
	int M = parts->M, M3 = M*M*M;
	construir_lista_celdas(K_c, idxs);
	for (int m = 0; m < M3; m++){
		mx = m/(M*M);
		my = (m/M) % M;
		mz = m % M;
		for (int k = 0; k < cantidad; k++){
			idx = (((mx+idxs[3*k]+M) % M)*M + (my+idxs[3*k+1]+M) % M)*M + (mz+idxs[3*k+2]+M) % M;
			i = parts->primero[m];
			if (idxs[3*k]<=1 && idxs[3*k+1]<=1 && idxs[3*k+2]<=1 && idxs[3*k]>=-1 && idxs[3*k+1]>=-1 && idxs[3*k+2]>=-1){
				while (i != -1){
					j = parts->primero[idx];
					while (j !=-1 && j != i){
						rsq = distancia(parts->q+3*i, parts->q+3*j, idxs+3*k, parts->l);
						interaction_sin_LUT(parts->type[i], parts->type[j], rsq, parts->p+3*i, parts->p+3*j, pot_tot);
						j = parts->siguiente[j];
					}
					i = parts->siguiente[i];
				}
			}else{
				while (i != -1){
					if(parts->type[i]/2 == 0){
						j = parts->primero[idx];
						while (j !=-1 && j != i){
							if (parts->type[j]/2 == 0){
								rsq = distancia(parts->q+3*i, parts->q+3*j, idxs+3*k, parts->l);
								interaction_sin_LUT(parts->type[i], parts->type[j], rsq, parts->p+3*i, parts->p+3*j, pot_tot);
							}
							j = parts->siguiente[j];
						}
					}
					i = parts->siguiente[i];
				}
			}
		}
	}
	// Energia cinetica
	for(int i = 0; i < parts->n; i++){
		for(int k = 0; k < 3; k++){
			kin = kin + (parts->p[3*i+k])*(parts->p[3*i+k]);
		}
	}
	parts->kinetic = 0.5*(float) (kin / (parts->mass * parts->n));
	parts->energy_panda_nn = (float) (pot_tot->energy_panda_nn / parts->n);
	parts->energy_panda_np = (float) (pot_tot->energy_panda_np / parts->n);
	parts->energy_coul = (float) (pot_tot->energy_coul / parts->n);
	parts->energy_qcnm = (float) (pot_tot->energy_qcnm / parts->n);
	parts->energy_pauli = (float) (pot_tot->energy_pauli / parts->n);
	return 0;
}

//----------------------- F/GUERZAS Y ENERGIA ----------------------------------

float fgorces_and_energy(struct Particles *parts, struct Interaction *pot_tot){
  double kin = 0, rsq;
  reset_energies(pot_tot);
  int M = parts->M, M3 = M*M*M, K_c = (int) ceil(pot_tot->coul->rcut / parts->l);
  int cantidad = ((2*K_c+1)*(2*K_c+1)*(2*K_c+1) + 1)/2;
  int i, j, mx, my, mz, idx, idxs[3*cantidad], neut_i, toca, primer_vecino = 0;
  construir_lista_celdas(K_c, idxs);
  float delta_q[3], delta_p[3], force_mod, gorce_mod = 0;
  for (int l = 0; l < 3*parts->n; l++){
    parts->forces[l] = 0;
    parts->gorces[l] = 0;
  }
  for (int m = 0; m < M3; m++){
    mx = m/(M*M);
    my = (m/M) % M;
    mz = m % M;
    i = parts->primero[m];
    for (int k = 0; k < cantidad; k++){
      idx = (((mx + idxs[3*k] + M)%M)*M + (my + idxs[3*k+1] + M)%M)*M + (mz + idxs[3*k+2] + M)%M;
			primer_vecino = (idxs[0]<=1 && idxs[1]<=1 && idxs[2]<=1 && idxs[0]>=-1 && idxs[1]>=-1 && idxs[2]>=-1);
      i = parts->primero[m];
      while (i != -1){
        neut_i = parts->type[i]/2;
        j = parts->primero[idx];
        while (j !=-1 && j != i){
          toca = primer_vecino || (!neut_i && (parts->type[j]/2 == 0));
          if(toca){
            rsq = 0;
            for(int l = 0; l < 3; l++){
              delta_q[l] = parts->q[3*i+l] - parts->q[3*j+l] - idxs[3*k+l]*parts->l;
              delta_p[l] = parts->p[3*i+l] - parts->p[3*j+l];
              rsq += delta_q[l]*delta_q[l];
            }
            interaction(parts->type[i], parts->type[j], rsq, parts->p+3*i, parts->p+3*j, pot_tot);
            force_mod = fgorce_mod(parts->type[i], parts->type[j], rsq, parts->p+3*i, parts->p+3*j, pot_tot, &gorce_mod);
            for (int l = 0; l < 3; l++){
              parts->forces[3*i+l] += force_mod*delta_q[l];
              parts->forces[3*j+l] -= force_mod*delta_q[l];
              parts->gorces[3*i+l] += gorce_mod*delta_p[l];
              parts->gorces[3*j+l] -= gorce_mod*delta_p[l];
            }
          }
          j = parts->siguiente[j];
        }
        i = parts->siguiente[i];
      }
    }
  }
  for(int i = 0; i < parts->n; i++){
    for(int k = 0; k < 3; k++){
      kin = kin + (parts->p[3*i+k])*(parts->p[3*i+k]);
    }
  }
  parts->kinetic = 0.5*(float)kin / (parts->mass * parts->n);
  parts->energy_panda_nn = (float) pot_tot->energy_panda_nn / parts->n;
  parts->energy_panda_np = (float) pot_tot->energy_panda_np / parts->n;
  parts->energy_coul = (float) pot_tot->energy_coul / parts->n;
  parts->energy_qcnm = (float) pot_tot->energy_qcnm / parts->n;
  parts->energy_pauli = (float) pot_tot->energy_pauli / parts->n;
  return 0;
}

float fgorces_and_energy_sin_LUT(struct Particles *parts, struct Interaction *pot_tot){
	double kin = 0, rsq;
	reset_energies(pot_tot);
	int M = parts->M, M3 = M*M*M, K_c = (int) ceil(pot_tot->coul->rcut / parts->l);
	int cantidad = ((2*K_c+1)*(2*K_c+1)*(2*K_c+1) + 1)/2;
	int i, j, mx, my, mz, idx, idxs[3*cantidad], neut_i, toca, primer_vecino = 0;
	construir_lista_celdas(K_c, idxs);
	float delta_q[3], delta_p[3], force_mod, gorce_mod = 0;
	for (int l = 0; l < 3*parts->n; l++){
		parts->forces[l] = 0;
		parts->gorces[l] = 0;
	}
	for (int m = 0; m < M3; m++){
		mx = m/(M*M);
		my = (m/M) % M;
		mz = m % M;
		i = parts->primero[m];
		for (int k = 0; k < cantidad; k++){
			idx = (((mx + idxs[3*k] + M)%M)*M + (my + idxs[3*k+1] + M)%M)*M + (mz + idxs[3*k+2] + M)%M;
			primer_vecino = (idxs[0]<=1 && idxs[1]<=1 && idxs[2]<=1 && idxs[0]>=-1 && idxs[1]>=-1 && idxs[2]>=-1);
			i = parts->primero[m];
			while (i != -1){
				neut_i = parts->type[i]/2;
				j = parts->primero[idx];
				while (j !=-1 && j != i){
					toca = primer_vecino || (!neut_i && (parts->type[j]/2 == 0));
					if(toca){
						rsq = 0;
						for(int l = 0; l < 3; l++){
							delta_q[l] = parts->q[3*i+l] - parts->q[3*j+l] - idxs[3*k+l]*parts->l;
							delta_p[l] = parts->p[3*i+l] - parts->p[3*j+l];
							rsq += delta_q[l]*delta_q[l];
						}
						interaction(parts->type[i], parts->type[j], rsq, parts->p+3*i, parts->p+3*j, pot_tot);
						force_mod = fgorce_mod_sin_LUT(parts->type[i], parts->type[j], rsq, parts->p+3*i, parts->p+3*j, pot_tot, &gorce_mod);
						for (int l = 0; l < 3; l++){
							parts->forces[3*i+l] += force_mod*delta_q[l];
							parts->forces[3*j+l] -= force_mod*delta_q[l];
							parts->gorces[3*i+l] += gorce_mod*delta_p[l];
							parts->gorces[3*j+l] -= gorce_mod*delta_p[l];
						}
					}
					j = parts->siguiente[j];
				}
				i = parts->siguiente[i];
			}
		}
	}
	for(int i = 0; i < parts->n; i++){
		for(int k = 0; k < 3; k++){
			kin = kin + (parts->p[3*i+k])*(parts->p[3*i+k]);
		}
	}
	parts->kinetic = 0.5*(float)kin / (parts->mass * parts->n);
	parts->energy_panda_nn = (float) pot_tot->energy_panda_nn / parts->n;
	parts->energy_panda_np = (float) pot_tot->energy_panda_np / parts->n;
	parts->energy_coul = (float) pot_tot->energy_coul / parts->n;
	parts->energy_qcnm = (float) pot_tot->energy_qcnm / parts->n;
	parts->energy_pauli = (float) pot_tot->energy_pauli / parts->n;
	return 0;
}
