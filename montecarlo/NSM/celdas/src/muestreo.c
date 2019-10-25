#include "general.h"
#include "celda.h"
#include "interaccion.h"
#include "energia.h"
#include "avanzar.h"
#include "muestreo.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

float presion_temperatura(struct Particles *parts, struct Interaction *pot_tot, float *PyT){
  int M = parts->M, M3 = M*M*M, K_c = (int) ceil(pot_tot->coul->rcut / parts->l);
  int cantidad = ((2*K_c+1)*(2*K_c+1)*(2*K_c+1) + 1)/2;
  int i, j, mx, my, mz, idx, idxs[3*cantidad], neut_i, toca, primer_vecino = 0;
  construir_lista_celdas(K_c, idxs);
  float delta_q[3], delta_p[3], force_mod, gorce_mod = 0, rsq, qF = 0, qp = 0;
  for (int m = 0; m < M3; m++){
    mx = m/(M*M);
    my = (m/M) % M;
    mz = m % M;
    for (int k = 0; k < cantidad; k++){
      idx = (((mx + idxs[3*k] + M)%M)*M + (my + idxs[3*k+1] + M)%M)*M + (mz + idxs[3*k+2] + M)%M;
      primer_vecino = (idxs[3*k]<=1 && idxs[3*k+1]<=1 && idxs[3*k+2]<=1 && idxs[3*k+0]>=-1 && idxs[3*k+1]>=-1 && idxs[3*k+2]>=-1);
      i = parts->primero[m];
      while (i != -1){
        neut_i = parts->type[i]/2;
        j = parts->primero[idx];
        while (j !=-1 && j != i){
          toca = primer_vecino || (!neut_i && parts->type[j]/2 == 0);
          if(toca){
            rsq = 0;
            for(int l = 0; l < 3; l++){
              delta_q[l] = parts->q[3*i+l] - parts->q[3*j+l] - idxs[3*k+l]*parts->l;
              delta_p[l] = parts->p[3*i+l] - parts->p[3*j+l];
              rsq += delta_q[l]*delta_q[l];
            }
            force_mod = fgorce_mod(parts->type[i], parts->type[j], rsq, parts->p+3*i, parts->p+3*j, pot_tot, &gorce_mod);
            for (int l = 0; l < 3; l++){
              qF += force_mod*delta_q[l]*delta_q[l];
              qp -= gorce_mod*delta_p[l]*delta_p[l];
            }
          }
          j = parts->siguiente[j];
        }
        i = parts->siguiente[i];
      }
    }
  }
  for (int l = 0; l < 3*parts->n; l++){
    qp += parts->p[l]*parts->p[l]/parts->mass;
  }
	float L = parts->M*parts->l;
  PyT[0] = (qF+qp)/(3*L*L*L);
  PyT[1] = qp/(3*parts->n);
  return 0;
}

float presion_temperatura_sin_LUT(struct Particles *parts, struct Interaction *pot_tot, float *PyT){
  int M = parts->M, M3 = M*M*M, K_c = (int) ceil(pot_tot->coul->rcut / parts->l);
  int cantidad = ((2*K_c+1)*(2*K_c+1)*(2*K_c+1) + 1)/2;
  int i, j, mx, my, mz, idx, idxs[3*cantidad], neut_i = 0, toca, primer_vecino = 0;
  construir_lista_celdas(K_c, idxs);
  float delta_q[3], delta_p[3], force_mod, gorce_mod = 0, rsq, qF = 0, qp = 0;
  for (int m = 0; m < M3; m++){
    mx = m/(M*M);
    my = (m/M) % M;
    mz = m % M;
    for (int k = 0; k < cantidad; k++){
      idx = (((mx + idxs[3*k] + M)%M)*M + (my + idxs[3*k+1] + M)%M)*M + (mz + idxs[3*k+2] + M)%M;
      primer_vecino = (idxs[3*k]<=1 && idxs[3*k+1]<=1 && idxs[3*k+2]<=1 && idxs[3*k+0]>=-1 && idxs[3*k+1]>=-1 && idxs[3*k+2]>=-1);
      i = parts->primero[m];
      while (i != -1){
        neut_i = parts->type[i]/2;
        j = parts->primero[idx];
        while (j !=-1 && j != i){
          toca = primer_vecino || (!neut_i && parts->type[j]/2 == 0);
          if(toca){
            rsq = 0;
            for(int l = 0; l < 3; l++){
              delta_q[l] = parts->q[3*i+l] - parts->q[3*j+l] - idxs[3*k+l]*parts->l;
              delta_p[l] = parts->p[3*i+l] - parts->p[3*j+l];
              rsq += delta_q[l]*delta_q[l];
            }
            force_mod = fgorce_mod_sin_LUT(parts->type[i], parts->type[j], rsq, parts->p+3*i, parts->p+3*j, pot_tot, &gorce_mod);
            for (int l = 0; l < 3; l++){
              qF += force_mod*delta_q[l]*delta_q[l];
              qp -= gorce_mod*delta_p[l]*delta_p[l];
            }
          }
          j = parts->siguiente[j];
        }
        i = parts->siguiente[i];
      }
    }
  }
  for (int l = 0; l < 3*parts->n; l++){
    qp += parts->p[l]*parts->p[l]/parts->mass;
  }
	float L = parts->M*parts->l;
	PyT[0] = (qF+qp)/(3*L*L*L);
  PyT[1] = qp/(3*parts->n);
  return 0;
}

int muestrear_termo(char *filename, struct Particles *parts, struct Interaction *pot_tot, struct Externos *params, int Nsamp, int factor_desc, int factor_term){
	int aceptados = 0, segs = time(NULL);
	// Termalizacion
	int N_tot = Nsamp*factor_desc + factor_term;
	printf("\rProgreso: 0.00%%");
	fflush(stdout);
	N_steps(parts, pot_tot, params, parts->n*factor_term);
	printf("\rProgreso: %.2f%%", (100.0*factor_term)/N_tot);
	fflush(stdout);
	// MUESTREO
	float ekin = 0, enn = 0, enp = 0, ecoul = 0, eqcnm = 0, epauli = 0, PyT_tot[2], PyT[2];
	PyT_tot[0] = 0, PyT_tot[1] = 0;
	for (int i = 0; i < Nsamp; i++){
		aceptados += N_steps(parts, pot_tot, params, parts->n*factor_desc);
		printf("\rProgreso: %.2f%%", (100.0*(factor_term + (i+1)*factor_desc))/N_tot);
		fflush(stdout);
		energia(parts, pot_tot);
		ekin += parts->kinetic / Nsamp;
		enn += parts->energy_panda_nn / Nsamp;
		enp += parts->energy_panda_np / Nsamp;
		ecoul += parts->energy_coul / Nsamp;
		eqcnm += parts->energy_qcnm / Nsamp;
		epauli += parts->energy_pauli / Nsamp;
		presion_temperatura_sin_LUT(parts, pot_tot, PyT);
		PyT_tot[0] += PyT[0]/Nsamp;
		PyT_tot[1] += PyT[1]/Nsamp;
	}
	printf("\n");
	segs = time(NULL) - segs;
	FILE* fp = fopen(filename, "a");
	fprintf(fp, "%f %f %f %f %f %f %f %f %f %d %.2f\n", params->T, ekin, enp, enn, eqcnm, ecoul, epauli, PyT_tot[0], PyT_tot[1], segs, 100.0*aceptados/(Nsamp*factor_desc*parts->n));
  fclose(fp);
	return aceptados;
}

float mean_radius(struct Particles *parts){
  float CM[3] = {0.0, 0.0, 0.0};
  for(int i = 0; i < parts->n; i++){
    for(int k = 0; k < 3; k++){
      CM[k] += parts->q[3*i+k]/parts->n;
    }
  }
  double R = 0;
  for(int i = 0; i < parts->n; i++){
    for(int k = 0; k < 3; k++){
      R += (CM[k] - parts->q[3*i+k])*(CM[k] - parts->q[3*i+k])/parts->n;
    }
  }
  return (float) R;
}

int muestrear_nucleo(char *filename, struct Particles *parts, struct Interaction *pot_tot, struct Externos *params, int Nsamp, int factor_desc, int factor_term){
	int aceptados = 0, segs = time(NULL);
	// Termalizacion
	//int N_tot = Nsamp*factor_desc + factor_term;
	//printf("\rProgreso: 0.00%%");
	//fflush(stdout);
	N_steps(parts, pot_tot, params, parts->n*factor_term);
	//printf("\rProgreso: %.2f%%", (100.0*factor_term)/N_tot);
	//fflush(stdout);
	// MUESTREO
	float ekin = 0, enn = 0, enp = 0, ecoul = 0, eqcnm = 0, epauli = 0, T = 0, R = 0, PyT[2];
	for (int i = 0; i < Nsamp; i++){
		aceptados += N_steps(parts, pot_tot, params, parts->n*factor_desc);
		//printf("\rProgreso: %.2f%%", (100.0*(factor_term + (i+1)*factor_desc))/N_tot);
		//fflush(stdout);
		energia(parts, pot_tot);
		ekin += parts->kinetic / Nsamp;
		enn += parts->energy_panda_nn / Nsamp;
		enp += parts->energy_panda_np / Nsamp;
		ecoul += parts->energy_coul / Nsamp;
		eqcnm += parts->energy_qcnm / Nsamp;
		epauli += parts->energy_pauli / Nsamp;
    R += mean_radius(parts) / Nsamp;
		presion_temperatura_sin_LUT(parts, pot_tot, PyT);
		T += PyT[1]/Nsamp;
	}
	//printf("\n");
	segs = time(NULL) - segs;
	FILE* fp = fopen(filename, "a");
	fprintf(fp, "%f %f %f %f %f %f %f %f %f %d %.2f\n", params->T, ekin, enp, enn, eqcnm, ecoul, epauli, sqrt(R), T, segs, 100.0*aceptados/(Nsamp*factor_desc*parts->n));
  fclose(fp);
	return aceptados;
}

float estimate_acceptance(struct Particles *parts, struct Interaction *pot_tot, struct Externos *params){
  int aceptados = N_steps(parts, pot_tot, params, 10*parts->n);
  return 1.0*aceptados/(10*parts->n);
}

float ajustar_deltas(struct Particles *parts, struct Interaction *pot_tot, struct Externos *params, float acept_min, float acept_max){
  float aceptacion = estimate_acceptance(parts, pot_tot, params);
  float prev_delta_q = params->delta_q, prev_delta_p = params->delta_p, f_max = 1, f_min = 1, f = 1;
  while(aceptacion < acept_min){
    f_min = f_min/2;
    params->delta_q = f_min*prev_delta_q;
    params->delta_p = f_min*prev_delta_p;
    aceptacion = estimate_acceptance(parts, pot_tot, params);
  }
  if(aceptacion <= acept_max) return f_min;
  //aceptacion(f_min) > acept_max
  while(aceptacion > acept_max){
    f_max = f_max*2;
    params->delta_q = f_max*prev_delta_q;
    params->delta_p = f_max*prev_delta_p;
    aceptacion = estimate_acceptance(parts, pot_tot, params);
  }
  if(aceptacion >= acept_min) return f_max;
  //aceptacion(f_max) < acept_min
  //aceptacion(f_min) > acept_max
  while(aceptacion < acept_min || aceptacion > acept_max){
    f = (f_min+f_max)/2;
    params->delta_q = f_max*prev_delta_q;
    params->delta_p = f_max*prev_delta_p;
    aceptacion = estimate_acceptance(parts, pot_tot, params);
    if(aceptacion < acept_min){
      f_max = f;
    }else{
      f_min = f;
    }
  }
	return f;
}

int ajustar_deltas_full(struct Particles *parts, struct Interaction *pot_tot, struct Externos *params, float acept_min, float acept_max){
  float prev_delta_p = params->delta_p;
  params->delta_p = 0;
  ajustar_deltas(parts, pot_tot, params, acept_min, acept_max);
  float prev_delta_q = params->delta_q;
  params->delta_p = prev_delta_p;
  params->delta_q = 0;
  ajustar_deltas(parts, pot_tot, params, acept_min, acept_max);
  params->delta_q = prev_delta_q;
  ajustar_deltas(parts, pot_tot, params, acept_min, acept_max);
	return 0;
}
