#include "general.h"
#include "celda.h"
#include "interaccion.h"
#include "inicializar.h"
#include "avanzar.h"
#include "virial.h"
#include "muestreo.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

int muestrear_termo(char *filename, struct Particles *parts, struct Interaction *pot_tot, struct Externos *params, int Nsamp, int factor_desc, int factor_term){
	int aceptados = 0, segs = time(NULL);
	// Termalizacion
	int N_tot = Nsamp*factor_desc + factor_term;
	printf("\rProgreso: 0%%");
	N_steps(parts, pot_tot, params, parts->n*factor_term);
	printf("\rProgreso: %d%%", (100*factor_term)/N_tot);
	// MUESTREO
	float ekin = 0, enuc = 0, epauli = 0, PyT_tot[2], PyT[2];
	PyT_tot[0] = 0, PyT_tot[1] = 0;
	for (int i = 0; i < Nsamp; i++){
		aceptados += N_steps(parts, pot_tot, params, parts->n*factor_desc);
		printf("\rProgreso: %d%%", (100*(factor_term + i*factor_desc))/N_tot);
		fflush(stdout);
		energia(parts, pot_tot);
		ekin += parts->kinetic/(parts->n*Nsamp);
		enuc += parts->energy_panda/(parts->n*Nsamp);
		epauli += parts->energy_pauli/(parts->n*Nsamp);
		presion_temperatura(parts, pot_tot, params, PyT);
		PyT_tot[0] += PyT[0]/Nsamp;
		PyT_tot[1] += PyT[1]/Nsamp;
	}
	printf("\n");
	segs = time(NULL) - segs;
	FILE* fp = fopen(filename, "a");
	fprintf(fp, "%f %f %f %f %f %f %d %.2f\n", params->T, ekin, enuc, epauli, PyT_tot[0], PyT_tot[1], segs, 100.0*aceptados/(Nsamp*factor_desc*parts->n));
  fclose(fp);
	return aceptados;
}

int muestrear_termo_QCNM(char *filename, struct Particles *parts, struct Interaction *pot_tot, struct Externos *params, int Nsamp, int factor_desc, int factor_term){
	int aceptados = 0, segs = time(NULL);
	// Termalizacion
	int N_tot = Nsamp*factor_desc + factor_term;
	printf("\rProgreso: 0%%");
	N_steps_QCNM(parts, pot_tot, params, parts->n*factor_term);
	printf("\rProgreso: %d%%", (100*factor_term)/N_tot);
	// MUESTREO
	float ekin = 0, enuc = 0, epauli = 0, PyT_tot[2], PyT[2];
	PyT_tot[0] = 0, PyT_tot[1] = 0;
	for (int i = 0; i < Nsamp; i++){
		aceptados += N_steps_QCNM(parts, pot_tot, params, parts->n*factor_desc);
		printf("\rProgreso: %d%%", (100*(factor_term + i*factor_desc))/N_tot);
		fflush(stdout);
		//energia_QCNM(parts, pot_tot);
		ekin += parts->kinetic/(parts->n*Nsamp);
		enuc += parts->energy_panda/(parts->n*Nsamp);
		epauli += parts->energy_pauli/(parts->n*Nsamp);
		presion_temperatura_QCNM(parts, pot_tot, params, PyT);
		PyT_tot[0] += PyT[0]/Nsamp;
		PyT_tot[1] += PyT[1]/Nsamp;
	}
	printf("\n");
	segs = time(NULL) - segs;
	FILE* fp = fopen(filename, "a");
	fprintf(fp, "%f %f %f %f %f %f %d %.2f\n", params->T, ekin, enuc, epauli, PyT_tot[0], PyT_tot[1], segs, 100.0*aceptados/(Nsamp*factor_desc*parts->n));
  fclose(fp);
	return aceptados;
}
