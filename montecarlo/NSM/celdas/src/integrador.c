#include "general.h"
#include "celda.h"
#include "interaccion.h"
#include "energia.h"
#include "integrador.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

//-------------------------------- VERLET --------------------------------------

int verlet_first_step(struct Particles *parts, float dt){
	return 0;
}

int verlet_second_step(struct Particles *parts, float dt){
	return 0;
}

int velocity_verlet(struct Particles *parts, struct Interaction *pot_tot, float dt){
	verlet_first_step(parts, dt);
	ajustar_celdas(parts);
	fgorces_and_energy(parts, pot_tot);
	verlet_second_step(parts, dt);
	ajustar_celdas(parts);
	fgorces_and_energy(parts, pot_tot);
	return 0;
}

//---------------------------- RUNGE-KUTTA 2 -----------------------------------

int RK2_half_step(struct Particles *parts, float dt){
	float m = parts->mass;
	for(int i = 0; i < 3*parts->n; i++){
		parts->q[i] += 0.5*(parts->p[i]/m - parts->gorces[i])*dt;
		parts->p[i] += 0.5*parts->forces[i]*dt;
	}
	return 0;
}

int RK2(struct Particles *parts, struct Interaction *pot_tot, float dt){
	RK2_half_step(parts, dt);
	ajustar_celdas(parts);
	fgorces_and_energy(parts, pot_tot);
	RK2_half_step(parts, dt);
	ajustar_celdas(parts);
	fgorces_and_energy(parts, pot_tot);
	return 0;
}

//---------------------------- MIDPOINT RULE -----------------------------------

int MidPoint_Rule_iteration(struct Particles *parts, float *qo, float *po, float dt){
	float m = parts->mass;
	for(int i = 0; i < 3*parts->n; i++){
		parts->q[i] = qo[i] + 0.5*(parts->p[i]/m - parts->gorces[i])*dt;
		parts->p[i] = po[i] + 0.5*parts->forces[i]*dt;
	}
	return 0;
}

int MidPoint_Rule(struct Particles *parts, struct Interaction *pot_tot, float dt, int K){
	float qo[3*parts->n];
	float po[3*parts->n];
	for(int i = 0; i < 3*parts->n; i++){
		qo[i] = parts->q[i];
		po[i] = parts->p[i];
	}
	for(int k = 0; k < K; k++){
		MidPoint_Rule_iteration(parts, qo, po, dt);
		ajustar_celdas(parts);
		fgorces_and_energy(parts, pot_tot);
	}
	for(int i = 0; i < 3*parts->n; i++){
		parts->q[i] = 2*parts->q[i] - qo[i];
		parts->p[i] = 2*parts->p[i] - po[i];
	}
	ajustar_celdas(parts);
	fgorces_and_energy(parts, pot_tot);
	return 0;
}
