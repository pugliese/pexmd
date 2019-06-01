#include "general.h"
#include "celda.h"
#include "interaccion.h"
#include "avanzar.h"
#include "inicializar.h"
#include "muestreo.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>


int main(int argc, char *argv[]){

  srand(1);
  float rho = 0.2;
  char opcion = 'g';
  if (argc >= 2){
    sscanf(argv[1], "%c\n", &opcion);
    if (argc >= 3){
      sscanf(argv[2], "%f\n", &rho);
    }
  }

  struct Particles parts;
// Potenciales
  struct Pauli pauli;
  pauli.qo = 1.644; // fm
  pauli.po = 120; // MeV/c
  pauli.D = 207; // MeV
  pauli.D = 0; // Por ahora sin Pauli
  pauli.scut2 = 10; // 0.7% del maximo
  pauli.shift = pauli.D*exp(-0.5*pauli.scut2);

  struct Panda_np panda_np;
  panda_np.mu_r = 1.7468; // fm^-1
  panda_np.mu_a = 1.6; // fm^-1
  panda_np.V_r = 3088.118; // MeV
  panda_np.V_a = 2666.647; // MeV
  panda_np.rcut = 5.4; // fm
  panda_np.rcut2 = 5.4*5.4; // fm^2
  panda_np.shift = 0;
  panda_np.shift = interaction_panda_np(panda_np.rcut, &panda_np);

  struct Panda_nn panda_nn;
  panda_nn.mu_o = 1.5; // fm^-1
  panda_nn.V_o = 373.118; // MeV
  panda_nn.rcut = 5.4; // fm
  panda_nn.rcut2 = 5.4*5.4; // fm^2
  panda_nn.shift = 0;
  panda_nn.shift = interaction_panda_nn(panda_nn.rcut, &panda_nn);

  struct Interaction pot_tot;
  pot_tot.pauli = &pauli;
  pot_tot.panda_nn = &panda_nn;
  pot_tot.panda_np = &panda_np;
	pot_tot.rcut = max(sqrt(pauli.scut2)*pauli.qo, panda_nn.rcut);
  int N_LUT = 100000;
  build_LUTs(&pot_tot, N_LUT, N_LUT, N_LUT);

// Parametros termodinamicos
  struct Externos params;
  int N = 20;
  params.L = N/pow(rho, 1.0/3); // mayor a 2*qo*scut, rho en fm^-3
  N = N*N*N;
  params.T = 4; // MeV



// MUESTREO EN TEMPERATURA
  if(opcion=='m'){
    // Particulas
    float mass = 938; // MeV/c^2
    int comps[4] = {N/4, N/4, N/4, N/4};
    inicializar(&parts, comps, 4, mass, params.L, params.T, &pot_tot);
    energia(&parts, &pot_tot);
    params.delta_q = 0.5*sqrt(params.T/4)*sqrt(0.05/rho); // fm
    params.delta_p = 0.5*sqrt(parts.mass*params.T); // MeV/c
    char filename_termo[255], filename_config[255];
    int factor_term = 14000, factor_desc = 10, Nsamp = 100;
    int append = 0; // Solo un frame para lammps
    N_steps(&parts, &pot_tot, &params, parts.n*factor_term);
    float dT = 0.1;
    sprintf(filename_termo, "data_sin_pauli/termo_%.3f.txt", rho);

    FILE* fp = fopen(filename_termo, "w");
    fclose(fp);
    for (int t = 0; t < 40; t++){
      energia(&parts, &pot_tot);
      params.delta_q = 0.5*sqrt(params.T/4)*sqrt(0.05/rho); // fm
      params.delta_p = 0.5*sqrt(parts.mass*params.T); // MeV/c
      printf("rho = %.2f fm^-3 | T = %.2f MeV\n", rho, params.T);
      muestrear_termo(filename_termo, &parts, &pot_tot, &params, Nsamp, factor_desc, factor_term);
      sprintf(filename_config, "data_sin_pauli/config_%.3f_%.3f.lammpstrj", rho, params.T);
      save_lammpstrj(filename_config, &parts, params.L, append);
      params.T -= dT;
    }
  }


// MUESTREO EN TEMPERATURA - 1000 particulas
  if(opcion=='n'){
    // Particulas
    N = 10;
    params.L = N/pow(rho, 1.0/3); // mayor a 2*qo*scut, rho en fm^-3
    N = N*N*N;
    float mass = 938; // MeV/c^2
    int comps[4] = {N/4, N/4, N/4, N/4};
    inicializar(&parts, comps, 4, mass, params.L, params.T, &pot_tot);
    energia(&parts, &pot_tot);
    params.delta_q = 0.5*pow(0.05/rho, 0.5); // fm
    params.delta_p = 0.5*sqrt(parts.mass*params.T); // MeV/c
    char filename_termo[255], filename_config[255];
    int factor_term = 25000, factor_desc = 50, Nsamp = 100;
    int append = 0; // Solo un frame para lammps
    N_steps(&parts, &pot_tot, &params, parts.n*factor_term);
    float dT = 0.1;
    sprintf(filename_termo, "data_sin_pauli_1000/termo_%.3f.txt", rho);

    FILE* fp = fopen(filename_termo, "w");
    fclose(fp);
    for (int t = 0; t < 40; t++){
      energia(&parts, &pot_tot);
      params.delta_q = 0.5*sqrt(params.T/4)*pow(0.05/rho, 0.5); // fm
      params.delta_p = 0.5*sqrt(parts.mass*params.T); // MeV/c
      printf("rho = %.2f fm^-3 | T = %.2f MeV\n", rho, params.T);
      muestrear_termo(filename_termo, &parts, &pot_tot, &params, Nsamp, factor_desc, factor_term);
      sprintf(filename_config, "data_sin_pauli_1000/config_%.3f_%.3f.lammpstrj", rho, params.T);
      save_lammpstrj(filename_config, &parts, params.L, append);
      params.T -= dT;
    }
  }

// MUESTREO DE UNICA TEMPERATURA - 1000 particulas con paso chico
  if(opcion=='c'){
    printf("Opcion paso chico\n");
    int factor_term = 3000, t;
    params.T = 0.55; // MeV
    float mass = 938; // MeV/c^2
    int comps[4] = {N/4, N/4, N/4, N/4};
    inicializar(&parts, comps, 4, mass, params.L, params.T, &pot_tot);
    energia(&parts, &pot_tot);
    params.delta_q = 0.25*pow(params.T/5, 0.5)*pow(0.05/rho, 0.5); // fm
    params.delta_p = 0.5*sqrt(parts.mass*params.T); // MeV/c
    char filename_config[255];
    sprintf(filename_config, "data_sin_pauli_1000_delta_chico/config_%.3f_%.3f.lammpstrj", rho, params.T);
    t = time(NULL);
    long int aceptados = 0;
    for(int i = 0; i < 100; i++){
      aceptados += N_steps(&parts, &pot_tot, &params, parts.n*factor_term);
      printf("\rProgreso: %d%%", i+1);
      fflush(stdout);
    }
    printf("\n%ld%% en %d segundos\n", (aceptados)/(parts.n*factor_term), (int)time(NULL)-t);
    save_lammpstrj(filename_config, &parts, params.L, 0);
  }
  if(opcion=='g'){
    printf("Opcion paso grande\n");
    int factor_term = 3000, t;
    params.T = 0.55; // MeV
    float mass = 938; // MeV/c^2
    int comps[4] = {N/4, N/4, N/4, N/4};
    inicializar(&parts, comps, 4, mass, params.L, params.T, &pot_tot);
    energia(&parts, &pot_tot);
    params.delta_q = 1*pow(params.T/5, 0.5)*pow(0.05/rho, 0.5); // fm
    params.delta_p = 0.5*sqrt(parts.mass*params.T); // MeV/c
    char filename_config[255];
    sprintf(filename_config, "data_sin_pauli_1000_delta_grande/config_%.3f_%.3f.lammpstrj", rho, params.T);
    t = time(NULL);
    long int aceptados = 0;
    for(int i = 0; i < 100; i++){
      printf("\rProgreso: %d%%", i);
      fflush(stdout);
      aceptados += N_steps(&parts, &pot_tot, &params, parts.n*factor_term);
    }
    printf("%ld%% en %d segundos\n", (aceptados)/(parts.n*factor_term), (int)time(NULL)-t);
    save_lammpstrj(filename_config, &parts, params.L, 0);
  }
  if(opcion=='s'){
    printf("Opcion simple cubica con/sin Pauli\n");
    float rho_min = 0.100, rho_max = 0.250, delta_rho = 2*0.005;
    float a_max = 15.0, a_min = 10.0, delta_a = 0.25;
    float b_max = 0.5, b_min = 0, delta_b = 0.025;

    params.T = 0.0; // MeV
    float mass = 938; // MeV/c^2
    int comps[4] = {N/4, N/4, N/4, N/4};
    char filename[255];
    params.L = pow(N/rho_min, 1.0/3); // mayor a 2*qo*scut, rho en fm^-3
    inicializar(&parts, comps, 4, mass, params.L, params.T, &pot_tot);

    pot_tot.pauli->D = 207; // MeV
    liberar_LUTs(&pot_tot);
    build_LUTs(&pot_tot, N_LUT, N_LUT, N_LUT);
    energia(&parts, &pot_tot);

    int idx_ref = 1000;
    float ref_np = pot_tot.panda_np->LUT[idx_ref], factor_np;
    float ref_nn = pot_tot.panda_nn->LUT[idx_ref], factor_nn;
    int frame = 0;
    for(rho = rho_min; rho <= rho_max; rho += delta_rho){
      printf("\nrho = %.3f\n", rho);
      params.L = pow(N/rho, 1.0/3); // mayor a 2*qo*scut, rho en fm^-3
      set_box(&parts, pot_tot.rcut, params.L);
      save_lammpstrj("SC/configs.lammpstrj", &parts, params.L, frame);
      frame++;
      factor_np = ref_np/pot_tot.panda_np->LUT[idx_ref];
      for(int i = 0; i < N_LUT; i++)  pot_tot.panda_np->LUT[i] *= factor_np;
      factor_nn = ref_nn/pot_tot.panda_nn->LUT[idx_ref];
      for(int j = 0; j < N_LUT; j++)  pot_tot.panda_nn->LUT[j] *= factor_nn;
      energia(&parts, &pot_tot);
      printf("%f %f\n", parts.energy_panda/N, parts.energy_pauli/N);
      sprintf(filename, "SC/data_%.3f.txt", rho);
      FILE *fp = fopen(filename, "w");
      fprintf(fp, "%f ", parts.energy_panda/N);
      for(float b = b_max; b >= b_min; b -= delta_b)  fprintf(fp, " %f", b);
      for(float a = a_min; a <= a_max; a += delta_a){
        factor_np = a*ref_np/pot_tot.panda_np->LUT[idx_ref];
        for(int i = 0; i < N_LUT; i++)  pot_tot.panda_np->LUT[i] *= factor_np;
        fprintf(fp, "\n%f", a);
        for(float b = b_max; b >= b_min; b -= delta_b){
          factor_nn = b*ref_nn/pot_tot.panda_nn->LUT[idx_ref];
          for(int j = 0; j < N_LUT; j++)  pot_tot.panda_nn->LUT[j] *= factor_nn;
          printf("a: %d%%  |||  b: %d%%\r", (int) floor(100*(a-a_min)/(a_max-a_min)), (int) floor(100*(b_max-b)/(b_max-b_min)));
          fflush(stdout);
          energia(&parts, &pot_tot);
          fprintf(fp, " %f", (parts.energy_panda+parts.energy_pauli)/N);
        }
      }
      fclose(fp);
    }
  }
  liberar(&parts);
  liberar_LUTs(&pot_tot);
  return 0;
}
