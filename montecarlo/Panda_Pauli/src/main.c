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
  if (argc >= 2){
    sscanf(argv[1], "%f\n", &rho);
  }

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
  int N_LUT = 10000;
  build_LUTs(&pot_tot, N_LUT, N_LUT, N_LUT);

// Parametros termodinamicos
  struct Externos params;
  int N = 20;
  params.L = N/pow(rho, 1.0/3); // mayor a 2*qo*scut, rho en fm^-3
  N = N*N*N;
  params.T = 4; // MeV


/*
// MUESTREO EN TEMPERATURA
  // Particulas
  struct Particles parts;
  float mass = 938; // MeV/c^2
  int comps[4] = {N/4, N/4, N/4, N/4};
  inicializar(&parts, comps, 4, mass, params.L, params.T, &pot_tot);
  energia(&parts, &pot_tot);
  params.delta_q = 0.7*pow(0.05/rho, 0.5); // fm
  params.delta_p = 0.5*sqrt(parts.mass*params.T); // MeV/c
  char filename_termo[255], filename_config[255];
  int factor_term = 20000, factor_desc = 50, Nsamp = 100;
  int append = 0; // Solo un frame para lammps
  N_steps(&parts, &pot_tot, &params, parts.n*factor_term);
  float dT = 0.1;
  sprintf(filename_termo, "data_sin_pauli/termo_%.3f.txt", rho);

  FILE* fp = fopen(filename_termo, "w");
  fclose(fp);
  for (int t = 0; t < 40; t++){
    energia(&parts, &pot_tot);
    params.delta_q = 0.7*pow(params.T/5, 0.45)*pow(0.05/rho, 0.5); // fm
    params.delta_p = 0.5*sqrt(parts.mass*params.T); // MeV/c
    printf("rho = %.2f fm^-3 | T = %.2f MeV\n", rho, params.T);
    muestrear_termo(filename_termo, &parts, &pot_tot, &params, Nsamp, factor_desc, factor_term);
    sprintf(filename_config, "data_sin_pauli/config_%.3f_%.3f.lammpstrj", rho, params.T);
    save_lammpstrj(filename_config, &parts, params.L, append);
    params.T -= dT;
  }
*/

/*
// MUESTREO EN TEMPERATURA - 1000 particulas
  // Particulas
  N = 10;
  params.L = N/pow(rho, 1.0/3); // mayor a 2*qo*scut, rho en fm^-3
  N = N*N*N;
  struct Particles parts;
  float mass = 938; // MeV/c^2
  int comps[4] = {N/4, N/4, N/4, N/4};
  inicializar(&parts, comps, 4, mass, params.L, params.T, &pot_tot);
  energia(&parts, &pot_tot);
  params.delta_q = 0.7*pow(0.05/rho, 0.5); // fm
  params.delta_p = 0.5*sqrt(parts.mass*params.T); // MeV/c
  char filename_termo[255], filename_config[255];
  int factor_term = 140000, factor_desc = 50, Nsamp = 200;
  int append = 0; // Solo un frame para lammps
  N_steps(&parts, &pot_tot, &params, parts.n*factor_term);
  float dT = 0.1;
  sprintf(filename_termo, "data_sin_pauli_1000/termo_%.3f.txt", rho);

  FILE* fp = fopen(filename_termo, "w");
  fclose(fp);
  for (int t = 0; t < 40; t++){
    energia(&parts, &pot_tot);
    params.delta_q = 0.7*sqrt(params.T/5)*pow(0.05/rho, 0.5); // fm
    params.delta_p = 0.5*sqrt(parts.mass*params.T); // MeV/c
    printf("rho = %.2f fm^-3 | T = %.2f MeV\n", rho, params.T);
    muestrear_termo(filename_termo, &parts, &pot_tot, &params, Nsamp, factor_desc, factor_term);
    sprintf(filename_config, "data_sin_pauli_1000/config_%.3f_%.3f.lammpstrj", rho, params.T);
    save_lammpstrj(filename_config, &parts, params.L, append);
    params.T -= dT;
  }
*/
/*
// MUESTREO EN TEMPERATURA - 1000 particulas - paso chico con checkpoint
  // Particulas
  N = 10;
  N = N*N*N;
  struct Particles parts;
  parts.mass = 938; // MeV/c^2
  load_lammpstrj("data_sin_pauli_1000/config_0.050_2.500.lammpstrj", &parts, &params.L, pot_tot.rcut);
  save_lammpstrj("que_mierda_cargo.lammpstrj", &parts, params.L, 0);
  for (int i = 0; i < N; i++){
    printf("%.2f ", parts.q[3*i+1]);
  }
  params.delta_q = 0.5*pow(0.05/rho, 0.5); // fm
  params.delta_p = 0.5*sqrt(parts.mass*params.T); // MeV/c
  char filename_termo[255], filename_config[255];
  int factor_term = 24000, factor_desc = 10, Nsamp = 100;
  int append = 0; // Solo un frame para lammps
  float dT = 0.1;
  rho = parts.n/(params.L*params.L*params.L);
  sprintf(filename_termo, "data_sin_pauli_1000_delta_chico/termo_%.3f.txt", rho);

  FILE* fp = fopen(filename_termo, "w");
  fclose(fp);
  params.T = 2.4; // MeV
  for (int t = 0; t < 24; t++){
    energia(&parts, &pot_tot);
    //params.delta_q = 0.5*sqrt(params.T/4)*pow(0.05/rho, 0.5); // fm
    params.delta_q = 0.25*sqrt(params.T/4)*pow(0.05/rho, 0.5); // fm
    params.delta_p = 0.5*sqrt(parts.mass*params.T); // MeV/c
    printf("rho = %.2f fm^-3 | T = %.2f MeV\n", rho, params.T);
    muestrear_termo(filename_termo, &parts, &pot_tot, &params, Nsamp, factor_desc, factor_term);
    sprintf(filename_config, "data_sin_pauli_1000_delta_chico/config_%.3f_%.3f.lammpstrj", rho, params.T);
    save_lammpstrj(filename_config, &parts, params.L, append);
    params.T -= dT;
  }
*/

// MUESTREO EN TEMPERATURA - 1000 particulas - paso grande con checkpoint
  // Particulas
  N = 10;
  N = N*N*N;
  struct Particles parts;
  parts.mass = 938; // MeV/c^2
  load_lammpstrj("data_sin_pauli_1000/config_0.050_2.500.lammpstrj", &parts, &params.L, pot_tot.rcut);
  save_lammpstrj("que_mierda_cargo.lammpstrj", &parts, params.L, 0);
  for (int i = 0; i < N; i++){
    printf("%.2f ", parts.q[3*i+1]);
  }
  params.delta_q = 0.5*pow(0.05/rho, 0.5); // fm
  params.delta_p = 0.5*sqrt(parts.mass*params.T); // MeV/c
  char filename_termo[255], filename_config[255];
  int factor_term = 50000, factor_desc = 50, Nsamp = 100;
  int append = 0; // Solo un frame para lammps
  float dT = 0.1;
  rho = parts.n/(params.L*params.L*params.L);
  sprintf(filename_termo, "data_sin_pauli_1000_delta_grande/termo_%.3f.txt", rho);

  FILE* fp = fopen(filename_termo, "w");
  fclose(fp);
  params.T = 2.4; // MeV
  for (int t = 0; t < 24; t++){
    energia(&parts, &pot_tot);
    //params.delta_q = 0.5*sqrt(params.T/4)*pow(0.05/rho, 0.5); // fm
    params.delta_q = 1.0*sqrt(params.T/4)*pow(0.05/rho, 0.5); // fm
    params.delta_p = 0.5*sqrt(parts.mass*params.T); // MeV/c
    printf("rho = %.2f fm^-3 | T = %.2f MeV\n", rho, params.T);
    muestrear_termo(filename_termo, &parts, &pot_tot, &params, Nsamp, factor_desc, factor_term);
    sprintf(filename_config, "data_sin_pauli_1000_delta_grande/config_%.3f_%.3f.lammpstrj", rho, params.T);
    save_lammpstrj(filename_config, &parts, params.L, append);
    params.T -= dT;
  }

/*
// MUESTREO DE UNICA TEMPERATURA
  int factor_term = 3000, t;
  params.T = 0.55; // MeV
  struct Particles parts;
  float mass = 938; // MeV/c^2
  int comps[4] = {N/4, N/4, N/4, N/4};
  inicializar(&parts, comps, 4, mass, params.L, params.T, &pot_tot);
  energia(&parts, &pot_tot);
  params.delta_q = 0.5*pow(params.T/5, 0.45)*pow(0.05/rho, 0.5); // fm
  params.delta_p = 0.4*sqrt(parts.mass*params.T); // MeV/c
  char filename_config[255];
  sprintf(filename_config, "data_sin_pauli/config_%.3f_%.3f.lammpstrj", rho, params.T);
  t = time(NULL);
  long int aceptados = 0;
  for(int i = 0; i < 100; i++){
    printf("\rProgreso: %d%%", i);
    fflush(stdout);
    aceptados += N_steps(&parts, &pot_tot, &params, parts.n*factor_term);
  }
  printf("%ld%% en %d segundos\n", (aceptados)/(parts.n*factor_term), (int)time(NULL)-t);
  save_lammpstrj(filename_config, &parts, params.L, 0);
*/
  liberar(&parts);
  liberar_LUTs(&pot_tot);
  return 0;
}
