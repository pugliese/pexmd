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

  srand(time(NULL));
  float rho = 0.2;
  int P = 1;
  if (argc >= 2){
    sscanf(argv[1], "%f\n", &rho);
    if (argc >= 3){
      sscanf(argv[2], "%d\n", &P);
    }
  }

// Potenciales
  struct Pauli pauli;
  pauli.qo = 1.644; // fm
  pauli.po = 120; // MeV/c
  pauli.D = 207; // MeV
  pauli.D = 0; // Por ahora sin Pauli
  pauli.scut2 = 10; // 0.7% del maximo
  pauli.shift = pauli.D*exp(-0.5*pauli.scut2) - pauli.shift;

  struct Panda_np panda_np;
  panda_np.mu_r = 1.7468; // fm^-1
  panda_np.mu_a = 1.6; // fm^-1
  panda_np.V_r = 3088.118; // MeV
  panda_np.V_a = 2666.647; // MeV
  panda_np.rcut = 5.4; // fm
  panda_np.shift = 0;
  panda_np.shift = interaction_panda_np(panda_np.rcut, &panda_np);

  struct Panda_nn panda_nn;
  panda_nn.mu_o = 1.5; // fm^-1
  panda_nn.V_o = 373.118; // MeV
  panda_nn.rcut = 5.4; // fm
  panda_nn.shift = 0;
  panda_nn.shift = interaction_panda_nn(panda_nn.rcut, &panda_nn);

  struct Interaction pot_tot;
  pot_tot.pauli = &pauli;
  pot_tot.panda_nn = &panda_nn;
  pot_tot.panda_np = &panda_np;
	pot_tot.rcut = max(sqrt(pauli.scut2)*pauli.qo, panda_nn.rcut);

// Parametros termodinamicos
  struct Externos params;
  int N = 20;
  params.L = N/pow(rho, 1.0/3); // mayor a 2*qo*scut, rho en fm^-3
  N = N*N*N;
  params.T = 5; // MeV

// Particulas
  struct Particles parts;
  float mass = 938; // MeV/c^2
  int comps[4] = {N/4, N/4, N/4, N/4};
  inicializar(&parts, comps, 4, mass, params.L, params.T, &pot_tot);
  energia(&parts, &pot_tot);

  params.delta_q = 0.7*pow(0.05/rho, 0.5); // fm
  params.delta_p = 0.5*sqrt(parts.mass*params.T); // MeV/c
/*
  int factor_term = 1000, factor_desc = 10, Nsamp = 100;
  save_lammpstrj("Inicial_random.lammpstrj", &parts, params.L, 0);
  N_steps(&parts, &pot_tot, &params, parts.n*factor_term);
  save_lammpstrj("Tras_1000_pasos_random.lammpstrj", &parts, params.L, 0);
  N_steps(&parts, &pot_tot, &params, parts.n*factor_term);
  save_lammpstrj("Tras_2000_pasos_random.lammpstrj", &parts, params.L, 0);
  N_steps(&parts, &pot_tot, &params, parts.n*factor_term);
  save_lammpstrj("Tras_3000_pasos_random.lammpstrj", &parts, params.L, 0);
  N_steps(&parts, &pot_tot, &params, parts.n*factor_term);
  save_lammpstrj("Tras_4000_pasos_random.lammpstrj", &parts, params.L, 0);
  N_steps(&parts, &pot_tot, &params, parts.n*factor_term);
  save_lammpstrj("Tras_5000_pasos_random.lammpstrj", &parts, params.L, 0);
*/
// MUESTREO
  char filename_termo[255], filename_config[255];
  int factor_term = 2000, factor_desc = 10, Nsamp = 100;
  int append = 0; // Solo un frame para lammps
  N_steps(&parts, &pot_tot, &params, parts.n*factor_term);
  float dT = 0.1;
  sprintf(filename_termo, "data_sin_pauli/termo_%.3f.txt", rho);

  FILE* fp = fopen(filename_termo, "w");
  fclose(fp);
  for (int t = 0; t < 50; t++){
    energia(&parts, &pot_tot);
    params.delta_q = 0.7*pow(params.T/5, 0.45)*pow(0.05/rho, 0.5); // fm
    params.delta_p = 0.5*sqrt(parts.mass*params.T); // MeV/c
    printf("rho = %.2f fm^-3 | T = %.2f MeV\n", rho, params.T);
    muestrear_termo(filename_termo, &parts, &pot_tot, &params, Nsamp, factor_desc, factor_term);
    sprintf(filename_config, "data_sin_pauli/config_%.3f_%.3f.lammpstrj", rho, params.T);
    save_lammpstrj(filename_config, &parts, params.L, append);
    params.T -= dT;
  }

  liberar(&parts);
  return 0;
}
