#include "general.h"
#include "celda.h"
#include "interaccion.h"
#include "tabla.h"
#include "avanzar.h"
#include "inicializar.h"
#include "energia.h"
#include "muestreo.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <string.h>


int main(int argc, char *argv[]){

  srand(1);
  float rho = 0.05;
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
  pauli.scut2 = 10; // 0.7% del maximo
  pauli.scut2 = 5.4*5.4/(pauli.qo*pauli.qo); // 0.7% del maximo
  pauli.shift = pauli.D*exp(-0.5*pauli.scut2);

  struct Panda_np panda_np;
  panda_np.mu_r = 1.7468 - (-3.613281250000000000e-01); // fm^-1
  panda_np.mu_a = 1.6 - (-3.609375000000000222e-01); // fm^-1
  panda_np.V_r = 3088.118*2.453125000000000000e+00; // MeV
  panda_np.V_a = 2666.647*2.621093750000000000e+00; // MeV
  panda_np.rcut = 5.4; // fm
  panda_np.rcut2 = 5.4*5.4; // fm^2
  panda_np.shift = 0;
  panda_np.shift = interaction_panda_np(panda_np.rcut, &panda_np);

  struct Panda_nn panda_nn;
  panda_nn.mu_o = 1.5 - (-3.484375000000000111e-01); // fm^-1
  panda_nn.V_o = 373.118*3.195312500000000000e+00; // MeV
  panda_nn.rcut = 5.4; // fm
  panda_nn.rcut2 = 5.4*5.4; // fm^2
  panda_nn.shift = 0;
  panda_nn.shift = interaction_panda_nn(panda_nn.rcut, &panda_nn);

  struct QCNM qcnm;
  qcnm.V_o = 25.93*1.496093750000000000e+00; //MeV
  qcnm.p_1 = 6.2 + 1.679687500000000000e-01;
  qcnm.p_2 = 3.0 + 6.542968750000000000e-01;
  qcnm.r_1 = 1.757 + 1.718749999999999445e-02; //fm
  qcnm.r_2 = 1.771 - 4.296875000000000000e-02; //fm
  qcnm.d = 3.35 - 1.066406249999999889e-01; //fm
  qcnm.a = 5.0/6.0 + 4.375000000000001110e-02; //fm
  qcnm.rcut = 5.4; //fm
  qcnm.rcut2 = 5.4*5.4; //fm*fm
  qcnm.shift = interaction_QCNM(qcnm.rcut, &qcnm);

  struct Coulomb coul;
  coul.q2 = 1.4403427984368629; // Mev*fm
  coul.lambda = 20; // fm
  coul.rcut = coul.lambda; // fm
  coul.rcut2 = coul.rcut*coul.rcut; // fm*fm
  coul.shift = 0;
  coul.shift = interaction_Coulomb(coul.rcut, &coul);

  struct Trap trap;
  trap.k = 0;
  trap.power = 6;

  struct Interaction pot_tot;
  pot_tot.pauli = &pauli;
  pot_tot.panda_nn = &panda_nn;
  pot_tot.panda_np = &panda_np;
  pot_tot.qcnm = &qcnm;
  pot_tot.coul = &coul;
  pot_tot.trap = &trap;
  pot_tot.rcut = 5.4; // fm
  reset_energies(&pot_tot);
  int N_LUT = 100000;
  panda_np.rcut = 0;
  panda_nn.rcut = 0;
  panda_np.V_r = 0;
  panda_np.V_a = 0;
  panda_nn.V_o = 0;
  //pauli.D = 0;
  build_LUTs(&pot_tot, 0, 0, 0, N_LUT, N_LUT, 0);
  //build_LUTs(&pot_tot, 0, 0, 0, N_LUT, 0, 0);
  char filename[255], nombre[4];
  sprintf(filename, "data/test/tabla_QCNM.txt");
  FILE* fp = fopen(filename, "w");
  for(int i = 0; i < N_LUT; i++){
    fprintf(fp, "%f %f\n", (i+1)*qcnm.dr2, qcnm.LUT[i]);
  }
  fclose(fp);

// Parametros termodinamicos
  struct Externos params;
  int N = 18; //16
  params.L = N/pow(rho, 1.0/3); // mayor a 2*qo*scut, rho en fm^-3
  N = N*N*N;
  params.T = 4; // MeV

  if(opcion == 't'){
    int comps[4] = {N/4, N/4, N/4, N/4};
    srand(1);

    inicializar(&parts, comps, 4, 938, params.L, params.T, &pot_tot, 0);

    float dT = 0.1;
    int term = 10000, desc = 25, samp = 100;
    char filename_termo[255], filename_config[255];
    sprintf(filename_termo, "data/test/termo_%.3f.txt", rho);
    FILE* fp = fopen(filename_termo, "w");
    fprintf(fp, "#Temper   Kinetic  Pandha_np Pandha_nn QCNM   Coulomb   Pauli Pressure Temp_Vi Time Acept\n");
    fclose(fp);
    for(int t = 0; t < 40; t++){
      params.delta_q = 0.4*sqrt(params.T/4)*sqrt(0.05/rho); // fm
      params.delta_p = 0.4*sqrt(parts.mass*params.T); // MeV/c
      ajustar_deltas_full(&parts, &pot_tot, &params, 0.3, 0.6);
      printf("rho=%.2f fm^-3 | T = %.1fMeV\n", rho, params.T);
      muestrear_termo(filename_termo, &parts, &pot_tot, &params, samp, desc, term);
      if(t % 5 == 0 || t > 34){
        sprintf(filename_config, "data/test/config_%.3f_%.3f.lammpstrj", rho, params.T);
        save_lammpstrj(filename_config, &parts, params.L, 0);
      }
      params.T -= dT;
    }

  }

  if(opcion == 'n'){
    int Np = 20, Nn = 20;
    sprintf(nombre, "Ca");
    float R;
    if (argc >= 3){
      sscanf(argv[2], "%s\n", nombre);
      if(argc >= 4){
        sscanf(argv[3], "%d\n", &Np);
        if(argc >= 5){
          sscanf(argv[4], "%d\n", &Nn);
          if(argc>=6){
            sscanf(argv[5], "%f\n", &R);
          }else{
            R = 2*cbrt(Nn+Np);
          }
        }
      }else{
        Nn = Np;
        R = 2*cbrt(Nn+Np);
      }
    }
    params.L = 150; // fm
    params.T = 1; // MeV
    liberar_LUTs(&pot_tot);
    trap.k = 1;
    trap.k = 1/interaction_trap(2*R, &trap);
    trap.rcut2 = params.L*params.L/4;
    coul.rcut = params.L/3;    coul.lambda = 1000; // fm, sin apantallar
    panda_np.rcut = params.L/3;    panda_np.rcut2 = panda_np.rcut*panda_np.rcut;
    panda_nn.rcut = params.L/3;    panda_nn.rcut2 = panda_nn.rcut*panda_nn.rcut;
    pauli.scut2 = params.L/(3*pauli.qo*pauli.qo);
    qcnm.rcut = params.L/3;
    pot_tot.rcut = params.L/3;
    build_LUTs(&pot_tot, 0, 0, 0, N_LUT, N_LUT, N_LUT);

    int comps[4] = {(Np+1)/2, Np/2, (Nn+1)/2, Nn/2};
    //int comps[4] = {Np, 0, Nn, 0};
    inicializar_nucleo(&parts, comps, 4, 938, params.L, params.T, R, &pot_tot, 0);

    int term = 150, desc = 5, samp = 10, N_temps = 100, trap_duration = N_temps/2;
    float dT = 1.0/N_temps, dk = trap.k/trap_duration;
    N_steps(&parts, &pot_tot, &params, term*parts.n);
    char filename_termo[255], filename_config[255];
    sprintf(filename_termo, "data/nucleos/proceso_%s_%d_%d.txt", nombre, Np, Nn);
    FILE* fp = fopen(filename_termo, "w");
    fprintf(fp, "#Temper   Kinetic  Pandha_np Pandha_nn QCNM    Coulomb   Pauli   Radius   Temp_Vi   Time  Accept\n");
    fclose(fp);
    for(int t = 0; t < 100; t++){
      params.delta_q = 0.2*sqrt(params.T)*(R*cbrt(40.0/parts.n)/7); // fm
      params.delta_p = 0.2*sqrt(parts.mass*params.T)*(R*cbrt(40.0/parts.n)/7); // MeV/c
      ajustar_deltas_full(&parts, &pot_tot, &params, 0.3, 0.6);
      muestrear_nucleo(filename_termo, &parts, &pot_tot, &params, samp, desc, term);
      if(trap.k < dk){
        trap.k = 0;
      }else{
        trap.k -= dk;
      }
      params.T -= dT;
    }
    sprintf(filename_config, "data/nucleos/config_nuc_%s_%d_%d.lammpstrj", nombre, Np, Nn);
    save_lammpstrj(filename_config, &parts, params.L, 0);
  }

  liberar(&parts);
  liberar_LUTs(&pot_tot);
  return 0;
}
