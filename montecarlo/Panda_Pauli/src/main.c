#include "general.h"
#include "celda.h"
#include "interaccion.h"
#include "avanzar.h"
#include "inicializar.h"
#include "muestreo.h"
#include "virial.h"
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
  pauli.scut2 = 10; // 0.7% del maximo
  pauli.scut2 = 5.4*5.4/(pauli.qo*pauli.qo); // 0.7% del maximo
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
    liberar_LUTs(&pot_tot);
    pauli.D = 207; // MeV
    panda_np.V_r = 2.5*3088.118; // MeV
    panda_np.V_a = 2.5*2666.647; // MeV
    panda_nn.V_o = 2*373.118; // MeV
    build_LUTs(&pot_tot, N_LUT, N_LUT, N_LUT);
    // Particulas
    float mass = 938; // MeV/c^2
    N = 512;
    params.L = pow(N/rho, 1.0/3); // mayor a 2*qo*scut, rho en fm^-3
    int comps[4] = {N/4, N/4, N/4, N/4};
    inicializar(&parts, comps, 4, mass, params.L, params.T, &pot_tot);
    energia(&parts, &pot_tot);
    params.delta_q = 0.5*sqrt(params.T/4)*sqrt(0.05/rho); // fm
    params.delta_p = 0.5*sqrt(parts.mass*params.T); // MeV/c
    char filename_termo[255], filename_config[255];
    int factor_term = 4000, factor_desc = 10, Nsamp = 100;
    int append = 0; // Solo un frame para lammps
    N_steps(&parts, &pot_tot, &params, parts.n*factor_term);
    float dT = 0.1;
    sprintf(filename_termo, "data_con_pauli/termo_2.5_2.5_2_%.3f.txt", rho);

    FILE* fp = fopen(filename_termo, "w");
    fclose(fp);
    for (int t = 0; t < 40; t++){
      energia(&parts, &pot_tot);
      params.delta_q = 0.5*sqrt(params.T/4)*sqrt(0.05/rho); // fm
      params.delta_p = 0.5*sqrt(parts.mass*params.T); // MeV/c
      printf("rho = %.2f fm^-3 | T = %.2f MeV\n", rho, params.T);
      muestrear_termo(filename_termo, &parts, &pot_tot, &params, Nsamp, factor_desc, factor_term);
      sprintf(filename_config, "data_con_pauli/config_2.5_2.5_2_%.3f_%.3f.lammpstrj", rho, params.T);
      if (t%10==0 || t>34) save_lammpstrj(filename_config, &parts, params.L, append);
      params.T -= dT;
    }
  }


// MUESTREO EN TEMPERATURA - 1000 particulas
  if(opcion=='n'){
    // Particulas
    N = 12;
    params.L = N/pow(rho, 1.0/3); // mayor a 2*qo*scut, rho en fm^-3
    N = N*N*N;
    float mass = 938; // MeV/c^2
    int comps[4] = {N/4, N/4, N/4, N/4};
    inicializar(&parts, comps, 4, mass, params.L, params.T, &pot_tot);
    energia(&parts, &pot_tot);
    params.delta_q = 0.5*pow(0.05/rho, 0.5); // fm
    params.delta_p = 0.5*sqrt(parts.mass*params.T); // MeV/c
    char filename_termo[255], filename_config[255];
    int factor_term = 10000, factor_desc = 25, Nsamp = 100;
    int append = 0; // Solo un frame para lammps
    N_steps(&parts, &pot_tot, &params, parts.n*factor_term);
    float dT = 0.1;
    sprintf(filename_termo, "data_con_pauli_1728/termo_%.3f.txt", rho);

    FILE* fp = fopen(filename_termo, "w");
    fclose(fp);
    for (int t = 0; t < 40; t++){
      energia(&parts, &pot_tot);
      params.delta_q = 0.5*sqrt(params.T/4)*pow(0.05/rho, 0.5); // fm
      params.delta_p = 0.5*sqrt(parts.mass*params.T); // MeV/c
      printf("rho = %.2f fm^-3 | T = %.2f MeV\n", rho, params.T);
      muestrear_termo(filename_termo, &parts, &pot_tot, &params, Nsamp, factor_desc, factor_term);
      sprintf(filename_config, "data_con_pauli_1728/config_%.3f_%.3f.lammpstrj", rho, params.T);
      if (t%5==0 || t>34) save_lammpstrj(filename_config, &parts, params.L, append);
      params.T -= dT;
    }
  }


  if(opcion=='s'){
    printf("Opcion para hallar fundamental en momentos para una SC fija\n\n");
    // Particulas
    int segs = time(NULL);
    char filename_config[255];
    sprintf(filename_config, "fundamental/config_redondear_%.3f_0.100.lammpstrj", rho);
    load_and_rep(filename_config, &parts, pot_tot.rcut, &params.L, 3);
    params.T = 0.1;
    parts.mass = 938;
    energia(&parts, &pot_tot);
    params.delta_q = 0; // fm
    params.delta_p = 0.8*sqrt(parts.mass*params.T); // MeV/c
    float E_tot_init = (parts.kinetic + parts.energy_panda + parts.energy_pauli)/parts.n, PyT[2], Teff, E_tot_final;
    printf("Inicial: %f + %f + %f = %f\n", parts.kinetic/parts.n, parts.energy_panda/parts.n, parts.energy_pauli/parts.n, E_tot_init);

    int factor_term = 5000, aceptados = 0, N_temps = 10;
    for(int t = 0; t < N_temps; t++){
      Teff = 0;
      aceptados = 0;
      params.delta_p = 0.8*sqrt(parts.mass*params.T); // MeV/c
      for(int i = 0; i < 100; i++){
        printf("Progress: %d%%\r", i);
        fflush(stdout);
        aceptados += N_steps(&parts, &pot_tot, &params, parts.n*factor_term/100);
        if (70 <= i){
          presion_temperatura(&parts, &pot_tot, &params, PyT);
          Teff += PyT[1]/30;
        }
      }
      E_tot_final = (parts.kinetic + parts.energy_panda + parts.energy_pauli)/parts.n;
      printf("T=%.3f: %f + %f + %f = %f\n", Teff, parts.kinetic/parts.n, parts.energy_panda/parts.n, parts.energy_pauli/parts.n, E_tot_final);
      printf("Aceptacion: %d%%\n", (100*aceptados)/(factor_term*parts.n));
      params.T -= 0.01;
    }
    sprintf(filename_config, "fundamental/config_1728_x1_%.3f_0.010.lammpstrj", rho);
    save_lammpstrj(filename_config, &parts, params.L, 0);
    /*
    E_tot_final = (parts.kinetic + parts.energy_panda + parts.energy_pauli)/parts.n;
    printf("  Final: %f) %f + %f + %f = %f | T = %.3f\n", parts.kinetic/parts.n, parts.energy_panda/parts.n, parts.energy_pauli/parts.n, E_tot_final, Teff);
    if(E_tot_init > E_tot_final)  printf("La energia total bajo\n");
    else  printf("La energia total subio\n\n");
    */
    segs = time(NULL) - segs;
    printf("Duracion: %d hs, %d min, %d segs\n", segs/(60*60), (segs/60)%60, segs%60);
  }
  if(opcion=='f'){
    float a = 1;
    if (argc >= 4){
      sscanf(argv[3], "%f\n", &a);
    }
    printf("Opcion para hallar fundamental con periodicidad SC con Vnp x%f\n", a);
    // Particulas
    float mass = 938; // MeV/c^2
    N = 8*8*8;
    N = 4*4*4;
    params.L = pow(N/rho, 1.0/3);
    int comps[4] = {N/4, N/4, N/4, N/4};
    inicializar(&parts, comps, 4, mass, params.L, params.T, &pot_tot);
    energia(&parts, &pot_tot);
    params.delta_q = 0.3*pow(params.T/4, 0.8)*sqrt(0.05/rho); // fm
    params.delta_p = 0.3*sqrt(4*parts.mass)*pow(params.T/4, 0.8); // MeV/c

    pot_tot.pauli->D = 207; // MeV
    panda_np.V_r = a*3088.118; // MeV
    panda_np.V_a = a*2666.647; // MeV
    liberar_LUTs(&pot_tot);
    build_LUTs(&pot_tot, N_LUT, N_LUT, N_LUT);
    energia(&parts, &pot_tot);

    char filename_termo[255], filename_config[255];
    int factor_term = 45000, factor_desc = 50, Nsamp = 100;
    int append = 0; // Solo un frame para lammps
    N_steps(&parts, &pot_tot, &params, parts.n*factor_term);
    float dT = 0.1;
    sprintf(filename_termo, "fundamental/termo_512_x%f_%.3f.txt", a, rho);
    sprintf(filename_termo, "fundamental/termo_fact_x1_%.3f.txt", rho);
    sprintf(filename_termo, "fundamental/termo_x1_%.3f.txt", rho);
    FILE* fp = fopen(filename_termo, "w");
    fclose(fp);
    for (int t = 0; t < 40; t++){
      energia(&parts, &pot_tot);
      params.delta_q = 0.3*pow(params.T/4, 0.8)*sqrt(0.05/rho); // fm
      params.delta_p = 0.25*sqrt(4*parts.mass)*pow(params.T/4, 0.8); // MeV/c
      printf("rho = %.2f fm^-3 | T = %.2f MeV\n", rho, params.T);
      muestrear_termo(filename_termo, &parts, &pot_tot, &params, Nsamp, factor_desc, factor_term);
      sprintf(filename_config, "fundamental/config_512_x%f_%.3f_%.3f.lammpstrj", a, rho, params.T);
      sprintf(filename_config, "fundamental/config_fact_x1_%.3f_%.3f.lammpstrj", rho, params.T);
      sprintf(filename_config, "fundamental/config_x1_%.3f_%.3f.lammpstrj", rho, params.T);
      if (t%10==0 || t>34) save_lammpstrj(filename_config, &parts, params.L, append);
      params.T -= dT;
    }
  }
  if(opcion=='l'){
    float a = 1.562, b = 1;
    if (argc >= 4){
      sscanf(argv[3], "%f\n", &a);
      if (argc >= 5){
        sscanf(argv[4], "%f\n", &b);
      }
    }
    char filename[255], filename_data[255];
    sprintf(filename_data, "fundamental/Evsrho_temp.txt");
    FILE *fp = fopen(filename_data, "w");
    float rhos[5] = {0.15, 0.155, 0.16, 0.165, 0.17};
    float Es[5];
    pot_tot.pauli->D = 207; // MeV
    panda_np.V_r = a*3088.118; // MeV
    panda_np.V_a = a*2666.647; // MeV
    panda_nn.V_o = b*373.118; // MeV
    liberar_LUTs(&pot_tot);
    build_LUTs(&pot_tot, N_LUT, N_LUT, N_LUT);
    fprintf(fp, "#  rho    E     Ekin      Enuc\n");
    printf("%f\n", pauli.shift);
    printf("#  rho    E     Ekin      Enuc\n");
    for(int r = 0; r < 5; r++){
      sprintf(filename, "fundamental/config_1728_x1_%.3f_0.010.lammpstrj", rhos[r]);
      //load_and_rep(filename, &parts, pot_tot.rcut, &params.L, 5);
      load_lammpstrj(filename, &parts, &params.L, pot_tot.rcut);
      parts.mass = 938;
      energia(&parts, &pot_tot);
      Es[r] = (parts.kinetic+parts.energy_panda+parts.energy_pauli)/parts.n;
      printf("%f %f %f %f %f\n", rhos[r], Es[r], parts.kinetic/parts.n, parts.energy_panda/parts.n, parts.energy_pauli/parts.n);
      energia_sin_LUT(&parts, &pot_tot);
      Es[r] = (parts.kinetic+parts.energy_panda+parts.energy_pauli)/parts.n;
      //Es[r] = parts.kinetic/parts.n;
      fprintf(fp, "%f %f %f %f %f\n", rhos[r], Es[r], parts.kinetic/parts.n, parts.energy_panda/parts.n, parts.energy_pauli/parts.n);
      printf("%f %f %f %f %f\n", rhos[r], Es[r], parts.kinetic/parts.n, parts.energy_panda/parts.n, parts.energy_pauli/parts.n);
    }
    fclose(fp);
  }
  if(opcion=='c'){
    printf("Opcion para confirmar fundamental con periodicidad SC\n");
    char filename[255];
    sprintf(filename, "fundamental/config_x1_%.3f_0.100.lammpstrj", rho);
    load_and_rep(filename, &parts, pot_tot.rcut, &params.L, 3);
    parts.mass = 938; // MeV/c^2
    pot_tot.pauli->D = 207; // MeV
    panda_np.V_r = 3088.118; // MeV
    panda_np.V_a = 2666.647; // MeV
    liberar_LUTs(&pot_tot);
    build_LUTs(&pot_tot, N_LUT, N_LUT, N_LUT);
    energia(&parts, &pot_tot);

    char filename_termo[255], filename_config[255];
    sprintf(filename_termo, "fundamental/termo_recalc_x1_%.3f.txt", rho);
    sprintf(filename_config, "fundamental/config_recalc_x1_%.3f.lammpstrj", rho);
    int factor_term = 20000, factor_desc = 50, Nsamp = 100;
    N_steps(&parts, &pot_tot, &params, parts.n*factor_term);
    params.T = 0.1;
    float dT = 0.1;
    FILE* fp = fopen(filename_termo, "w");
    fclose(fp);
    for (int t = 0; t < 9; t++){
      energia(&parts, &pot_tot);
      params.delta_q = 0.3*pow(params.T/4, 0.8)*sqrt(0.05/rho); // fm
      params.delta_p = 0.25*sqrt(4*parts.mass)*pow(params.T/4, 0.8); // MeV/c
      printf("rho = %.2f fm^-3 | T = %.2f MeV\n", rho, params.T);
      muestrear_termo(filename_termo, &parts, &pot_tot, &params, Nsamp, factor_desc, factor_term);
      save_lammpstrj(filename_config, &parts, params.L, t);
      params.T += (2*(t<4) - 1)*dT;
    }
  }
  if(opcion=='r'){
    printf("Opcion para redondear red del fundamental con periodicidad SC\n");
    parts.mass = 938;
    pot_tot.pauli->D = 207; // MeV
    panda_np.V_r = 3088.118; // MeV
    panda_np.V_a = 2666.647; // MeV
    liberar_LUTs(&pot_tot);
    build_LUTs(&pot_tot, N_LUT, N_LUT, N_LUT);

    char filename[255];
    sprintf(filename, "fundamental/config_redondear_%.3f_0.100.lammpstrj", rho);
    sprintf(filename, "fundamental/config_x1_%.3f_0.100.lammpstrj", rho);
    load_lammpstrj(filename, &parts, &params.L, pot_tot.rcut);
    energia(&parts, &pot_tot);
    printf("%.3f)  Inicial  %f + %f + %f = %f\n", rho, parts.kinetic/parts.n, parts.energy_panda/parts.n, parts.energy_pauli/parts.n, (parts.kinetic+parts.energy_panda+parts.energy_pauli)/parts.n);

    redondear_SC(&parts, pot_tot.rcut, params.L);
    energia(&parts, &pot_tot);

    params.T = 0.1;
    params.delta_q = 0;
    params.delta_p = 2.00*sqrt(4*parts.mass)*pow(params.T/4, 0.8); // MeV/c
    int factor_term = 1000000;
    factor_term = 10000;
    long int aceptados = N_steps(&parts, &pot_tot, &params, parts.n*factor_term);
    energia(&parts, &pot_tot);

    printf("%.3f) Final (%ld%%) %f + %f + %f = %f\n", rho, 100*aceptados/(parts.n*factor_term), parts.kinetic/parts.n, parts.energy_panda/parts.n, parts.energy_pauli/parts.n, (parts.kinetic+parts.energy_panda+parts.energy_pauli)/parts.n);
    sprintf(filename, "fundamental/config_redondear_512_%.3f_0.100.lammpstrj", rho);
    save_lammpstrj(filename, &parts, params.L, 0);
  }

  if(opcion=='d'){
    printf("Opcion para calcular saturaciones en SC via dilatacion/contraccion\n");
    char filename[255];
    sprintf(filename, "fundamental/config_redondear_%.3f_0.100.lammpstrj", rho);
    load_and_rep(filename, &parts, pot_tot.rcut, &params.L, 5);
    parts.mass = 938;
    pot_tot.pauli->D = 207; // MeV
    panda_np.V_r = 3088.118; // MeV
    panda_np.V_a = 2666.647; // MeV
    liberar_LUTs(&pot_tot);
    build_LUTs(&pot_tot, N_LUT, N_LUT, N_LUT);
    energia(&parts, &pot_tot);
    printf("%f + %f + %f = %f\n", parts.kinetic/parts.n, parts.energy_panda/parts.n, parts.energy_pauli/parts.n, (parts.kinetic+parts.energy_panda+parts.energy_pauli)/parts.n);

    float f = cbrt(0.15/rho);
    printf("x%.3f:  ", f);
    for(int i = 0; i < 3*parts.n; i++){
      parts.q[i] *= f;
    }
    params.L *= f;
    parts.l *= f;
    energia(&parts, &pot_tot);
    printf("%f + %f + %f = %f\n", parts.kinetic/parts.n, parts.energy_panda/parts.n, parts.energy_pauli/parts.n, (parts.kinetic+parts.energy_panda+parts.energy_pauli)/parts.n);

    f = cbrt(0.17/rho);
    printf("x%.3f:  ", f);
    for(int i = 0; i < 3*parts.n; i++){
      parts.q[i] *= f;
    }
    params.L *= f;
    parts.l *= f;
    energia(&parts, &pot_tot);
    printf("%f + %f + %f = %f\n", parts.kinetic/parts.n, parts.energy_panda/parts.n, parts.energy_pauli/parts.n, (parts.kinetic+parts.energy_panda+parts.energy_pauli)/parts.n);
  }
  if(opcion=='g'){
    printf("Calcular la energia de una configuracion de Guille\n");
    char filename[255];
    float PyT[2];
    sprintf(filename, "fundamental/config_redondear_%.3f_0.100.lammpstrj", rho);
    //sprintf(filename, "pandha_pauli_spin_mc_n64_rho18_delta001.lammpstrj");
    //sprintf(filename, "pandha_pauli_mc.lammpstrj");
    //load_lammpstrj_guille(filename, &parts, &params.L, pot_tot.rcut);
    load_lammpstrj(filename, &parts, &params.L, pot_tot.rcut);
    //sprintf(filename, "prueba.lammpstrj");
    //save_lammpstrj(filename, &parts, params.L, 0);
    parts.mass = 938;
    pot_tot.pauli->D = 207; // MeV
    panda_np.V_r = 3088.118; // MeV
    panda_np.V_a = 2666.647; // MeV
    liberar_LUTs(&pot_tot);
    build_LUTs(&pot_tot, N_LUT, N_LUT, N_LUT);
    energia(&parts, &pot_tot);
    presion_temperatura(&parts, &pot_tot, &params, PyT);
    printf("%f\n", PyT[1]);
    printf("%f + %f + %f = %f\n", parts.kinetic/parts.n, parts.energy_panda/parts.n, parts.energy_pauli/parts.n, (parts.kinetic+parts.energy_panda+parts.energy_pauli)/parts.n);
  }

  if(opcion=='2'){
    printf("Calcular la energia de una configuracion de 2 particulas\n");
    char filename[255];
    sprintf(filename, "pandha_pauli_initial.lammpstrj");
    load_lammpstrj_guille(filename, &parts, &params.L, pot_tot.rcut);
    int t0, t1;
    if (argc >= 3){
      sscanf(argv[2], "%d\n", &t0);
      parts.type[0] = t0;
      if (argc >= 4){
        sscanf(argv[3], "%d\n", &t1);
        parts.type[1] = t1;
      }
    }
    parts.mass = 938;
    pot_tot.pauli->D = 207; // MeV
    panda_np.V_r = 3088.118; // MeV
    panda_np.V_a = 2666.647; // MeV
    liberar_LUTs(&pot_tot);
    build_LUTs(&pot_tot, N_LUT, N_LUT, N_LUT);
    energia(&parts, &pot_tot);
    int delta_idx[3] = {0, 0, 0};
    float rsq = distancia(parts.q, parts.q+3, delta_idx, 0.0);
    float s2 = rsq/(pauli.qo*pauli.qo) + distancia_p(parts.p, parts.p+3)/(pauli.po*pauli.po);
    printf("c0 = %d | c1 = %d  |||  dq = %f  |||  dp = %f |||  ds = %f\n", parts.celda[0], parts.celda[1], sqrt(rsq), sqrt(distancia_p(parts.p, parts.p+3)), sqrt(s2));
    printf("%f + %f + %f = %f\n", parts.kinetic/parts.n, parts.energy_panda/parts.n, parts.energy_pauli/parts.n, (parts.kinetic+parts.energy_panda+parts.energy_pauli)/parts.n);
    float PyT[2];
    presion_temperatura(&parts, &pot_tot, &params, PyT);
    printf("%f %f\n", PyT[0], PyT[1]);
  }

  if(opcion=='p'){
    printf("Opcion para encontrar fundamental en p via punto fijo\n");
    parts.mass = 938;
    pot_tot.pauli->D = 207; // MeV
    panda_np.V_r = 3088.118; // MeV
    panda_np.V_a = 2666.647; // MeV
    liberar_LUTs(&pot_tot);
    build_LUTs(&pot_tot, N_LUT, N_LUT, N_LUT);

    char filename[255];
    sprintf(filename, "fundamental/config_redondear_%.3f_0.100.lammpstrj", rho);
    load_lammpstrj(filename, &parts, &params.L, pot_tot.rcut);
    energia(&parts, &pot_tot);
    printf("%f ", (parts.kinetic+parts.energy_panda+parts.energy_pauli)/parts.n);

    sprintf(filename, "fundamental/punto_fijo_%.3f.txt", rho);
    FILE *fp = fopen(filename, "w");
    float *gorces = (float *) malloc(3*parts.n*sizeof(float));
    int N_iter = 100000, idx;
    printf("%f\n", parts.mass);
    fprintf(fp, "%f ", (parts.kinetic+parts.energy_panda+parts.energy_pauli)/parts.n);
    for(int i = 0; i < N_iter; i++){
      punto_fijo(&parts, &pot_tot, &params, gorces);
      for(int l = 0; l < 3*parts.n; l++)  parts.p[l] = gorces[l];
      idx = min_vec(parts.p, parts.n);
      if(i%1000==0)  printf("%f  ", sqrt(norma(parts.p+3*idx)));
      idx = max_vec(parts.p, parts.n);
      if(i%1000==0)  printf("%f\n", sqrt(norma(parts.p+3*idx)));
      energia(&parts, &pot_tot);
      fprintf(fp, "%f ", (parts.kinetic+parts.energy_panda+parts.energy_pauli)/parts.n);
    }
    fclose(fp);
    free(gorces);
    float *aux = parts.q;
    parts.q = parts.p;
    parts.p = aux;
    sprintf(filename, "prueba.lammpstrj");
    save_lammpstrj(filename, &parts, params.L, 0);
  }
  liberar(&parts);
  liberar_LUTs(&pot_tot);
  return 0;
}
