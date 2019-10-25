import numpy as np
import matplotlib.pylab as plt
import matplotlib.animation as animation
import matplotlib.gridspec as gridspec
import sys
from GA import *
from mediador import *
plt.ion()

red = "B1"
option = sys.argv[1]
option_split = option.split("-")
option = option_split[0]
if(len(option_split) > 1):
  red = option_split[1]

D_estrella = 207*938/(120**2)
otro_po = 120
formato1 = "Pauli_puro/B1/config_maruyama_dqx0.0_%.3f_0.001.lammpstrj"
formato2 = "Pauli_puro/B2/config_B2_fija_%.3f_0.001.lammpstrj"
formato3 = "Pauli_puro/B3/config_B3_fija_%.3f_0.001.lammpstrj"
formato = formato1*(red=="B1") + formato2*(red=="B2") + formato3*(red=="B3")

pauli = Pauli(qo = 1.644, po = otro_po, D = D_estrella*(otro_po**2)/938, scut2 = (5.4/1.644)**2)
panda_nn = Panda_nn(mu_o = 1.5, V_o = 373.118, rcut = 5.4)
panda_np = Panda_np(mu_r = 1.7468, mu_a = 1.6, V_r = 3088.118, V_a = 2666.647, rcut = 5.4)
qcnm = QCNM(V_o=25.93, p_1=6.2, p_2=3.0, r_1=1.757, r_2=1.771, d=3.35, a=5./6, rcut=5.4)
pot_tot = TotalPotential(pauli, panda_nn, panda_np, qcnm)

rhos = np.array([0.150, 0.153, 0.155, 0.158, 0.160, 0.162, 0.165, 0.167, 0.170])
rhos = np.array([0.150,        0.155, 0.158, 0.160, 0.162, 0.165,        0.170])
rhos = np.array([0.150,        0.155, 0.158, 0.160, 0.162, 0.165, 0.167, 0.170])
E_NM_pandha = np.zeros_like(rhos)
E_NM_pauli = np.zeros_like(rhos)
E_NM_kinetic = np.zeros_like(rhos)
vec_parts = [0 for rho in rhos]
L = np.array([3.0], dtype = np.float32)
for i in range(len(rhos)):
  vec_parts[i] = Particles(1, 938)
  vec_parts[i].load_lammpstrj(formato1 %rhos[i], L, 5.4)
  vec_parts[i].energy(pot_tot)
  E_NM_pandha[i] = vec_parts[i].energy_panda/vec_parts[i].n
  vec_parts[i] = Particles(1, 938)
  vec_parts[i].load_lammpstrj(formato %rhos[i], L, 5.4)
  vec_parts[i].energy(pot_tot)
  E_NM_pauli[i] = vec_parts[i].energy_pauli/vec_parts[i].n
  E_NM_kinetic[i] = vec_parts[i].kinetic*(otro_po/120)**2/vec_parts[i].n
params_fit = np.polyfit(rhos, E_NM_pauli + E_NM_kinetic, 2)
E_resto_fit = np.polyval(params_fit, rhos)
params_fit = np.polyfit(rhos, E_NM_pandha, 2)

def equiv_params_NM(params_fit):
  Eo = np.polyval(params_fit, 0.16)
  rho_min = -params_fit[1]/(2*params_fit[0])
  K = 2*params_fit[0]*(3*.16)**2
  return Eo, rho_min, K

Eo, rho_min, K = [-16.003, 0.16, 283]
E_NM_vec = Eo + K*(rhos/rho_min - 1)**2/18
V_E = np.var(E_NM_vec)
Chi2_VE = np.mean((E_NM_vec-E_NM_pandha)**2)/V_E

if (option=="test"):
  option_func = sys.argv[2]
  # Only one maximum at (0,0)
  function1 = lambda X: np.exp(-(X[0]**2+X[1]**2))
  # Discrete degeneration of maximum: at (0,0.1) and (0,-1.1)
  function2 = lambda X: np.exp(-9*(X[0]**2+(X[1]-0.1)**2)) + np.exp(-9*(X[0]**2+(X[1]+1.1)**2))
  # Continous degeneration of maximum: at x^2+y^2 = 1/4
  function3 = lambda X: np.exp(-9*(np.sqrt(X[0]**2+X[1]**2)-0.5)**2)
  # One maximum and a continous degeneration of maximum: at x^2+y^2 = 1/4 and (0,0)
  function4 = lambda X: np.exp(-25*(np.sqrt(X[0]**2+X[1]**2)-0.5)**2) + np.exp(-25*(X[0]**2+X[1]**2))
  # Two continous degeneration of maximum: at x^2+y^2 = 1/4 and x^2+y^2 = 1
  function5 = lambda X: np.exp(-25*(np.sqrt(X[0]**2+X[1]**2)-0.5)**2) + np.exp(-25*(np.sqrt(X[0]**2+X[1]**2)-1)**2)
  # Pandharipande fit for SC lattice
  l_bounds = [-1, -.5]
  u_bounds = [1, 1.5]
  if(option_func=="pandha" or option_func=="full"):
    def curva_energia_pandha(params):
      global pot_tot, parts, panda_nn, panda_np, L
      panda_np.V_a = params[0]*2666.647
      panda_np.V_r = params[0]*3088.118
      panda_nn.V_o = params[1]*373.118
      Es = np.zeros_like(rhos)
      for i in range(len(rhos)):
        vec_parts[i].energy(pot_tot)
        Es[i] = vec_parts[i].energy_panda/vec_parts[i].n
      return Es

    def functionpandha(params, alpha=2*np.pi/3):
      Es = curva_energia_pandha(params)
      Chi2 = np.mean((Es - E_NM_vec)**2)
      return np.sqrt((1.0 + Chi2_VE)/(1.0 + Chi2/(V_E)))

    def curva_energia_full(params):
      global pot_tot, parts, panda_nn, panda_np, L
      panda_np.V_a = params[0]*2666.647
      panda_np.V_r = params[0]*3088.118
      panda_nn.V_o = params[1]*373.118
      Es = np.zeros_like(rhos)
      for i in range(len(rhos)):
        vec_parts[i].energy(pot_tot)
        #Es[i] = vec_parts[i].energy_panda/vec_parts[i].n + E_NM_pauli[i] + E_NM_kinetic[i]
        Es[i] = vec_parts[i].energy_panda/vec_parts[i].n + E_resto_fit[i]
      return Es

    def functionfull(params, alpha=2*np.pi/3):
      Es = curva_energia_full(params)
      Chi2 = np.mean((Es - E_NM_vec)**2)
      return np.sqrt(1.0/(1.0 + Chi2/(V_E)))

    l_bounds = [0.5, 0.5]
    u_bounds = [2.5, 2.5]

  # Chosen function
  function = eval("function"+option_func)

  coding = GenCoding([10, 10], l_bounds, u_bounds)
  scheme = Schemata('Ranking', 2, 0.005, True)
  N = 20
  solutions = Poblacion(N, coding)
  F = solutions.eval_fitness(function)
  N_steps = 100
  mean_F = np.zeros(N_steps+1)
  max_F = np.zeros(N_steps+1)
  mean_F[0] = F
  max_F[0] = solutions.best_fitness

  fig1 = plt.figure(figsize=(23.8, 12.4))
  gs = gridspec.GridSpec(1, 2, width_ratios=[2,1])
  #ax1 = fig1.add_subplot(gs[0], autoscale_on=False, xlim=(-1, 1), ylim=(-1.5, 0.5))
  ax1 = fig1.add_subplot(gs[0], autoscale_on=False, xlim=(l_bounds[0], u_bounds[0]), ylim=(l_bounds[1], u_bounds[1]))
  ax1.grid()
  line2, = ax1.plot([], [], "ro", markersize=15)
  linefittest, = ax1.plot([], [], "g*", markersize=15)
  line1, = ax1.plot([], [], "bs")
  ax2 = fig1.add_subplot(gs[1], autoscale_on=False, xlim=(0, N_steps), ylim=(-0.05, 1.1))
  ax2.grid()
  ax2.set_xlabel("Generacion")
  ax2.set_ylabel("Fitness")
  line3, = ax2.plot([], [], "r-")
  line4, = ax2.plot([], [], "b-")
  ax2.legend(["Fitness promedio", "Mejor fitness"], loc = 'upper center')
  template = 'Generacion %d'
  text1 = ax1.text(0.35, 0.05, '', transform=ax1.transAxes, fontsize=20)
  if(option_func == "1"):
    line2.set_data([0], [0])
  if(option_func == "2"):
    line2.set_data([0, 0], [0.1, -1.1])
  theta = np.linspace(0, 2*np.pi, 100)
  if(option_func == "3"):
    line2.set_data(0.5*np.cos(theta), 0.5*np.sin(theta))
  if(option_func == "4"):
    line2.set_data(np.concatenate([0.5*np.cos(theta),[0]]), np.concatenate([0.5*np.sin(theta),[0]]))
  if(option_func == "5"):
    line2.set_data(np.concatenate([0.5*np.cos(theta),np.cos(theta)]), np.concatenate([0.5*np.sin(theta),np.sin(theta)]))

  def init():
    line1.set_data([], [])
    line3.set_data([], [])
    line4.set_data([], [])
    linefittest.set_data([], [])
    text1.set_text('')
    return line1, text1, line3, line4, linefittest

  def animate(i):
    global mean_F, max_F, solutions
    Ps = solutions.all_params
    Xs, Ys = Ps[0, :], Ps[1, :]
    line1.set_data(Xs, Ys)
    text1.set_text(template %(i+1))
    mean_F[i+1] = solutions.advance(scheme, function)
    max_F[i+1] = solutions.best_fitness
    line3.set_data(np.arange(i+1), mean_F[:i+1])
    line4.set_data(np.arange(i+1), max_F[:i+1])
    linefittest.set_data([solutions.best_params[0]], [solutions.best_params[1]])
    return line1, text1, line3, line4, linefittest

  anim = animation.FuncAnimation(fig1, animate, range(N_steps), interval=300, blit=True, init_func = init, repeat = False)
  plt.show()
  if(len(sys.argv)>3):
    Writer = animation.writers['avconv']
    writer = Writer(fps=3, metadata=dict(artist='Me'), bitrate=1800)
    anim.save('GA_conv'+option_func+'.mp4', writer=writer)

"""
Options are:
'fl,fu': Where fl,fu are floats; the bounds of the (free) parameter
'i': Where i is int; the index of the parameter it mimics
'=f': Where f is float; the fixed value of the parameter
Defaulted to 0.5,2
"""
idxs_free = []; l_bounds = []; u_bounds = []; idxs_fixed = [];
fixed_values = [];  idxs_dep = []; idxs_dep_of = []
for i in range(2, len(sys.argv)):
  data = (sys.argv[i]).split(',')
  if(len(data)==2):
    idxs_free.append(i-2)
    l_bounds.append(float(data[0]))
    u_bounds.append(float(data[1]))
  else:
    if(data[0][0]=="="):
      idxs_fixed.append(i-2)
      fixed_values.append(float(data[0][1:]))
    else:
      idxs_dep.append(i-2)
      idxs_dep_of.append(int(data[0][0]))
n_options = 8 + (option[:4]=="qcnm")
for i in range(len(sys.argv), n_options):
  idxs_free.append(i-2)
  l_bounds.append(0.5)
  u_bounds.append(2)
cant_params = len(idxs_free)
l_bounds = np.array(l_bounds)
u_bounds = np.array(u_bounds)
assert(set(idxs_dep_of).issubset(set(idxs_free)))

if (option=="posta"):

  def nombre_parametro(i):
    return "V"*(i<3) + "\mu"*(i>=3) + "_" + "a"*(i%3==0) + "r"*(i%3==1) + "o"*(i%3==2)

  def real_parameters(params):
    real_params = np.ones(6)
    real_params[idxs_fixed] = fixed_values
    real_params[idxs_free] = params
    real_params[idxs_dep] = real_params[idxs_dep_of]
    return real_params

  def set_params(params):
    real_params = real_parameters(params)
    if (len(params)==6):
      real_params = params
    panda_np.V_a = real_params[0]*2666.647
    panda_np.V_r = real_params[1]*3088.118
    panda_nn.V_o = real_params[2]*373.118
    panda_np.mu_a = 1.6 - real_params[3]
    panda_np.mu_r = 1.7468 - real_params[4]
    panda_nn.mu_o = 1.5 - real_params[5]
    rcut = 5.4
    panda_np.shift = (panda_np.V_r*np.exp(-panda_np.mu_r*rcut) - panda_np.V_a*np.exp(-panda_np.mu_a*rcut))/rcut
    panda_nn.shift = panda_nn.V_o*np.exp(-panda_nn.mu_o*rcut)/rcut

  def curva_energia_full(params):
    global pot_tot, parts, panda_nn, panda_np, L
    set_params(params)
    Es = np.zeros_like(rhos)
    for i in range(len(rhos)):
      vec_parts[i].energy(pot_tot)
      Es[i] = vec_parts[i].energy_panda/vec_parts[i].n + E_resto_fit[i]
    return Es

  def dame_params_NM(params):
    Es = curva_energia_full(params)
    params_fit_sol = np.polyfit(rhos, Es, 2)
    Eo_sol, rho_min_sol, K_sol = equiv_params_NM(params_fit_sol)
    return Eo_sol, rho_min_sol, K_sol

  tols = [1, 0.007, 300]
  def function(params):
    """
    Es = curva_energia_full(params)
    Chi2 = np.mean((Es - E_NM_vec)**2)
    return np.sqrt((1.0 + Chi2_VE)/(1.0 + Chi2/(V_E)))
    """
    Eo_sol, rho_min_sol, K_sol = dame_params_NM(params)
    new_Chi2 = ((Eo - Eo_sol)/tols[0])**2 + ((rho_min - rho_min_sol)/tols[1])**2 + ((400 - K_sol)/tols[2])**2 + 1000000*((K_sol < 100) + (K_sol > 1000))# + (Eo_sol < -20) + (rho_min_sol > 0.18) + (rho_min_sol < 0.14))
    return np.sqrt(1.0/(1.0 + new_Chi2))
    #"""

  def convergencia(chromosomes):
    return 1 - 2*np.mean(np.std(chromosomes, axis=0))

  def comparar(sols):
    params_NM = np.zeros((len(sols), 3))
    n = len(sols)
    colores = ["b", "g", "r", "c", "k", "y"]
    markers = ["o", "v", "^", "s", "d", "p", "<", "*", "X", ">", "D", "8"]
    plt.figure()
    plt.plot(rhos, E_NM_pandha, colores[0]+markers[0]+"-", markersize=10)
    for i,params in zip(range(n), sols):
      params_NM[i,:] = dame_params_NM(params)
      Es = curva_energia_full(params)
      #plt.plot(rhos, Es, colores[(i+1)%6]+markers[i+1]+"-", markersize=8)
      plt.plot(rhos, Es, marker="$%d$" %i, linestyle="-", markersize=10)
    plt.xlabel(r"$\rho$ [fm$^{-3}$]", fontsize=18)
    plt.ylabel("E [MeV]", fontsize=18)
    plt.xticks([0.15, 0.155, 0.16, 0.165, 0.17])
    plt.xlim(0.15, 0.17)
    #plt.legend(["Pandha puro $-$ $20$MeV"], loc="upper left")
    plt.figure()
    #plt.plot([rho_min], [K], "bo", markersize=15)
    plt.plot([rho_min], [Eo], "bo", markersize=15)
    for i in range(n):
      #plt.plot(params_NM[i,1], params_NM[i,2], marker="$%d$" %i, markersize=15)
      plt.plot(params_NM[i,1], params_NM[i,0], marker="$%d$" %i, markersize=15)
    #plt.axis([0.14, 0.165, 150, 290])
    plt.axis([0.155, 0.171, -17.4, -15.1])
    #plt.yticks(np.concatenate([params_NM[:,2], [K]]))
    #plt.yticks(np.concatenate([params_NM[:,0], [Eo]]))
    plt.xlabel(r"$\rho_o$ [fm$^{-3}$]", fontsize=18)
    #plt.ylabel(r"$\kappa$ [MeV]", fontsize=18)
    plt.ylabel(r"$E_o$ [MeV]", fontsize=18)
    plt.title("Ajustes", fontsize=18)
    plt.grid()
    plt.figure()
    ver_pandha_np([1, 1, 1, 0, 0, 0])
    for i in range(n):
      Rs, Vs = ver_pandha_np(sols[i])
      plt.plot(Rs[0:1000:100], Vs[0:1000:100], marker="$%d$" %i, markersize=10, linestyle="")
    plt.xlim(0.9, 5.4)
    plt.ylim(-45, 75)

  def comparar_redes(params):
    set_params(params)
    E = np.zeros((3, len(rhos)))
    leyenda = []
    for i in range(3):
      formato = eval("formato"+str(i+1))
      for j in range(len(rhos)):
        parts = Particles(1, 938)
        parts.load_lammpstrj(formato %rhos[j], L, 5.4)
        parts.energy(pot_tot)
        E[i,j] = parts.kinetic/parts.n + parts.energy_panda/parts.n + parts.energy_pauli/parts.n
      plt.plot(rhos, E[i,:], "o--")
      leyenda.append("B"+str(i+1))
    plt.plot(rhos, E_NM_pandha, "k-")
    leyenda.append("Pandha puro")
    plt.legend(leyenda)
    plt.xlabel(r"Densidad [fm$^{-3}$]")
    plt.ylabel("Energia [MeV]")

  def ver_pandha_np(params, ro = 1, n = 1000, rc = 5.4):
    Rs = np.linspace(ro, rc, n)
    real_params = real_parameters(params)
    V_a = real_params[0]*2666.647
    V_r = real_params[1]*3088.118
    #V_o = real_params[2]*373.118
    mu_a = 1.6 - real_params[3]
    mu_r = 1.7468 - real_params[4]
    #mu_o = 1.5 - real_params[5]
    Vc = (V_r*np.exp(-mu_r*rc) - V_a*np.exp(-mu_a*rc))/rc
    Vs = (V_r*np.exp(-mu_r*Rs) - V_a*np.exp(-mu_a*Rs))/Rs - Vc
    plt.plot(Rs, Vs)
    print(Rs[np.where(Vs == np.min(Vs))])
    plt.show()
    return Rs, Vs


  coding = GenCoding(11*np.ones(cant_params, dtype=np.int), l_bounds, u_bounds)
  mutacion = 0.05/cant_params
  scheme = Schemata('Ranking', 2, mutacion, True)
  N = 20*cant_params
  solutions = Poblacion(N, coding)
  F = solutions.eval_fitness(function)
  N_steps = 100
  mean_F = np.zeros(N_steps+1)
  max_F = np.zeros(N_steps+1)
  conv = np.zeros(N_steps+1)
  mean_F[0] = F
  max_F[0] = solutions.best_fitness

  fig1 = plt.figure(figsize=(23.8, 12.4))
  ax1 = fig1.add_subplot(121, autoscale_on=False, xlim=(0.145, 0.175), ylim=(-20, -12))
  ax1.grid()
  ax1.set_xlabel("Densidad")
  ax1.set_ylabel("Energia")
  lineopt, = ax1.plot(rhos, E_NM_pandha, "ro--")
  linefittest, = ax1.plot([], [], "bo--")
  ax1.legend(["Pandha puro", "Mejor ajuste Pandha"], loc = 'upper left')
  GLs = 'GL: '
  for i in idxs_free:
    GLs += '$' + nombre_parametro(i) + '$|'
  valores = ''
  for i in range(6):
    valores += '$' + nombre_parametro(i) +  ' *'*(i<3) + ' -'*(i>=3) + ' %.3f' + '$       ' + '\n'*(i%3==2)
  valores += '\n'
  template = valores+GLs+'\t Generacion %d'
  text1 = ax1.text(0.2, 0.03, '', transform=ax1.transAxes, fontsize=20)
  template2 = r'$E_o=%.2f$MeV   $\rho_o=%.3f$fm$^{-3}$   $\kappa=%.1f$MeV'
  text2 = ax1.text(0.5, 0.9, '', transform=ax1.transAxes, fontsize=20, ha="center", va="center")
  ax2 = fig1.add_subplot(122, autoscale_on=False, xlim=(0, N_steps), ylim=(-0.05, 1.1))
  ax2.grid()
  ax2.set_xlabel("Generacion")
  ax2.set_ylabel("Fitness")
  line3, = ax2.plot([], [], "r-")
  line4, = ax2.plot([], [], "b-")
  line5, = ax2.plot([], [], "g-")
  ax2.legend(["Fitness promedio", "Mejor fitness", "Convergencia"], loc = 'upper center')

  def init():
    line3.set_data([], [])
    line4.set_data([], [])
    line5.set_data([], [])
    linefittest.set_data([], [])
    text1.set_text('')
    text2.set_text('')
    return text1, text2, line3, line4, linefittest, line5

  def animate(i):
    global mean_F, max_F, solutions
    mean_F[i+1] = solutions.advance(scheme, function)
    max_F[i+1] = solutions.best_fitness
    conv[i+1] = convergencia(solutions.chromosomes)
    real_fittest = list(real_parameters(solutions.best_params))
    params_NM_fittest = list(dame_params_NM(real_fittest))
    real_fittest.append(i+1)
    text1.set_text(template %tuple(real_fittest))
    text2.set_text(template2 %tuple(params_NM_fittest))
    linefittest.set_data(rhos, curva_energia_full(solutions.best_params))
    line3.set_data(np.arange(i+1), mean_F[:i+1])
    line4.set_data(np.arange(i+1), max_F[:i+1])
    line5.set_data(np.arange(i+1), conv[:i+1])
    return text1, text2, line3, line4, linefittest, line5

  anim = animation.FuncAnimation(fig1, animate, range(N_steps), interval=300, blit=True, init_func = init, repeat = False)
  plt.show()

#-------------------------------------------------------------------------------#
if (option[:4]=="qcnm"):

  def nombre_parametro(i):
    if(i == 0):
      return "V_o"
    elif ((i-1)//2 == 0):
      return "p_" + str(i)
    elif ((i-1)//2 == 1):
      return "r_" + str(i-2)
    elif (i == 5):
      return "d"
    elif (i == 6):
      return "a"
    return ""

  def real_parameters(params):
    real_params = np.ones(7)
    real_params[idxs_fixed] = fixed_values
    real_params[idxs_free] = params
    real_params[idxs_dep] = real_params[idxs_dep_of]
    return real_params

  def set_params(params):
    global pot_tot, vec_parts, panda_nn, panda_np, L, qcnm
    real_params = real_parameters(params)
    if (len(params)==7):
      real_params = params
    qcnm.V_o = 25.93*real_params[0]
    qcnm.p_1 = 6.2 - real_params[1]
    qcnm.p_2 = 3.0 - real_params[2]
    qcnm.r_1 = 1.757 - real_params[3]
    qcnm.r_2 = 1.771 - real_params[4]
    qcnm.d = 3.35 - real_params[5]
    qcnm.a = 5./6 - real_params[6]
    qcnm.rcut = 5.4
    qcnm.actualizar_shift()

  def curva_energia_full(params):
    global pot_tot, vec_parts, panda_nn, panda_np, L, qcnm
    set_params(params)
    Es = np.zeros_like(rhos)
    for i in range(len(rhos)):
      vec_parts[i].energy_QCNM(pot_tot)
      Es[i] = vec_parts[i].energy_panda/vec_parts[i].n + (len(option)==4)*E_resto_fit[i]
    return Es

  def dame_params_NM(params):
    Es = curva_energia_full(params)
    params_fit_sol = np.polyfit(rhos, Es, 2)
    Eo_sol, rho_min_sol, K_sol = equiv_params_NM(params_fit_sol)
    return Eo_sol, rho_min_sol, K_sol

  tols = [1, 0.007, 200]
  def function(params):
    """
    Es = curva_energia_full(params)
    Chi2 = np.mean((Es - E_NM_vec)**2)
    return np.sqrt((1.0 + Chi2_VE)/(1.0 + Chi2/(V_E)))
    """
    Eo_sol, rho_min_sol, K_sol = dame_params_NM(params)
    new_Chi2 = ((Eo - Eo_sol)/tols[0])**2 + ((rho_min - rho_min_sol)/tols[1])**2 + ((300 - K_sol)/tols[2])**2 + 1000000*((K_sol < 100) + (K_sol > 1000))# + (Eo_sol < -20) + (rho_min_sol > 0.18) + (rho_min_sol < 0.14))
    return np.sqrt(1.0/(1.0 + new_Chi2))
    #"""

  def convergencia(chromosomes):
    return 1 - 2*np.mean(np.std(chromosomes, axis=0))

  def comparar_Evsrho(sols):
    formato_params = r"$E_o=%.2f$MeV | $\rho_o=%.3f$fm$^{-3}$ | $\kappa=%.0f$MeV"
    plt.plot(rhos, E_NM_pandha, "-")
    #plt.text(0.16, np.min(E_NM_pandha), formato_params %(Eo, rho_min, K), ha="center", va="top")
    params_NM = []
    for params in sols:
      Es = curva_energia_full(params)
      params_NM.append(dame_params_NM(params))
      plt.plot(rhos, Es, "o--")
      plt.text(params_NM[-1][1], np.min(Es), formato_params %tuple(params_NM[-1]), ha="center", va="top", fontsize=14)
    plt.xlabel(r"$\rho$ [fm$^{-3}$]", fontsize=18)
    plt.ylabel("E [MeV]", fontsize=18)
    plt.xticks([0.15, 0.155, 0.16, 0.165, 0.17])
    plt.xlim(0.15, 0.17)
    #plt.legend(["Pandha puro"])
    plt.legend([formato_params %(Eo, rho_min, K)])
    plt.show()

  def comparar_redes(params):
    set_params(params)
    E = np.zeros((3, len(rhos)))
    leyenda = []
    for i in range(3):
      formato = eval("formato"+str(i+1))
      for j in range(len(rhos)):
        parts = Particles(1, 938)
        parts.load_lammpstrj(formato %rhos[j], L, 5.4)
        parts.energy_QCNM(pot_tot)
        E[i,j] = parts.kinetic/parts.n + parts.energy_panda/parts.n + parts.energy_pauli/parts.n
      plt.plot(rhos, E[i,:], "o--")
      leyenda.append("B"+str(i+1))
    plt.plot(rhos, E_NM_pandha, "k-")
    plt.xlim(0.15, 0.17)
    leyenda.append("Pandha puro")
    plt.legend(leyenda, loc="center right")
    plt.xlabel(r"Densidad [fm$^{-3}$]")
    plt.ylabel("Energia [MeV]")
    return E


  def ver_QCNM(params, ro = 1, n = 1000, rc = 5.4):
    Rs = np.linspace(ro, rc, n)
    real_params = real_parameters(params)
    qcnm_temp = QCNM(real_params[0]*25.93, 6.2 - real_params[1], 3.0 - real_params[2], 1.757 - real_params[3], 1.771 - real_params[4], 3.35 - real_params[5], 5./6 - real_params[6], 5.4)
    V = np.array([qcnm_temp.potential(r) for r in Rs])
    plt.plot(Rs, V, "-")
    V_min = np.min(V)
    r_min = Rs[np.where(V == V_min)]
    plt.plot([r_min], [V_min], "o")
    plt.xlabel("r [fm]")
    plt.ylabel("QCNM [MeV]")
    print(r_min, V_min)
    plt.show()


  coding = GenCoding(10*np.ones(cant_params, dtype=np.int), l_bounds, u_bounds)
  mutacion = 0.05/cant_params
  scheme = Schemata('Ranking', 2, mutacion, True)
  N = 20*cant_params
  solutions = Poblacion(N, coding)
  F = solutions.eval_fitness(function)
  N_steps = 100
  mean_F = np.zeros(N_steps+1)
  max_F = np.zeros(N_steps+1)
  conv = np.zeros(N_steps+1)
  mean_F[0] = F
  max_F[0] = solutions.best_fitness

  fig1 = plt.figure(figsize=(23.8, 12.4))
  ax1 = fig1.add_subplot(121, autoscale_on=False, xlim=(0.145, 0.175), ylim=(-20, -12))
  ax1.grid()
  ax1.set_xlabel("Densidad")
  ax1.set_ylabel("Energia")
  lineopt, = ax1.plot(rhos, E_NM_pandha, "ro--")
  linefittest, = ax1.plot([], [], "bo--")
  ax1.legend(["Pandha puro", "Mejor ajuste QCNM"], loc = 'upper left')
  GLs = 'GL: '
  for i in idxs_free:
    GLs += '$' + nombre_parametro(i) + '$|'
  valores = ''
  for i in range(7):
    valores += '$' + nombre_parametro(i) +  ' *'*(i==0) + ' -'*(i!=0) + ' %.3f' + '$       ' + '\n'*(i%3==2)
  valores += '\n'
  template = valores+GLs+'\t Generacion %d'
  text1 = ax1.text(0.2, 0.03, '', transform=ax1.transAxes, fontsize=20)
  template2 = r'$E_o=%.2f$MeV   $\rho_o=%.3f$fm$^{-3}$   $\kappa=%.1f$MeV'
  text2 = ax1.text(0.5, 0.9, '', transform=ax1.transAxes, fontsize=20, ha="center", va="center")
  ax2 = fig1.add_subplot(122, autoscale_on=False, xlim=(0, N_steps), ylim=(-0.05, 1.1))
  ax2.grid()
  ax2.set_xlabel("Generacion")
  ax2.set_ylabel("Fitness")
  line3, = ax2.plot([], [], "r-")
  line4, = ax2.plot([], [], "b-")
  line5, = ax2.plot([], [], "g-")
  ax2.legend(["Fitness promedio", "Mejor fitness", "Convergencia"], loc = 'upper center')

  def init():
    line3.set_data([], [])
    line4.set_data([], [])
    line5.set_data([], [])
    linefittest.set_data([], [])
    text1.set_text('')
    text2.set_text('')
    return text1, text2, line3, line4, linefittest, line5

  def animate(i):
    global mean_F, max_F, solutions
    mean_F[i+1] = solutions.advance(scheme, function)
    max_F[i+1] = solutions.best_fitness
    conv[i+1] = convergencia(solutions.chromosomes)
    real_fittest = list(real_parameters(solutions.best_params))
    params_NM_fittest = list(dame_params_NM(real_fittest))
    real_fittest.append(i+1)
    text1.set_text(template %tuple(real_fittest))
    text2.set_text(template2 %tuple(params_NM_fittest))
    linefittest.set_data(rhos, curva_energia_full(solutions.best_params))
    line3.set_data(np.arange(i+1), mean_F[:i+1])
    line4.set_data(np.arange(i+1), max_F[:i+1])
    line5.set_data(np.arange(i+1), conv[:i+1])
    return text1, text2, line3, line4, linefittest, line5

  anim = animation.FuncAnimation(fig1, animate, range(N_steps), interval=300, blit=True, init_func = init, repeat = False)
  plt.show()
