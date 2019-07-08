import numpy as np
import matplotlib.pylab as plt
import matplotlib.animation as animation
import matplotlib.gridspec as gridspec
import ctypes as ct
import sys
import copy
plt.ion()

# ---------------------------------------------------------------------------- #
# ------------------------ Pandha+Pauli -------------------------------------- #
# ---------------------------------------------------------------------------- #
C_funcs = ct.CDLL('../Panda_Pauli/lib.so')

class Pauli(ct.Structure):
  """
  Pauli interaction class to match struct from C
  """
  _fields_ = [("qo", ct.c_float), ("po", ct.c_float), ("D", ct.c_float), ("scut2", ct.c_float), ("shift", ct.c_float), ("LUT", ct.c_voidp), ("ds2", ct.c_float)]

  def __init__(self, qo, po, D, scut2):
    self.qo = qo
    self.po = po
    self.D = D
    self.scut2 = scut2
    self.shift = D*np.exp(-0.5*scut2)

class Panda_np(ct.Structure):
  """
  Pandharipande potential for neutron-proton
  interaction to match struct from C
  """
  _fields_ = [("mu_r", ct.c_float), ("mu_a", ct.c_float), ("V_r", ct.c_float), ("V_a", ct.c_float), ("rcut", ct.c_float), ("rcut2", ct.c_float),
            ("shift", ct.c_float), ("LUT", ct.c_voidp), ("dr2", ct.c_float)]

  def __init__(self, mu_r, mu_a, V_r, V_a, rcut):
    self.mu_r = mu_r
    self.mu_a = mu_a
    self.V_r = V_r
    self.V_a = V_a
    self.rcut = rcut
    self.rcut2 = rcut**2
    self.shift = (V_r*np.exp(-mu_r*rcut) - V_a*np.exp(-mu_a*rcut))/rcut

class Panda_nn(ct.Structure):
  """
  Pandharipande potential for neutron-neutron or
  proton-proton interaction to match struct from C
  """
  _fields_ = [("mu_o", ct.c_float), ("V_o", ct.c_float), ("rcut", ct.c_float), ("rcut2", ct.c_float), ("shift", ct.c_float), ("LUT", ct.c_voidp), ("dr2", ct.c_float)]

  def __init__(self, mu_o, V_o, rcut):
    self.mu_o = mu_o
    self.V_o = V_o
    self.rcut = rcut
    self.rcut2 = rcut**2
    self.shift = V_o*np.exp(-mu_o*rcut)/rcut

class TotalPotential(ct.Structure):
  """
  Total potential of interaction to match struct from C
  """
  _fields_ = [("pauli", ct.POINTER(Pauli)), ("panda_nn", ct.POINTER(Panda_nn)), ("panda_np", ct.POINTER(Panda_np)),  ("rcut", ct.c_float)]

  def __init__(self, pauli, panda_nn, panda_np):
    self.pauli = ct.pointer(pauli)
    self.panda_nn = ct.pointer(panda_nn)
    self.panda_np = ct.pointer(panda_np)
    self.rcut = max([np.sqrt(pauli.scut2)*pauli.qo, panda_nn.rcut, panda_np.rcut])

load_lammpstrj_c = C_funcs.load_lammpstrj
load_lammpstrj_c.argtypes = [ct.c_char_p, ct.c_voidp, ct.c_voidp, ct.c_float]
load_lammpstrj_c.restype = ct.c_int

energia_sin_LUT_c = C_funcs.energia_sin_LUT
energia_sin_LUT_c.argtypes = [ct.c_voidp, ct.c_voidp]
energia_sin_LUT_c.restype = ct.c_int

class Particles(ct.Structure):
  """
  Particles class to match struct from C
  """
  _fields_ = [("n", ct.c_int), ("type", ct.c_voidp), ("q", ct.c_voidp), ("p", ct.c_voidp), ("mass", ct.c_float), ("siguiente", ct.c_voidp),
                ("anterior", ct.c_voidp), ("primero", ct.c_voidp), ("celda", ct.c_voidp), ("M", ct.c_int), ("l", ct.c_float),
                ("energy_pauli", ct.c_float), ("energy_panda", ct.c_float), ("kinetic", ct.c_float)]

  def __init__(self, n, mass):
    self.n = n
    self.mass = np.float32(mass)
    self.q = np.zeros(3*n, dtype = np.float32).ctypes.data_as(ct.c_voidp)
    self.p = np.zeros(3*n, dtype = np.float32).ctypes.data_as(ct.c_voidp)
    self.kinetic = 0
    self.energy_panda = 0
    self.energy_pauli = 0

  def load_lammpstrj(self, filename, L, rcut):
    filename_p = ct.c_char_p(filename.encode('utf-8'))
    #filename_p = ct.pointer(filename.encode('utf-8'))
    L_p = L.ctypes.data_as(ct.c_voidp)
    rcut_32 = np.float32(rcut)
    return load_lammpstrj_c(filename_p, ct.pointer(self), L_p, rcut_32)

  def energy(self, pot_tot):
    return energia_sin_LUT_c(ct.pointer(self), ct.pointer(pot_tot))

# ---------------------------------------------------------------------------- #
# ------------------------ GA specific --------------------------------------- #
# ---------------------------------------------------------------------------- #
class GenCoding():
  """
  Class for chromosome representation
  """
  def __init__(self, n_bits, lower_bounds, upper_bounds):
    """
    Parameters
    ----------
    n_bits: 1D array of ints
        Number of bits for the representation of each parameter of the solution

    lower_bounds: 1D array of floats
        Lower bound of parameters

    upper_bounds: 1D array of floats
        Upper bound of parameters
    """
    if(np.isscalar(n_bits)):
      self.n_params = 1
    else:
      self.n_params = len(n_bits)
    if len(lower_bounds) != self.n_params:
      msg = "Trying to set {0} lower bounds for {1} parameters"
      raise ValueError(msg.format(len(lower_bounds), self.n_params))
    if len(upper_bounds) != self.n_params:
      msg = "Trying to set {0} upper bounds for {1} parameters"
      raise ValueError(msg.format(len(upper_bounds), self.n_params))
    self.bits = n_bits
    self.lower = lower_bounds
    self.upper = upper_bounds

  def params(self, chromosome):
    """
    Returns the parameter values for a given chromosome

    Parameters
    ----------
    chromosome: 1D array of ints
        Array of 0's and 1's for the bit representation of parameters
    """
    if len(chromosome) != np.sum(self.bits):
      msg = "Chromosome of {0} bits for an {1} bit representation of {2} parameters"
      raise ValueError(msg.format(len(upper_bounds), np.sum(self.bits), self.n_params))
    res = np.zeros(self.n_params)
    idx = 0
    for i in range(self.n_params):
      for j in range(self.bits[i]):
        res[i] = (res[i] + chromosome[idx+j])*0.5
      res[i] = res[i]*(self.upper[i]-self.lower[i])+self.lower[i]
      idx += self.bits[i]
    return res

class Schemata():
  """
  Class for schemata representation: selection, crossover and mutation
  """
  def __init__(self, selection, crossover, mutation, elitist=False):
    """
    Parameters
    ----------
    selection: {'FairTournament', 'RandomTournament', 'FPS', 'Ranking'} for now
        Selection method for pool parents

    crossover: int
          Number of crossovers between parents

    mutation: float
        Probability of gen 'flipping' (for each gen), tipically 0.001
    """
    self.selection = selection
    self.crossover = crossover
    self.mutation = mutation
    self.elitist = elitist

  def select_parents(self, fitness):
    """
    Return pool of parent selected (indexes) for reproduction
    Depends on self.selection
    -----------
    In a tournament, we choose pairs of individuals and add to the pool the
    fittest of them. Fair and Random defines the criteria used for pairing
    Random: Random pairs with replacement
    Fair: Random pairs without replacement, 2 tournaments
    -----------
    In FPS and Ranking, we assing a probability of being added to the pool
    based on the fitness of the individual.
    FPS: Probability proportional to fitness
    Ranking: Fixed probability based of fitness order
    """
    N = len(fitness)
    pool = np.zeros(N, dtype=np.int)
    if(self.selection=='FairTournament' or self.selection=='RandomTournament'):
      if(self.selection=='RandomTournament'):
        matchings = np.random.randint(0, N, 2*N)
        for k in range(N):
          if (matchings[2*k] == matchings[2*k+1]):
            new = np.random.randint(0, N-1)
            matchings[2*k+1] = new + (new >= matchings[2*k])
      elif(self.selection=='FairTournament'):
        idxs1 = np.arange(0, N)
        np.random.shuffle(idxs1)
        idxs2 = np.arange(0, N)
        np.random.shuffle(idxs2)
        matchings = np.concatenate([idxs1, idxs2])
      for k in range(N):
        ii = matchings[2*k]
        jj = matchings[2*k+1]
        if (fitness[ii] > fitness[jj]):
          pool[k] = ii
        else:
          pool[k] = jj
    elif (self.selection=='FPS' or self.selection=='Ranking'):
      if (self.selection=='FPS'):
        probs = fitness/np.sum(fitness)
      elif (self.selection=='Ranking'):
        temp = np.argsort(fitness)
        ranks = np.empty_like(temp)
        ranks[temp] = np.arange(N)
        weights = 1.0/(N-ranks)
        probs = weights/np.sum(weights)
      dP = 1.0/N
      start = np.random.uniform(0, dP)
      suma = 0.0
      i = -1
      for k in range(N):
        while (suma < (start + k*dP)):
          suma += probs[i]
          i += 1
        pool[k] = i
    return pool

  def children(self, parents, chromosomes):
    """
    Given a set of parents (N indexes) with ther corresponding chromosomes
    (2D array of N x total_bits), creates unmutated offspring.
    Depends on self.crossover; the number of subchromosomes in which each
    parent is separated to recombine (alternating) for the children.
    """
    N = len(parents)
    N_genes = len(chromosomes[0,:])
    assert(N % 2 == 0)
    new_chromosomes = np.zeros_like(chromosomes)
    np.random.shuffle(parents)
    for i in range(N//2):
      ps = [parents[2*i], parents[2*i+1]]
      cuts = np.concatenate([[0], np.sort(np.random.randint(0, N, self.crossover)), [N]])
      for g in range(len(cuts)-1):
        new_chromosomes[2*i, cuts[g]:cuts[g+1]] = chromosomes[ps[g%2], cuts[g]:cuts[g+1]]
        new_chromosomes[2*i+1, cuts[g]:cuts[g+1]] = chromosomes[ps[(g+1)%2], cuts[g]:cuts[g+1]]
    return new_chromosomes

  def mutate(self, chromosomes):
    """
    Returns the mutated chromosomes.
    Each gene (bit) mutates with chance self.mutation
    """
    N = len(chromosomes[:,0])
    N_genes = len(chromosomes[0,:])
    new_chromosomes = np.zeros_like(chromosomes)
    for i in range(N):
      for g in range(N_genes):
        p = np.random.uniform()
        new_chromosomes[i, g] = (chromosomes[i, g]+(p < self.mutation))%2
    return new_chromosomes

  def new_generation(self, fitness, chromosomes, fittest = 0):
    """
    Performs the 3 steps from above and returns a new generation.
    If Elitist, a random child is replaced with the actual fittest.
    Possible improvement: replace worst child?
    """
    pool = self.select_parents(fitness)
    new_chromosomes = self.children(pool, chromosomes)
    new_chromosomes = self.mutate(new_chromosomes)
    if(self.elitist):
      idx = np.random.randint(0,N)
      new_chromosomes[idx, :] = chromosomes[fittest, :]
    return new_chromosomes

class Poblacion():
  """
  Class for population (GenCoding+chromosomes)
  """
  def __init__(self, N, coding):
    """
    Parameters
    ----------
    N: int
        Number of individuals (chromosomes)

    coding: class GenCoding
          GenCoding of the individuals for chromosome interpretation
    """
    total_bits = np.sum(coding.bits)
    self.size = N
    self.chromosomes = np.random.randint(0, 2, (N, total_bits))
    self.coding = coding
    self.fitness = np.zeros(N)
    self.fittest = 0

  def eval_fitness(self, fitness_function):
    """
    Evaluates fitness for each individual and returns
    total fitness (for mean fitness calculation)
    """
    for i in range(self.size):
      self.fitness[i] = fitness_function(self.coding.params(self.chromosomes[i,:]))
    self.fittest = np.argmax(self.fitness)
    return np.mean(self.fitness)

  def advance(self, schemata, function):
    """
    Advances to the new generation and recalculates fitness, returns mean
    """
    self.chromosomes = schemata.new_generation(self.fitness, self.chromosomes, self.fittest)
    return self.eval_fitness(function)

  def params(self, i):
    """
    Returns parameters of a single individual (solution) i
    """
    return self.coding.params(self.chromosomes[i,:])

  @property
  def all_params(self):
    """
    Returns the whole population parameters
    """
    res = np.zeros((self.coding.n_params, self.size))
    for i in range(self.size):
      res[:, i] = self.params(i)
    return res

  @property
  def best_params(self):
    """
    Return the fittest parameters
    """
    return self.params(self.fittest)

  @property
  def best_fitness(self):
    """
    Returns the fittest fitness
    """
    return self.fitness[self.fittest]

  @property
  def mean_fitness(self):
    """
    Returns the mean fitness
    """
    return np.mean(self.fitness)

option = sys.argv[1]

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
  if(option_func=="pandha"):
    pauli = Pauli(qo = 1.644, po = 120, D = 207*(option_func!="pandha"), scut2 = (5.4/1.644)**2)
    panda_nn = Panda_nn(1.5, 373.118, 5.4)
    panda_np = Panda_np(mu_r = 1.7468, mu_a = 1.6, V_r = 3088.118, V_a = 2666.647, rcut = 5.4)
    pot_tot = TotalPotential(pauli, panda_nn, panda_np)

    rhos = np.array([0.140, 0.150, 0.153, 0.155, 0.158, 0.160, 0.162, 0.165, 0.167, 0.170, 0.180])
    rhos = np.array([0.140, 0.150, 0.160, 0.170, 0.180])
    E_NM_pandha = np.zeros_like(rhos)
    vec_parts = [0 for rho in rhos]
    L = np.array([3.0], dtype = np.float32)
    for i in range(len(rhos)):
      vec_parts[i] = Particles(64, 938)
      vec_parts[i].load_lammpstrj("fundamental/config_1728_x1_%.3f_0.010.lammpstrj" %rhos[i], L, 5.4)
#      vec_parts[i].load_lammpstrj("fundamental/config_redondear_%.3f_0.100.lammpstrj" %rhos[i], L, 5.4)
      vec_parts[i].energy(pot_tot)
      E_NM_pandha[i] = vec_parts[i].energy_panda/vec_parts[i].n
    params_fit = np.polyfit(rhos, E_NM_pandha, 2)
    E_NM = lambda rho: np.polyval(params_fit, rhos)
    E_NM_vec = E_NM(rhos)
    V_E = np.var(E_NM_vec)
    sqrtChi2_VE = np.sqrt(np.mean((E_NM_vec-E_NM_pandha)**2)/V_E)

    def curva_energia(params):
      global pot_tot, parts, panda_nn, panda_np, L
      panda_np.V_a = params[0]*2666.647
      panda_np.V_r = params[0]*3088.118
      panda_nn.V_o = params[1]*373.118
      Es = np.zeros_like(rhos)
      for i in range(len(rhos)):
        vec_parts[i].energy(pot_tot)
        Es[i] = vec_parts[i].energy_panda/vec_parts[i].n
      return Es

    def functionpandha(params):
      Es = curva_energia(params)
      Chi2 = np.mean((Es - E_NM_vec)**2)
      return (1.0 + sqrtChi2_VE)/(1.0 + np.sqrt(Chi2/(V_E)))

  # Chosen function
  function = eval("function"+option_func)

  #coding = GenCoding([10,10], [-1, -1.5], [1, 0.5])
  coding = GenCoding([10, 10], [0.5, 0.5], [2.5, 2.5])
  scheme = Schemata('Ranking', 2, 0.001, True)
  N = 100
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
  ax1 = fig1.add_subplot(gs[0], autoscale_on=False, xlim=(0.5, 2.5), ylim=(0.5, 2.5))
  ax1.grid()
  line2, = ax1.plot([], [], "ro", markersize=15)
  linefittest, = ax1.plot([], [], "g*", markersize=15)
  line1, = ax1.plot([], [], "bs")
  ax2 = fig1.add_subplot(gs[1], autoscale_on=False, xlim=(0, N_steps-1), ylim=(0, 1.1))
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
pauli = Pauli(qo = 1.644, po = 120, D = 207, scut2 = (5.4/1.644)**2)
panda_nn = Panda_nn(1.5, 373.118, 5.4)
panda_np = Panda_np(mu_r = 1.7468, mu_a = 1.6, V_r = 3088.118, V_a = 2666.647, rcut = 5.4)
pot_tot = TotalPotential(pauli, panda_nn, panda_np)

parts = Particles(64, 938)
L = np.array([3.0], dtype = np.float32)
parts.load_lammpstrj("fundamental/config_1728_x1_0.160_0.010.lammpstrj", L, 5.4)
A = ct.cast(parts.p, ct.POINTER(ct.c_float))
#print(0.5*np.sum(np.array([A[i] for i in range(3*parts.n)])**2/938)/parts.n)
print(parts.kinetic/parts.n, parts.energy_panda/parts.n, parts.energy_pauli/parts.n)
a = parts.energy(pot_tot)
print(parts.kinetic/parts.n, parts.energy_panda/parts.n, parts.energy_pauli/parts.n)
#print(a)
"""
