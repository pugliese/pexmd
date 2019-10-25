import numpy as np

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
  def __init__(self, selection, crossover, mutation, elitist=False, crossover_rate=0.8):
    """
    Parameters
    ----------
    selection: {'FairTournament', 'RandomTournament', 'FPS', 'Ranking'} for now
        Selection method for pool parents

    crossover: int
          Number of crossovers between parents

    mutation: float
        Probability of gen 'flipping' (for each gen), tipically 0.001

    crossover_rate: float
          Percentage of renewal per generation (successful crossovers)
    """
    self.selection = selection
    self.crossover = crossover
    self.crossover_rate = crossover_rate
    self.mutation = mutation
    self.elitist = elitist

  def select_parents(self, fitness):
    """
    Return pool of parents selected (indexes) for reproduction.
    Depends on self.selection.
    -----------
    In a tournament, we choose pairs of individuals and add to the pool the
    fittest of them. Fair and Random defines the criteria used for pairing
    Random: Random pairs with replacement
    Fair: Random pairs without replacement, 2 tournaments
    -----------
    In FPS and Ranking, we assing a probability of being added to the pool
    based on the fitness of the individual.
    FPS: Probability proportional to fitness
    Ranking: Fixed probability based on fitness order
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
    Also depends on self.crossover_rate, the chance of actual crossover of
    parents (if not, they replicate themselves exactly)
    """
    N = len(parents)
    N_genes = len(chromosomes[0,:])
    assert(N % 2 == 0)
    new_chromosomes = np.zeros_like(chromosomes)
    np.random.shuffle(parents)
    for i in range(N//2):
      ps = [parents[2*i], parents[2*i+1]]
      if(np.random.random() < self.crossover_rate):
        cuts = np.concatenate([[0], np.sort(np.random.randint(0, N, self.crossover)), [N]])
      else:
        cuts = np.array([0,N])
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
    N = len(chromosomes[:,0])
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
