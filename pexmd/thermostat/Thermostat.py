import numpy as np
import itertools as it
import ctypes as ct

class Thermostat():
  """
  Base Thermostat class.
  """
  def __init__(self, T):
      self._T = T

class Andersen(Thermostat):

  def __init__(self, temperature, mass, seed=None):
    self.mass = mass
    self.rng = np.random.RandomState(seed)
    self._seed = seed
    self.temperature = temperature

  @property
  def temperature(self):
    return self._temperature

  @temperature.setter
  def temperature(self, value):
    """
    Set temperature of thermostat.
    Parameters
    ----------
    value : float
        The new temperature
    """
    if value > 0:
      self._temperature = value
    else:
      msg = "Trying to set negative or zero ({0}) temperature"
      raise ValueError(msg.format(value))

  @property
  def mass(self):
    return self._mass

  @mass.setter
  def mass(self, value):
    """
    Set mass of thermostat
    Parameters
    ----------
    value : float
        The new mass
    """
    if value >= 0 and value <= 1:
      self._mass = value
    else:
      msg = "Trying to set mass {0} outside range [0,1]"
      raise ValueError(msg.format(value))

  def step(self, v, m):
    N = len(v)
    s = np.sqrt(self.temperature/m)
    for i in range(N):
      if (self.rng.random_sample() <= self.mass):
        v[i,:] = self.rng.normal(loc = 0, scale = s[i], size = 3)
    return v
