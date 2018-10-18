"""
Main Particles module.
"""

import numpy as np

class Base(object):
  """
  Base Particles class. It is abstract and we should specify which
  type of particle we actually want in order to fill it
  """
  def __init__(self, n):
    self.n = n
    self._x = np.zeros((n, 3), dtype=np.float32)
    self._v = np.zeros((n, 3), dtype=np.float32)
    self._f = np.zeros((n, 3), dtype=np.float32)
    self._t = np.zeros(n, dtype=np.int32)
    self._mass = np.zeros(n, dtype=np.float32)
    self.idx = np.arange(n)

  @property
  def x(self):
    return self._x

  @x.setter
  def x(self, value):
    """
    Set positions of particles.

    Parameters
    ----------

    value : 2D NumPy array
        The new positons of particles in an Nx3 array
    """
    value = np.array(value, dtype=np.float32)
    number = np.shape(value)[0]
    if self.n == number:
      self._x = value
    else:
      msg = "Trying to set {0} positions for a system with {1} particles"
      raise ValueError(msg.format(number, self.n))

  @property
  def v(self):
    return self._v

  @v.setter
  def v(self, value):
    """
    Set velocities of particles.

    Parameters
    ----------

    value : 2D NumPy array
        The new velocities of particles in an Nx3 array
    """
    value = np.array(value, dtype=np.float32)
    number = np.shape(value)[0]
    if self.n == number:
      self._v = value
    else:
      msg = "Trying to set {0} velocities for a system with {1} particles"
      raise ValueError(msg.format(number, self.n))

  @property
  def f(self):
    return self._f

  @f.setter
  def f(self, value):
    """
    Set forces of particles.

    Parameters
    ----------

    value : 2D NumPy array
        The new forces of particles in an Nx3 array
    """
    value = np.array(value, dtype=np.float32)
    number = np.shape(value)[0]
    if self.n == number:
      self._f = value
    else:
      msg = "Trying to set {0} forces for a system with {1} particles"
      raise ValueError(msg.format(number, self.n))

  @property
  def a(self):
    return self._f/self._mass[:, np.newaxis]

  @property
  def p(self):
    return self.v * self._mass[:, np.newaxis]

  @property
  def t(self):
    return self._t

  @t.setter
  def t(self, value):
    """
    Set types of particles.

    Parameters
    ----------

    value : 1D NumPy array or integer
        The new types of particles in an Nx3 array
    """
    if np.isscalar(value):
      self._t = np.array([value]*self.n, dtype=np.int32)
    else:
      value = np.array(value)
      number = np.shape(value)[0]
      if self.n == number:
        self._t = value
      else:
        msg = "Trying to set {0} types for a system with {1} particles"
        raise ValueError(msg.format(number, self.n))

  @property
  def mass(self):
    return self._mass

  @mass.setter
  def mass(self, value):
    """
    Set types of particles.

    Parameters
    ----------

    value : 1D NumPy array or integer
        The new types of particles in an Nx3 array
    """
    if np.isscalar(value):
      self._mass = np.array([value]*self.n, dtype=np.float32)
    else:
      value = np.array(value)
      number = np.shape(value)[0]
      if self.n == number:
        self._mass = value
      else:
        msg = "Trying to set {0} masses for a system with {1} particles"
        raise ValueError(msg.format(number, self.n))

class PointParticles(Base):
  """
  PointParticles class.
  Typical point particles used, for example, in LJ potential.
  """


class PointFermions(object):

  def __init__(self, n):
    self.n = n
    self._x = np.zeros((n, 3), dtype=np.float32)
    self._p = np.zeros((n, 3), dtype=np.float32)
    self._f = np.zeros((n, 3), dtype=np.float32)
    self._g = np.zeros((n, 3), dtype=np.float32)
    self._t = np.zeros(n, dtype=np.int32)
    self._mass = np.zeros(n, dtype=np.float32)
    self.idx = np.arange(n)

  @property
  def x(self):
    return self._x

  @x.setter
  def x(self, value):
    """
    Set positions of particles.

    Parameters
    ----------

    value : 2D NumPy array
        The new positons of particles in an Nx3 array
    """
    value = np.array(value, dtype=np.float32)
    number = np.shape(value)[0]
    if self.n == number:
      self._x = value
    else:
      msg = "Trying to set {0} positions for a system with {1} particles"
      raise ValueError(msg.format(number, self.n))

  @property
  def p(self):
    return self._v

  @p.setter
  def p(self, value):
    """
    Set momenta of particles.

    Parameters
    ----------

    value : 2D NumPy array
        The new velocities of particles in an Nx3 array
    """
    value = np.array(value, dtype=np.float32)
    number = np.shape(value)[0]
    if self.n == number:
      self._p = value
    else:
      msg = "Trying to set {0} momenta for a system with {1} particles"
      raise ValueError(msg.format(number, self.n))

  @property
  def f(self):
    return self._f

  @f.setter
  def f(self, value):
    """
    Set forces of particles.

    Parameters
    ----------

    value : 2D NumPy array
        The new forces of particles in an Nx3 array
    """
    value = np.array(value, dtype=np.float32)
    number = np.shape(value)[0]
    if self.n == number:
      self._f = value
    else:
      msg = "Trying to set {0} forces for a system with {1} particles"
      raise ValueError(msg.format(number, self.n))

  @property
  def g(self):
    return self._g

  @g.setter
  def g(self, value):
    """
    Set gorces of particles.

    Parameters
    ----------

    value : 2D NumPy array
        The new forces of particles in an Nx3 array
    """
    value = np.array(value, dtype=np.float32)
    number = np.shape(value)[0]
    if self.n == number:
      self._g = value
    else:
      msg = "Trying to set {0} gorces for a system with {1} particles"
      raise ValueError(msg.format(number, self.n))

  @property
  def v(self):
    return self._p/self._mass[:, np.newaxis] - self._g

  @property
  def t(self):
    return self._t

  @t.setter
  def t(self, value):
    """
    Set types of particles.

    Parameters
    ----------

    value : 1D NumPy array or integer
        The new types of particles in an Nx3 array
    """
    if np.isscalar(value):
      self._t = np.array([value]*self.n, dtype=np.int32)
    else:
      value = np.array(value)
      number = np.shape(value)[0]
      if self.n == number:
        self._t = value
      else:
        msg = "Trying to set {0} types for a system with {1} particles"
        raise ValueError(msg.format(number, self.n))

  @property
  def mass(self):
    return self._mass

  @mass.setter
  def mass(self, value):
    """
    Set types of particles.

    Parameters
    ----------

    value : 1D NumPy array or integer
        The new types of particles in an Nx3 array
    """
    if np.isscalar(value):
      self._mass = np.array([value]*self.n, dtype=np.float32)
    else:
      value = np.array(value)
      number = np.shape(value)[0]
      if self.n == number:
        self._mass = value
      else:
        msg = "Trying to set {0} masses for a system with {1} particles"
        raise ValueError(msg.format(number, self.n))

  @property
  def Tid(self):
    return self.kinetic_energy/(1.5*self.n)

  @property
  def kinetic_energy(self):
    return 0.5*np.sum(self.p**2/self._mass[:, np.newaxis])

  @property
  def Teff_corr(self):
    return -np.sum(self._p*self._g)/(3*self.n)

  def set_pos_box(self, L):
    Npart = self.n
    x = np.zeros((Npart, 3), dtype=np.float32)
    n3 = int(np.ceil(Npart**(1.0/3)))
    i = 0
    for p in it.product(range(n3),range(n3),range(n3)):
      if Npart <= i:
        break
        self._x[i, :] = np.array(p)*L/n3
        i += 1
