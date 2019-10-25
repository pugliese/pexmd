import numpy as np
import ctypes as ct

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

interaction_QCNM_c = C_funcs.interaction_QCNM
interaction_QCNM_c.argtypes = [ct.c_float, ct.c_voidp]
interaction_QCNM_c.restype = ct.c_float

class QCNM(ct.Structure):
  """
  QCNM potential for nucleon-nucleon interaction
  """
  _fields_ = [("V_o", ct.c_float), ("p_1", ct.c_float), ("p_2", ct.c_float), ("r_1", ct.c_float), ("r_2", ct.c_float), ("d", ct.c_float), ("a", ct.c_float),
                  ("rcut", ct.c_float), ("rcut2", ct.c_float), ("shift", ct.c_float), ("LUT", ct.c_voidp), ("dr2", ct.c_float)]


  def __init__(self, V_o, p_1, p_2, r_1, r_2, d, a, rcut):
    self.V_o = V_o
    self.p_1 = p_1
    self.p_2 = p_2
    self.r_1 = r_1
    self.r_2 = r_2
    self.d = d
    self.a = a
    self.rcut = rcut
    self.rcut2 = rcut**2
    self.shift = V_o*((r_1/rcut)**p_1 - (r_2/rcut)**p_2)/(1 + np.exp((rcut - d)/a))

  def actualizar_shift(self):
    self.shift = self.V_o*((self.r_1/self.rcut)**self.p_1 - (self.r_2/self.rcut)**self.p_2)/(1 + np.exp((self.rcut - self.d)/self.a))

  def potential(self, r):
    #return self.V_o*((self.r_1/r)**self.p_1 - (self.r_2/r)**self.p_2)/(1 + np.exp((r - self.d)/self.a)) - self.shift
    return interaction_QCNM_c(r, ct.pointer(self))

class TotalPotential(ct.Structure):
  """
  Total potential of interaction to match struct from C
  """
  _fields_ = [("pauli", ct.POINTER(Pauli)), ("panda_nn", ct.POINTER(Panda_nn)), ("panda_np", ct.POINTER(Panda_np)), ("qcnm", ct.POINTER(QCNM)),  ("rcut", ct.c_float)]

  def __init__(self, pauli, panda_nn, panda_np, qcnm):
    self.pauli = ct.pointer(pauli)
    self.panda_nn = ct.pointer(panda_nn)
    self.panda_np = ct.pointer(panda_np)
    self.qcnm = ct.pointer(qcnm)
    self.rcut = max([np.sqrt(pauli.scut2)*pauli.qo, panda_nn.rcut, panda_np.rcut, qcnm.rcut])

load_lammpstrj_c = C_funcs.load_lammpstrj
load_lammpstrj_c.argtypes = [ct.c_char_p, ct.c_voidp, ct.c_voidp, ct.c_float]
load_lammpstrj_c.restype = ct.c_int

energia_sin_LUT_c = C_funcs.energia_sin_LUT
energia_sin_LUT_c.argtypes = [ct.c_voidp, ct.c_voidp]
energia_sin_LUT_c.restype = ct.c_int

energia_con_QCNM_c = C_funcs.energia_QCNM
energia_con_QCNM_c.argtypes = [ct.c_voidp, ct.c_voidp]
energia_con_QCNM_c.restype = ct.c_int

liberar_c = C_funcs.liberar
liberar_c.argtypes = [ct.c_voidp]
liberar_c.restype = ct.c_int

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

  def __del__(self):
    liberar_c(ct.pointer(self))

  def load_lammpstrj(self, filename, L, rcut):
    filename_p = ct.c_char_p(filename.encode('utf-8'))
    #filename_p = ct.pointer(filename.encode('utf-8'))
    L_p = L.ctypes.data_as(ct.c_voidp)
    rcut_32 = np.float32(rcut)
    return load_lammpstrj_c(filename_p, ct.pointer(self), L_p, rcut_32)

  def energy(self, pot_tot):
    return energia_sin_LUT_c(ct.pointer(self), ct.pointer(pot_tot))

  def energy_QCNM(self, pot_tot):
    return energia_con_QCNM_c(ct.pointer(self), ct.pointer(pot_tot))


if (True):
  h_bar = 197.327
  m = 938
  deg = lambda E: 4*np.pi*np.sqrt(2*m**3*E)*(L/h)**3
  long_term = lambda T: (2*np.pi*h_bar**2/(m*T))**0.5
  Ef_rho = h_bar**2/(2*m)*(6*np.pi**2)**(2/3)

  data = np.loadtxt("../FD_fit/LUT_F32.txt")
  z32 = data[0,:]
  F32 = data[1,:]
  F32_A = lambda z: (4/(3*np.pi**0.5))*(np.log(z)**1.5)*(1 + np.pi**2/(8*np.log(z)**2) + (7*np.pi**4/640)*np.log(z)**(-4))

  data = np.loadtxt("../FD_fit/LUT_F52.txt")
  z52 = data[0,:]
  F52 = data[1,:]

  F52_A = lambda z: (8/(15*np.pi**0.5))*(np.log(z)**2.5)*(1 + 5*np.pi**2/(8*np.log(z)**2) - (5*3/16)*(7*np.pi**4/360)*np.log(z)**-4)

  def dame_z(Y):
    F32_A = lambda z: (4/(3*np.pi**0.5))*(np.log(z)**1.5)*(1 + np.pi**2/(8*np.log(z)**2) + (7*np.pi**4/640)*np.log(z)**(-4))
    if (30<=Y):
      z_inf = 30.0
      z_sup = 2*z_inf
      while(F32_A(z_sup) < Y):
        z_sup *= 2
      z_med = (z_inf+z_sup)/2
      if(z_sup==np.inf):
        return np.exp((3*np.pi**0.5*Y/4)**(2/3))
      while (Y*1E-4 < np.abs(F32_A(z_med)-Y)):
        z_med = (z_inf+z_sup)/2
        if(F32_A(z_med)<Y):
          z_inf = z_med
        else:
          z_sup = z_med
      return z_med
    else:
      inf = 0
      sup = len(z32)-1
      while (inf<sup-1):
        med = (inf+sup)//2
        if(F32[med]<Y):
          inf = med
        else:
          sup = med
      m = (F32[inf+1]-F32[inf])/(z32[inf+1]-z32[inf])
      return z32[inf]+(Y-F32[inf])/m

  def F_32(fug):
    if (z32[-1]<=fug):
      return F32_A(fug)
    else:
      inf = 0
      sup = len(z32)-1
      while (inf<sup-1):
        med = (inf+sup)//2
        if(z32[med]<fug):
          inf = med
        else:
          sup = med
      m = (F32[inf+1]-F32[inf])/(z32[inf+1]-z32[inf])
      return F32[inf] + m*(fug-z32[inf])

  def F_52(fug):
    if (z52[-1]<=fug):
      return F52_A(fug)
    else:
      inf = 0
      sup = len(z52)-1
      while (inf<sup-1):
        med = (inf+sup)//2
        if(z52[med]<fug):
          inf = med
        else:
          sup = med
      m = (F52[inf+1]-F52[inf])/(z52[inf+1]-z52[inf])
      return F52[inf] + m*(fug-z52[inf])

  def P_rhoT(A):
    return F_52(dame_z(A))/A

  def kinetic_energy(rho, T):
    n = len(rho)
    m = len(T)
    res = np.zeros((n, m))
    for r in range(n):
      for t in range(m):
        A = long_term(T[t])**3*rho[r]
        res[r,t] = 1.5*P_rhoT(A)*T[t]
    return res
