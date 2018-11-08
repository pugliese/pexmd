import numpy as np
import pexmd
import matplotlib.pylab as plt
import sys
import itertools as it
import time
from datetime import datetime

# Parametros
# Pauli: qo = 6.00fm , po = 2.067 MeV*1E-22 s/fm , D = 34.32MeV, h_bar = 6.582119 MeV*1E-22s
# Coulomb: k = e^2 = h_bar * c / 137 = 1.4403427984368629 Mev*fm
# Nuclear: Vo = 25.93 , r1 = 1.757 fm , p1 = 6.2 , r2 = 1.771 , p2 = 3 , d = 3.350fm , a = 5/6fm
# Masa proton: mp = 938.27203 MeV/c*c = 938.27203 MeV /(30)**2 = 1.0422

mp = 1.0422
scut = 2000
# Pauli
qo = 6.00
po = 2.067
h_barra = 6.582119
D = 34.32
pauli = pexmd.interaction.Pauli(scut, D, qo, po)

# Coulomb
e2 = 1.4403427984368629
coul = pexmd.interaction.Coulomb(scut, e2)

# Nuclear
Vo = 25.93
r1 = 1.757
p1 = 6.2
r2 = 1.771
p2 = 3
d = 3.350
a = 5.0/6.0
nuclear = pexmd.interaction.QCNM(scut, Vo, r1, p1, r2, p2, d, a)


def FixedPoint(parts, pairs_pauli, pairs_coul, pairs_nuc, dt, k=2, Niter = 5):
  Z_x = parts.x
  Z_p = parts.p
  for l in range(Niter):
    Z_x = parts.x + .5*dt*(Z_p/parts.mass[0] - parts.g)
    Z_p = parts.p + .5*dt*parts.f
    parts.f, parts.g, e_pauli, e_coul, e_nuc, ecaja = calc_fgorces(Z_x, Z_p, pairs_pauli, pairs_coul, pairs_nuc, k)
  parts.x = 2*Z_x - parts.x
  parts.p = 2*Z_p - parts.p
  parts.f, parts.g, e_pauli, e_coul, e_nuc, e_caja = calc_fgorces(parts.x, parts.p, pairs_pauli, pairs_coul, pairs_nuc)
  return e_pauli, e_coul, e_nuc, ecaja

def calc_fgorces(X, P, pairs_pauli, pairs_coul, pairs_nuc, k=2):
  trampa = lambda x: (-k*x, 0.5*k*np.sum(x**2))
  f_pau, g, e_pauli = pauli.fgorces(X, P, pairs_pauli)
  f_coul, e_coul = coul.forces(X, P, pairs_coul)
  f_nuc, e_nuc = nuclear.forces(X, P, pairs_nuc)
  fcaja, ecaja = trampa(X)
  f_tot = f_pau+f_coul+f_nuc+fcaja
  return f_tot, g, e_pauli, e_coul, e_nuc, ecaja

def armar_atomo(N_prot, N_neut, To=0.2, L=4):
  Npart = N_prot + N_neut
  parts = pexmd.particles.PointFermions(Npart)
  parts.x = parts.set_pos_box(L) - 0.5*L
  #parts.p = np.random.normal(0, np.sqrt(To*1.5*mp), (Npart, 3))
  parts.mass = mp
  parts.t = [0 for i in range(N_prot)] + [1 for i in range(N_neut)]
  pairs_nuc = np.array(list(it.combinations(range(Npart), 2)), dtype=np.int64)
  '''
  pairs_pauli = list(it.combinations(range(N_prot//2), 2))
  pairs_pauli += list(it.combinations(range(N_prot//2, N_prot), 2))
  pairs_pauli += list(it.combinations(range(N_prot, N_prot+N_neut//2), 2))
  pairs_pauli += list(it.combinations(range(N_prot+N_neut//2, Npart), 2))
  '''
  pairs_pauli = list(it.combinations(range(N_prot), 2))
  pairs_pauli += list(it.combinations(range(N_prot, Npart), 2))
  pairs_pauli = np.array(pairs_pauli, dtype=np.int64)
  pairs_coul = np.array(list(it.combinations(range(N_prot), 2)), dtype=np.int64)
  return parts, pairs_pauli, pairs_nuc, pairs_coul


def apagar_trampa(Nstep, parts, pairs_pauli, pairs_coul, pairs_nuc, dt, k=2, Niter = 5):
  e_pauli = np.zeros(Nstep)
  e_coul = np.zeros(Nstep)
  e_nuc = np.zeros(Nstep)
  ecaja = np.zeros(Nstep)
  ecin = np.zeros(Nstep)
  for i in range(Nstep):
    e_pauli[i], e_coul[i], e_nuc[i], ecaja[i] = FixedPoint(parts, pairs_pauli, pairs_coul, pairs_nuc, h, k*(1.0 - 1.0*i/Nstep))
    ecin[i] = parts.kinetic_energy
  return e_pauli, e_coul, e_nuc, ecaja, ecin


def dist_fases(parts):
  q = []
  p = []
  for i in range(1, parts.n):
    for j in range(i):
      q.append(np.linalg.norm(parts.x[i] - parts.x[j]))
      p.append(np.linalg.norm(parts.p[i] - parts.p[j]))
  return np.array(q), np.array(p)

def energy(parts):
  parts.f, parts.g, e_pauli, e_coul, e_nuc, e_caja = calc_fgorces(parts.x, parts.p, pairs_pauli, pairs_coul, pairs_nuc)
  return parts.kinetic_energy + e_pauli + e_coul + e_nuc

def save_checkpoint(parts, filename='checkpoint_nucleo.txt'):
  data = []
  for i in range(parts.n):
    data.append([])
    data[i].append(parts.t[i])
    data[i].append(parts.mass[i])
    data[i] = data[i] + list(parts.x[i, :])
    data[i] = data[i] + list(parts.p[i, :])
    data[i] = data[i] + list(parts.f[i, :])
    data[i] = data[i] + list(parts.g[i, :])
  data = np.array(data)
  np.savetxt(filename, data, fmt="%1d %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e")


def load_checkpoint(filename='checkpoint_nucleo.txt'):
  data = np.loadtxt(filename, unpack=True)
  parts = pexmd.particles.PointFermions(len(data[0]))
  parts.t = data[0]
  parts.mass = data[1]
  parts.x = np.transpose(data[2:5])
  parts.p = np.transpose(data[5:8])
  parts.f = np.transpose(data[8:11])
  parts.g = np.transpose(data[11:14])
  return parts

def pal_vmd(filename):
  new_filename = filename[:-3] + "lammpstrj"
  parts = load_checkpoint(filename)
  f = open(filename, "r")
  data = f.read()
  f.close()
  f = open(new_filename, "w")
  lims = [min(parts.x[:, 0]), max(parts.x[:, 0]), min(parts.x[:, 1]), max(parts.x[:,1]), min(parts.x[:,2]), max(parts.x[:, 2])]
  header1 = "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n{0}\n".format(parts.n)
  header2 = "ITEM: BOX BOUNDS pp pp pp\n{0} {1}\n{2} {3}\n{4} {5}\n".format(lims[0], lims[1], lims[2], lims[3], lims[4], lims[5])
  header3 = "ITEM: ATOMS type mass x y z vx vy vz fx fy fz gx gy gz\n"
  f = open(new_filename, "w")
  f.write(header1 + header2 + header3 + data)
  f.close()


def max_dist(x):
  return np.sqrt(max(sum(x**2)))

def muestrear_N(parts, N=int(1E4)):
  E_real = np.zeros(N)
  T_eff = np.zeros(N)
  for i in range(N):
    e_pauli, e_coul, e_nuc, e_caja = FixedPoint(parts, pairs_pauli, pairs_coul, pairs_nuc, h, 0)
    e_cin = parts.kinetic_energy
    E_real[i] = e_coul + e_nuc + e_cin + e_pauli
    T_eff[i] = e_cin/(1.5*parts.n) + parts.Teff_corr
  return E_real, T_eff

def enfriar(parts, N_resc, N_steps, N_ramp, k, factor = 0.9):
  N_tot = 2*N_ramp + N_resc*N_steps
  E_real = np.zeros(N_tot)
  E_tot = np.zeros(N_tot)
  E_cin = np.zeros(N_tot)
  idx = 0
  for i in range(N_ramp):
    e_pauli, e_coul, e_nuc, e_caja = FixedPoint(parts, pairs_pauli, pairs_coul, pairs_nuc, h, k*(1.0+i)/N_ramp)
    E_cin[idx] = parts.kinetic_energy
    E_real[idx] = e_coul + e_nuc + E_cin[idx] + e_pauli
    E_tot[idx] = E_real[idx] + e_caja
    idx += 1
  for j in range(N_resc):
    parts.p = factor*parts.p
    for i in range(N_steps):
      e_pauli, e_coul, e_nuc, e_caja = FixedPoint(parts, pairs_pauli, pairs_coul, pairs_nuc, h, k)
      E_cin[idx] = parts.kinetic_energy
      E_real[idx] = e_coul + e_nuc + E_cin[idx] + e_pauli
      E_tot[idx] = E_real[idx] + e_caja
      idx += 1
  for i in range(N_ramp):
    e_pauli, e_coul, e_nuc, e_caja = FixedPoint(parts, pairs_pauli, pairs_coul, pairs_nuc, h, k - k*(1.0+i)/N_ramp)
    E_cin[idx] = parts.kinetic_energy
    E_real[idx] = e_coul + e_nuc + E_cin[idx] + e_pauli
    E_tot[idx] = E_real[idx] + e_caja
    idx += 1
  if(E_real[0]<E_real[-1]):
    print('La energia aumenta')
  return E_real, E_tot, E_cin

# Atomo
parts, pairs_pauli, pairs_nuc, pairs_coul = armar_atomo(20, 20)
h = 1E-3
