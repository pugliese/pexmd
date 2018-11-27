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
# Masa proton: mp = 938.27203 MeV/c*c = 938.27203 MeV /(29.98fm/s * 1E22)**2 = 1.043916

mp = 1.043916
scut = 2000
# Pauli
qo = 6.00
po = 2.067
h_barra = 6.582119
D = 34.32
pauli = pexmd.interaction.Pauli(scut, D*(h_barra/(qo*po))**3, qo, po)
#pauli = pexmd.interaction.Pauli(scut, 6.7934, qo, po)
#pauli = pexmd.interaction.Pauli(scut, D, qo, po)

# Coulomb
e2 = 1.439965 # Wolfram
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


# ------------------- GENERALES ----------------------- #

def FixedPoint(parts, pairs_pauli, pairs_coul, pairs_nuc, dt, Niter = 5):
  Z_x = parts.x
  Z_p = parts.p
  for l in range(Niter):
    Z_x = parts.x + .5*dt*(Z_p/parts.mass[0] - parts.g)
    Z_p = parts.p + .5*dt*parts.f
    parts.f, parts.g, e_pauli, e_coul, e_nuc = calc_fgorces(Z_x, Z_p, pairs_pauli, pairs_coul, pairs_nuc)
  parts.x = 2*Z_x - parts.x
  parts.p = 2*Z_p - parts.p
  parts.f, parts.g, e_pauli, e_coul, e_nuc = calc_fgorces(parts.x, parts.p, pairs_pauli, pairs_coul, pairs_nuc)
  return e_pauli, e_coul, e_nuc

def calc_fgorces(X, P, pairs_pauli, pairs_coul, pairs_nuc):
  f_pau, g, e_pauli = pauli.fgorces(X, P, pairs_pauli)
  f_coul, e_coul = coul.forces(X, P, pairs_coul)
  f_nuc, e_nuc = nuclear.forces(X, P, pairs_nuc)
  f_tot = f_pau+f_coul+f_nuc
  return f_tot, g, e_pauli, e_coul, e_nuc


def pares_interaccion(t):
  protones = np.where(t==0)[0]
  neutrones = np.where(t==1)[0]
  pairs_pauli = list(it.combinations(protones, 2)) + list(it.combinations(neutrones, 2))
  pairs_nuc = list(it.combinations(range(len(t)), 2))
  pairs_coul = list(it.combinations(protones, 2))
  return np.array(pairs_pauli), np.array(pairs_nuc), np.array(pairs_coul)

def unir_particulas(parts1, parts2):
  parts = pexmd.particles.PointFermions(parts1.n + parts2.n)
  parts.x = np.concatenate([parts1.x, parts2.x])
  parts.p = np.concatenate([parts1.p, parts2.p])
  parts.t = np.concatenate([parts1.t, parts2.t])
  parts.mass = np.concatenate([parts1.mass, parts2.mass])
  return parts


# filename: configuracion inicial de cada nucleo
# d: distancia inicial (idealmente d-2r>scut o algo asi)
# p: impulso inicial neto
def colision_frontal(filename, d, p, eje=0):
  parts1 = load_checkpoint(filename)
  parts2 = load_checkpoint(filename)
  parts1.x[:, eje] -= d/2
  parts1.p[:, eje] -= -p/2
  parts2.x[:, eje] += d/2
  parts2.p[:, eje] += -p/2
  parts = unir_particulas(parts1, parts2)
  pairs_pauli, pairs_nuc, pairs_coul = pares_interaccion(parts.t)
  parts.f, parts.g, e_pauli, e_coul, e_nuc = calc_fgorces(parts.x, parts.p, pairs_pauli, pairs_coul, pairs_nuc)
  return parts, pairs_pauli, pairs_nuc, pairs_coul

def dist_fases(parts):
  q = []
  p = []
  for i in range(1, parts.n):
    for j in range(i):
      q.append(np.linalg.norm(parts.x[i] - parts.x[j]))
      p.append(np.linalg.norm(parts.p[i] - parts.p[j]))
  return np.array(q), np.array(p)

def energy(parts):
  parts.f, parts.g, e_pauli, e_coul, e_nuc  = calc_fgorces(parts.x, parts.p, pairs_pauli, pairs_coul, pairs_nuc)
  return parts.kinetic_energy + e_pauli + e_coul + e_nuc

def energies(parts):
  parts.f, parts.g, e_pauli, e_coul, e_nuc  = calc_fgorces(parts.x, parts.p, pairs_pauli, pairs_coul, pairs_nuc)
  return parts.kinetic_energy, e_pauli, e_coul, e_nuc


# ---------------- CHECKPOINTS ----------------------- #

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


# ---------------- VISUALIZACION ----------------------- #

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

def video(filename, parts, N_frames, N_skip = 10):
  f = open(filename+".lammpstrj", 'w')
  formato = "%1d %.18e %.18e %.18e\n"
  for i in range(N_frames):
    if (i%(N_frames/10) == 0):
      print('{0}%'.format(10*i/N_frames))
    for j in range(N_skip):
      e_pauli, e_coul, e_nuc = FixedPoint(parts, pairs_pauli, pairs_coul, pairs_nuc, h)
    lims = [min(parts.x[:, 0]), max(parts.x[:, 0]), min(parts.x[:, 1]), max(parts.x[:,1]), min(parts.x[:,2]), max(parts.x[:, 2])]
    header1 = "ITEM: TIMESTEP\n{0}\nITEM: NUMBER OF ATOMS\n{1}\n".format(i*N_skip, parts.n)
    header2 = "ITEM: BOX BOUNDS pp pp pp\n{0} {1}\n{2} {3}\n{4} {5}\n".format(lims[0], lims[1], lims[2], lims[3], lims[4], lims[5])
    header3 = "ITEM: ATOMS type x y z\n"
    f.write(header1 + header2 + header3)
    for j in range(parts.n):
      f.write(formato %(parts.t[j], parts.x[j, 0], parts.x[j, 1], parts.x[j, 2]))
  f.close()


h = 1E-3
d = 30
p = 10

if (len(sys.argv)>=2):
  d = float(sys.argv[1])
  if (len(sys.argv)>=3):
    p = float(sys.argv[2])
parts, pairs_pauli, pairs_nuc, pairs_coul = colision_frontal('checkpoint_Ca.txt', d, p)
