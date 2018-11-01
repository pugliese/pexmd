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
e2 = 0#1.4403427984368629
coul = pexmd.interaction.Coulomb(scut, e2)

# Nuclear
Vo = 0#25.93
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
  pairs_pauli = np.array(list(it.combinations(range(N_prot), 2))+list(it.combinations(range(N_prot, Npart), 2)), dtype=np.int64)
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

def save_checkpoint(parts, filename='checkpoint_nucleo.txt'):
  f = open(filename, 'w')
  f.write(str(parts.n)+'\n')
  for i in range(parts.n):
    f.write(str(parts.x[i])[1:-1]+' \n')
  for i in range(parts.n):
    f.write(str(parts.p[i])[1:-1]+' \n')
    #f.write('{0}, {1}, {2}'.format(parts.p[i]))
  for i in range(parts.n):
    f.write(str(parts.f[i])[1:-1]+' \n')
    #f.write('{0}, {1}, {2}'.format(parts.f[i]))
  for i in range(parts.n):
    f.write(str(parts.g[i])[1:-1]+' \n')
    #f.write('{0}, {1}, {2}'.format(parts.g[i]))
  for i in range(parts.n):
    f.write('{0} {1} \n'.format(parts.mass[i], parts.t[i]))
  f.close()

def load_checkpoint(filename='checkpoint_nucleo.txt'):
  f = open(filename, 'r')
  line = f.readline()
  n = int(line)
  parts = pexmd.particles.PointFermions(n)
  for i in range(n):
    line = f.readline()
    data = line.split(' ')
    print(data)
    parts.x[i, 0] = float(data[0])
    parts.x[i, 1] = float(data[1])
    parts.x[i, 2] = float(data[2])
  for i in range(parts.n):
    line = f.readline()
    data = line.split(' ')
    parts.p[i, 0] = float(data[0])
    parts.p[i, 1] = float(data[1])
    parts.p[i, 2] = float(data[2])
  for i in range(parts.n):
    line = f.readline()
    data = line.split(' ')
    parts.f[i, 0] = float(data[0])
    parts.f[i, 1] = float(data[1])
    parts.f[i, 2] = float(data[2])
  for i in range(parts.n):
    line = f.readline()
    data = line.split(' ')
    parts.g[i, 0] = float(data[0])
    parts.g[i, 1] = float(data[1])
    parts.g[i, 2] = float(data[2])
  for i in range(parts.n):
    line = f.readline()
    data = line.split(' ')
    parts.t[i] = int(data[0])
    parts.mass[i] = float(data[1])
  f.close()
  return parts

def max_dist(x):
  return np.sqrt(max(sum(x**2)))


parts, pairs_pauli, pairs_nuc, pairs_coul = armar_atomo(3, 3)

h = 1E-3
N = int(1/h)
N_resc = 50
fact_resc = 0.8
e_pauli = np.zeros((N_resc+10)*N)
e_coul = np.zeros((N_resc+10)*N)
e_nuc = np.zeros((N_resc+10)*N)
ecaja = np.zeros((N_resc+10)*N)
ecin = np.zeros((N_resc+10)*N)
tam = np.zeros((N_resc+10)*N)
dists = np.zeros(((N_resc+10)*N, 6))
qs = np.zeros((N_resc+4, 15))
ps = np.zeros((N_resc+4, 15))
k_init = 15
parts.f, parts.g, e_pauli[0], e_coul[0], e_nuc[0], ecaja[0] = calc_fgorces(parts.x, parts.p, pairs_pauli, pairs_coul, pairs_nuc)
for i in range(N):
  e_pauli[i], e_coul[i], e_nuc[i], ecaja[i] = FixedPoint(parts, pairs_pauli, pairs_coul, pairs_nuc, h, k_init)
  ecin[i] = parts.kinetic_energy
  tam[i] = max_dist(parts.x)
  dists[i] = np.sqrt(sum(np.transpose(parts.x)**2))
qs[0, :], ps[0, :] = dist_fases(parts)
for j in range(1, N_resc+1):
    parts.p = fact_resc*parts.p
    for i in range(N):
      e_pauli[i+j*N], e_coul[i+j*N], e_nuc[i+j*N], ecaja[i+j*N] = FixedPoint(parts, pairs_pauli, pairs_coul, pairs_nuc, h, k_init)
      ecin[i+j*N] = parts.kinetic_energy
      tam[i+j*N] = max_dist(parts.x)
      dists[i+j*N] = np.sqrt(sum(np.transpose(parts.x)**2))
    qs[j, :], ps[j, :] = dist_fases(parts)
for i in range(3*N):
  e_pauli[i+(N_resc+1)*N], e_coul[i+(N_resc+1)*N], e_nuc[i+(N_resc+1)*N], ecaja[i+(N_resc+1)*N] = FixedPoint(parts, pairs_pauli, pairs_coul, pairs_nuc, h, k_init*(3.0-1.0*i/N)/3)
  ecin[i+(N_resc+1)*N] = parts.kinetic_energy
  tam[i+(N_resc+1)*N] = max_dist(parts.x)
  dists[i+(N_resc+1)*N] = np.sqrt(sum(np.transpose(parts.x)**2))
qs[N_resc+2, :], ps[N_resc+2, :] = dist_fases(parts)
for i in range(6*N):
  e_pauli[i+(N_resc+4)*N], e_coul[i+(N_resc+4)*N], e_nuc[i+(N_resc+4)*N], ecaja[i+(N_resc+4)*N] = FixedPoint(parts, pairs_pauli, pairs_coul, pairs_nuc, h, 0)
  ecin[i+(N_resc+4)*N] = parts.kinetic_energy
  tam[i+(N_resc+4)*N] = max_dist(parts.x)
  dists[i+(N_resc+4)*N] = np.sqrt(sum(np.transpose(parts.x)**2))
qs[N_resc+3, :], ps[N_resc+3, :] = dist_fases(parts)

plt.figure()
plt.plot(e_pauli, 'k')
plt.plot(e_coul, 'b')
plt.plot(e_nuc, 'r')
plt.plot(ecaja, 'y')
plt.plot(ecin, 'r--')
plt.plot(ecin+e_pauli + e_coul + ecaja + e_nuc, 'c--')
plt.legend(['Pauli', 'Coulomb', 'Nuclear', 'Trampa', 'Cinetica', 'Total'])
plt.figure()
for i in range(N_resc+4):
  plt.plot(qs[i], ps[i], '.')
plt.xlabel('q')
plt.ylabel('p')

'''
plt.figure()
plt.plot(tam, 'k')
plt.plot(dists[:, 0], 'k')
plt.plot(dists[:, 1], 'k')
plt.plot(dists[:, 2], 'k')
plt.plot(dists[:, 3], 'k')
plt.plot(dists[:, 4], 'k')
plt.plot(dists[:, 5], 'k')
'''

plt.show()
