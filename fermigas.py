import numpy as np
import pexmd
import matplotlib.pylab as plt
import sys
import itertools as it
import time
from datetime import datetime

def FixedPoint(parts, pairs, pauli, dt, Niter = 5):
  Z_x = parts.x
  Z_p = parts.p
  #fuerzas_caja = lambda x: -5*(x/3)**5
  for k in range(Niter):
    parts.f, parts.g, e = pauli.fgorces(Z_x, Z_p)
    fcaja = 0#fuerzas_caja(parts.x)
    Z_x = parts.x + .5*dt*parts.v
    Z_p = parts.p + .5*dt*(parts.f+fcaja)
  return 2*Z_x - parts.x, 2*Z_p - parts.p

def avanzar_fp(parts, pairs, pauli, dt):
  parts.x, parts.p = FixedPoint(parts, pairs, pauli, dt)
  parts.f, parts.g, e = pauli.fgorces(parts.x, parts.p, pairs)
  parts.x, parts.p = bx.wrap_boundary(parts.x, parts.p)

def Dist_Fases(parts, pairs):
  qs = np.zeros(len(pairs))
  ps = np.zeros(len(pairs))
  for k in range(len(pairs)):
    i = pairs[k, 0]
    j = pairs[k, 1]
    qs[k] = np.linalg.norm(parts.x[i] - parts.x[j])
    ps[k] = np.linalg.norm(parts.p[i] - parts.p[j])
  return qs, ps


# Interaction
qo = 1.664
po = 120
h_barra = 196.727394
D = 207
scut = 2000
pauli = pexmd.interaction.Pauli(scut, D, qo, po)
# Box
L = 2#5*qo
bx = pexmd.box.Box(-L*np.ones(3)/2, L*np.ones(3)/2, 'Fixed')
# Thermostat
Q = 0.02
To = [1]#np.concatenate((np.linspace(8, 2, 4), np.linspace(1.0, 0.2, 5)))
therm = pexmd.thermostat.Andersen(To[0], 0)
# Particles
Npart = 50#30
m = 0.5109989461
parts = pexmd.particles.PointFermions(Npart)
parts.set_pos_box(L)
parts.x = parts.x - 0.5*L
parts.p = np.random.normal(0, np.sqrt(To[0]*1.5*m), (Npart, 3))
parts.mass = m
pairs = np.array(list(it.combinations(range(Npart), 2)), dtype=np.int64)
# Integrator
h = 0.0005
integ = pexmd.integrator.VelVerlet(h)
# Potencial de caja

def tiempo(Nparts):
  tiempos = np.zeros(len(Nparts))
  i = 0
  for n in Nparts:
    print(n)
    m = 0.5109989461
    parts = pexmd.particles.PointFermions(n)
    parts.set_pos_box(L)
    parts.x = parts.x - 0.5*L
    parts.p = np.random.normal(0, np.sqrt(1.5*m), (n, 3))
    parts.mass = m
    pairs = np.array(list(it.combinations(range(n), 2)), dtype=np.int64)
    t0 = datetime.now()
    for j in range(1000):
      avanzar_fp(parts, pairs, pauli, h)
    t1 = datetime.now()
    tiempos[i]= (t1-t0).total_seconds()
    i += 1
  return tiempos


"""
Nterm = 50000
Ntemp = len(To)
Epot = np.zeros(Ntemp*Nterm)
Ecin = np.zeros(Ntemp*Nterm)
Tcorr = np.zeros(Ntemp*Nterm)
rmax = np.zeros(Ntemp*Nterm)
t = time.time()
i = 0
print("Termalizacion inicial")
for k in range(Nterm):
  avanzar_fp(parts, pairs, pauli, h)
for j in range(Ntemp):
  therm = pexmd.thermostat.Andersen(To[j], 0)
  print("Temperatura:", To[j])
  for k in range(Nterm):
    avanzar_fp(parts, pairs, pauli, h)
    Ecin[i] = parts.kinetic_energy
    i += 1
  #parts.v = parts.v*np.sqrt(1.5*To[j]*Npart/Ecin[i-1])
t = time.time() - t
T = Ecin/(1.5*Npart)
Teff = T + Tcorr
Etot = Ecin + Epot
E = [np.mean(Etot[i*Nterm:(i+1)*Nterm]) for i in range(Ntemp)]
Estd = [np.std(Etot[i*Nterm:(i+1)*Nterm]) for i in range(Ntemp)]
Teff_mean = [np.mean(Teff[i*Nterm:(i+1)*Nterm]) for i in range(Ntemp)]
Teff_std = [np.std(Teff[i*Nterm:(i+1)*Nterm]) for i in range(Ntemp)]
print("%d Pasos en %d min, %d segs" %(Nterm*(Ntemp+1), t/60, t%60))

qs, ps = Dist_Fases(parts, pairs)
plt.plot(qs, ps, "b.")
plt.xlabel("q")
plt.ylabel("p")
plt.axis([0, np.max(qs), 0, np.max(ps)])
plt.figure()
plt.plot(Etot, "b-")
plt.show()
"""
