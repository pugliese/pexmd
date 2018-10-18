import numpy as np
import pexmd
import matplotlib.pylab as plt
import sys
import itertools as it
import time

def FixedPoint(parts, pairs, pauli, dt, Niter = 3):
  Z_x = parts.x
  Z_p = parts.vp
  for k in range(Niter):
    parts.f, parts.g, e = pauli.fgorces(Z_x, Z_p)
    Z_x = parts.x + .5*dt*parts.v
    Z_p = parts.p + .5*dt*parts.f
  return 2*Z_x - parts.x, 2*Z_p - parts.p

def avanzar_fp(parts, pairs, pauli, integ, caja, therm):
  parts.x, parts.v = FixedPoint(parts, pairs, pauli, integ.dt, caja, 5)
  parts.f, parts.g, e = pauli.fgorces(parts.x, parts.p, pairs)
  parts.x, parts.v = bx.wrap_boundary(parts.x, parts.v)

# Interaction
qo = 1.664
po = 120
h_barra = 196.727394
D = 207
scut = 2000
pauli = pexmd.interaction.Pauli(scut, D, qo, po)
# Box
L = 2#5*qo
bx = pexmd.box.Box(np.zeros(3), L*np.ones(3), 'Fixed')
# Thermostat
Q = 0.02
To = np.concatenate((np.linspace(8, 2, 4), np.linspace(1.0, 0.2, 5)))
therm = pexmd.thermostat.Andersen(To[0], 0)
# Particles
Npart = 50#30
m = 0.5109989461
parts = pexmd.particles.PointParticles(Npart)
parts.set_pos_box(L)
parts.p = np.random.normal(0, np.sqrt(To[0]*1.5/m), (Npart, 3))
parts.mass = m
pairs = np.array(list(it.combinations(range(Npart), 2)), dtype=np.int64)
# Integrator
h = 0.0005
integ = pexmd.integrator.VelVerlet(h)
# Potencial de caja

Nterm = 5000
Ntemp = len(To)
Epot = np.zeros(Ntemp*Nterm)
Ecin = np.zeros(Ntemp*Nterm)
Tcorr = np.zeros(Ntemp*Nterm)
rmax = np.zeros(Ntemp*Nterm)
t = time.time()
i = 0
print("Termalizacion inicial")
for k in range(Nterm):
  parts.x, parts.v, a, b, c, d = avanzar_fp(parts, pairs, pauli, integ, bx, therm)
for j in range(Ntemp):
  therm = pexmd.thermostat.Andersen(To[j], 0)
  print("Temperatura:", To[j])
  for k in range(Nterm):
    parts.x, parts.v, Epot[i], Epotcaja[i], rmax[i], Tcorr[i] = avanzar_fp(parts, pairs, pauli, integ, bx, therm)
    Ecin[i] = energia_cinetica(parts.v, parts.mass)
    i += 1
  parts.v = parts.v*np.sqrt(1.5*To[j]*Npart/Ecin[i-1])
t = time.time() - t
T = Ecin/(1.5*Npart)
Teff = T + Tcorr
Etot = Ecin + Epot
E = [np.mean(Etot[i*Nterm:(i+1)*Nterm]) for i in range(Ntemp)]
Estd = [np.std(Etot[i*Nterm:(i+1)*Nterm]) for i in range(Ntemp)]
Teff_mean = [np.mean(Teff[i*Nterm:(i+1)*Nterm]) for i in range(Ntemp)]
Teff_std = [np.std(Teff[i*Nterm:(i+1)*Nterm]) for i in range(Ntemp)]
print("%d min, %d segs" %(t/60, t%60))
