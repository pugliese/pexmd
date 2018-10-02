import numpy as np
import pexmd
import matplotlib.pylab as plt
import sys
import itertools as it
import time

def FixedPoint(parts, interact, dt, caja, Niter = 3):
  Z_x = parts.x
  Z_v = parts.v
  for k in range(Niter):
    forces, e = pauli.forces(Z_x, parts.mass[0]*Z_v)
    gorces, e = pauli.gorces(Z_x, parts.mass[0]*Z_v)
    fuerza_aux, ecaja = caja(Z_x)
    Z_x = parts.x + .5*dt*(Z_v - gorces)
    Z_v = parts.v + .5*dt*(forces+fuerza_aux)/parts.mass[0]
  return 2*Z_x - parts.x, 2*Z_v - parts.v

def particulas(Npart,L):
  x = np.zeros((Npart,3), dtype=np.float32)
  n3 = int(np.ceil(Npart**(1.0/3)))
  i = 0
  for p in it.product(range(n3),range(n3),range(n3)):
    if Npart <= i:
      break
    x[i, :] = np.array(p)*L/n3
    i += 1
  return x

def energia_cinetica(v, m):
  return 0.5*sum(sum(v**2))*m[0]

def Corr_temp_eff(p, gorces):
  Npart = len(p[:,0])
  return -sum(sum(p*gorces))/(3*Npart)

def fuerzas_caja(x, L, V, k = 20):
  fuerzas = V*(0.5*L - x)**(k-1)
  epot = sum(sum((0.5*L - x)**k))*V/k
  return fuerzas, epot

def avanzar_fp(parts, pauli, integ, caja, therm):
  parts.x, parts.v = FixedPoint(parts, pauli, integ.dt, caja, 5)

  parts.f, ecaja = caja(parts.x)
  gorces, e = pauli.gorces(parts.x, parts.p)

  #parts.x, parts.v = bx.wrap_boundary(parts.x, parts.v)
  parts.v = therm.step(parts.v, parts.mass)
  Tcorr = Corr_temp_eff(parts.mass[0]*parts.v, gorces)
  rmax = max([sum((parts.x[i,:] - 2)**2) for i in range(len(parts.x[:,0]))])
  return parts.x, parts.v, e, ecaja, np.sqrt(rmax), Tcorr

# Interaction
qo = 1.664
po = 120
h_barra = 196.727394
D = 207*(h_barra/(po*qo))**3
scut = 20
pauli = pexmd.interaction.Pauli(scut, D, qo, po)
# Box
L = 2#5*qo
bx = pexmd.box.Box(np.zeros(3), L*np.ones(3), 'Fixed')
# Thermostat
Q = 0.02
To = np.concatenate((np.linspace(8, 2, 4), np.linspace(1.0, 0.2, 5)))
therm = pexmd.thermostat.Andersen(To[0], Q)
# Particles
Npart = 50#30
m = 0.5109989461
parts = pexmd.particles.PointParticles(Npart)
parts.x = particulas(Npart, L)
parts.v = np.random.normal(0,np.sqrt(To[0]*1.5), (Npart,3))
parts.mass = m
# Integrator
h = 0.0005
integ = pexmd.integrator.VelVerlet(h)
# Potencial de caja
bx = lambda x: fuerzas_caja(x, 0.5*L, 10, 20)

Nterm = 5000
Ntemp = len(To)
Epot = np.zeros(Ntemp*Nterm)
Epotcaja = np.zeros(Ntemp*Nterm)
Ecin = np.zeros(Ntemp*Nterm)
Tcorr = np.zeros(Ntemp*Nterm)
rmax = np.zeros(Ntemp*Nterm)
t = time.time()
i = 0
print("Termalizacion inicial")
for k in range(Nterm):
  parts.x, parts.v, a, b, c, d = avanzar_fp(parts, pauli, integ, bx, therm)
for j in range(Ntemp):
  therm = pexmd.thermostat.Andersen(To[j], 0)
  print("Temperatura:", To[j])
  for k in range(Nterm):
    parts.x, parts.v, Epot[i], Epotcaja[i], rmax[i], Tcorr[i] = avanzar_fp(parts, pauli, integ, bx, therm)
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
plt.figure()
plt.plot(rmax)
plt.figure()
plt.plot(T, "r-")
plt.plot(Tcorr, "b--")
plt.plot(T+Tcorr, "k-")
plt.plot(T)
plt.figure()
plt.plot(Ecin, 'r')
plt.plot(Epot, 'b')
plt.plot(Epotcaja, 'b--')
plt.plot(Ecin+Epot, 'k')
plt.figure()
plt.errorbar(To, Teff_mean, Teff_std)
plt.show()
