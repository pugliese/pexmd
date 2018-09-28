import numpy as np
import pexmd
import matplotlib.pylab as plt
import sys
import itertools as it

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

def fuerzas_caja(x, L, k):
  fuerzas = k*(0.5*L - x)
  epot = 0.5*k*sum(sum((0.5*L - x)**2))
  return fuerzas, epot

def avanzarRK(parts, pauli, integ, bx, therm):
  parts.f, e = pauli.forces(parts.x, parts.p)
  gorces, e = pauli.gorces(parts.x, parts.p)
  x_temp, v_temp = integ.first_step(parts.x, parts.v, gorces, parts.a)
  parts.f, e = pauli.forces(x_temp, v_temp)
  gorces, e = pauli.gorces(x_temp, parts.mass[0]*v_temp)
  parts.x, parts.v = integ.last_step(parts.x, parts.v, gorces, parts.a)
  parts.x, parts.v = bx.wrap_boundary(parts.x, parts.v)
  parts.v = therm.step(parts.v, parts.mass)
  return parts.x, parts.v, e

def avanzar(parts, pauli, integ, bx, therm):
  parts.f, e = pauli.forces(parts.x, parts.p)
  force_aux, ecaja = fuerzas_caja(parts.x, 4, 10)
  gorces, e = pauli.gorces(parts.x, parts.p)
  parts.x, parts.v = integ.step(parts.x, parts.v, gorces, parts.a + force_aux/parts.mass[0])
  #parts.x, parts.v = bx.wrap_boundary(parts.x, parts.v)
  parts.v = therm.step(parts.v, parts.mass)
  rmax = max([sum((parts.x[i,:]-2)**2) for i in range(len(parts.x[:,0]))])
  return parts.x, parts.v, e, ecaja, np.sqrt(rmax)

# Interaction
qo = 1.664
po = 120
h_barra = 196.727394
D = -207*(h_barra/(po*qo))**3
scut = 20
pauli = pexmd.interaction.Pauli(scut, D, qo, po)
# Box
L = 4#5*qo
bx = pexmd.box.Box(np.zeros(3), L*np.ones(3), 'Fixed')
# Thermostat
Q = 0#0.05
To = .5
therm = pexmd.thermostat.Andersen(To, Q)
# Particles
Npart = 30
m = 0.5109989461
parts = pexmd.particles.PointParticles(Npart)
parts.x = particulas(Npart, L)
parts.v = np.random.normal(0,np.sqrt(To*1.5), (Npart,3))
parts.mass = m
# Integrator
h = 0.00005
integ = pexmd.integrator.Euler(h)

Nterm = 100000
Epot = np.zeros(Nterm)
Epotcaja = np.zeros(Nterm)
Ecin = np.zeros(Nterm)
rmax = np.zeros(Nterm)
for i in range(Nterm):
  if (i%1000==0):
    print(i)
  parts.x, parts.v, Epot[i], Epotcaja[i], rmax[i] = avanzar(parts, pauli, integ, bx, therm)
  Ecin[i] = energia_cinetica(parts.v, parts.mass)
T = Ecin/(1.5*Npart)

plt.figure()
plt.plot(rmax)
plt.figure()
plt.plot(Ecin, 'r')
plt.plot(Epot, 'b')
plt.plot(Epotcaja, 'b--')
plt.plot(Ecin+Epot+Epotcaja, 'g')
plt.show()
