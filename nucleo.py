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
a = 5/6
nuclear = pexmd.interaction.QCNM(scut, Vo, r1, p1, r2, p2, d, a)

# Trampa 

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

def avanzar(parts, pairs, pauli, dt):
  parts.x, parts.p = FixedPoint(parts, pairs, pauli, dt)
  parts.f, parts.g, e = pauli.fgorces(parts.x, parts.p, pairs)
  parts.x, parts.p = bx.wrap_boundary(parts.x, parts.p)
