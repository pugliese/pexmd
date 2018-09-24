import numpy as np
import pexmd
import matplotlib.pylab as plt

# Interaction
qo = 1
po = 1
scut = 200
D = 100000
pauli = pexmd.interaction.Pauli(scut, D, qo, po)
# Particles
m = 1
parts = pexmd.particles.PointParticles(2)
parts.x = np.array([[5.0, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype = np.float32)
parts.v = np.array([[-5.0, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype = np.float32)
parts.mass = m
# Integrator
h = 0.001
integ = pexmd.integrator.Euler(h)
# Box
# box = box.Box(-50*np.ones(3), 50*np.ones(3))
# Neighbour list

Nstep = 3000
q = np.zeros(Nstep+1)
p = np.zeros(Nstep+1)
fuerzas = np.zeros(Nstep+1)
guerzas = np.zeros(Nstep+1)
q[0] = parts.x[0, 0] - parts.x[1, 0]
p[0] = parts.mass[0]*(parts.v[0, 0] - parts.v[1, 0])
for i in range(Nstep):
    parts.f, e = pauli.forces(parts.x, parts.p)
    gorces, e = pauli.gorces(parts.x, parts.p)
    parts.x, parts.v = integ.step(parts.x, parts.v, gorces, parts.a)
    q[i+1] = parts.x[0, 0] - parts.x[1, 0]
    p[i+1] = parts.mass[0]*(parts.v[0, 0] - parts.v[1, 0])
    fuerzas[i+1] = parts.f[0, 0]
    guerzas[i+1] = gorces[0, 0]

plt.figure()
plt.plot(q, p)
plt.figure()
plt.plot(q)
plt.plot(p)
plt.figure()
plt.plot(fuerzas)
plt.figure()
plt.plot(guerzas)
plt.show()
