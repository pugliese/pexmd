import sys
import time
import numpy as np
import matplotlib.pylab as plt
import pexmd
import itertools as it
import time

# Interaction
scut = 200
opcion = 1
if (opcion==1):
  po = 1
  qo = 1
  DD = 10000
  Nstep = 2500
else:
  qo = 1.664
  po = 120
  h_barra = 196.727394
  DD = 207*(h_barra/(po*qo))**3
  Nstep = 5000
pauli = pexmd.interaction.Pauli(scut, DD, qo, po)
# Particles
m = 1
parts = pexmd.particles.PointParticles(2)
parts.x = np.array([[6.0, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype = np.float32)
parts.v = np.array([[-4.0, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype = np.float32)
parts.mass = m
# Integrator
h = 0.0001
integ = pexmd.integrator.Euler(h)
# Box
# box = box.Box(-50*np.ones(3), 50*np.ones(3))
# Neighbour list

if (sys.argv[1] == "e"):
  h = 0.0001
  if (len(sys.argv)==3):
    Nstep = int(Nstep*h/float(sys.argv[2]))
    h = float(sys.argv[2])
  integ = pexmd.integrator.Euler(h)
  q = np.zeros(Nstep+1)
  p = np.zeros(Nstep+1)
  pot = np.zeros(Nstep)
  Ecin = np.zeros(Nstep)
  q[0] = parts.x[0, 0] - parts.x[1, 0]
  p[0] = parts.mass[0]*(parts.v[0, 0] - parts.v[1, 0])
  for i in range(Nstep):
    parts.f, pot[i] = pauli.forces(parts.x, parts.p)
    gorces, pot[i] = pauli.gorces(parts.x, parts.p)
    parts.x, parts.v = integ.step(parts.x, parts.v, gorces, parts.a)
    q[i+1] = parts.x[0, 0] - parts.x[1, 0]
    p[i+1] = parts.mass[0]*(parts.v[0, 0] - parts.v[1, 0])
    Ecin[i] = (parts.v[0, 0]**2 + parts.v[1, 0]**2)/2

  e_std = np.std(Ecin+pot)
  e_mean = np.mean(Ecin+pot)
  print(e_mean, e_std, e_std/e_mean)
  plt.figure()
  plt.plot(q, p)
  plt.xlabel(r"$\Delta q$")
  plt.ylabel(r"$\Delta p$")
  plt.figure()
  plt.plot(Ecin, "r-")
  plt.plot(pot, "b-")
  plt.plot(Ecin+pot, "k-")
  plt.xlabel("Paso")
  plt.title("Energia")
  plt.legend(["Cinetica", "Potencial", "Total"], loc=6)
  plt.show()


if (sys.argv[1] == "erk"):
  h = 0.0001
  if (len(sys.argv)==3):
    Nstep = int(Nstep*h/float(sys.argv[2]))
    h = float(sys.argv[2])
  integ = pexmd.integrator.RK2(h)
  q = np.zeros(Nstep+1)
  p = np.zeros(Nstep+1)
  pot = np.zeros(Nstep)
  Ecin = np.zeros(Nstep)
  q[0] = parts.x[0, 0] - parts.x[1, 0]
  p[0] = parts.mass[0]*(parts.v[0, 0] - parts.v[1, 0])
  for i in range(Nstep):
    parts.f, pot[i] = pauli.forces(parts.x, parts.p)
    gorces, pot[i] = pauli.gorces(parts.x, parts.p)
    x_temp, v_temp = integ.first_step(parts.x, parts.v, gorces, parts.a)
    parts.f, pot[i] = pauli.forces(x_temp, v_temp)
    gorces, pot[i] = pauli.gorces(x_temp, parts.mass[0]*v_temp)
    parts.x, parts.v = integ.last_step(parts.x, parts.v, gorces, parts.a)
    q[i+1] = parts.x[0, 0] - parts.x[1, 0]
    p[i+1] = parts.mass[0]*(parts.v[0, 0] - parts.v[1, 0])
    Ecin[i] = (parts.v[0, 0]**2 + parts.v[1, 0]**2)/2

  e_std = np.std(Ecin+pot)
  e_mean = np.mean(Ecin+pot)
  print(e_mean, e_std, e_std/e_mean)
  plt.figure()
  plt.plot(q, p)
  plt.xlabel(r"$\Delta q$")
  plt.ylabel(r"$\Delta p$")
  plt.figure()
  plt.plot(Ecin, "r-")
  plt.plot(pot, "b-")
  plt.plot(Ecin+pot, "k-")
  plt.xlabel("Paso")
  plt.title("Energia")
  plt.legend(["Cinetica", "Potencial", "Total"], loc=6)
  plt.show()

def trayectoria(qi, pi, pauli, integ, parts):
  assert(pi*qi<=0)
  parts.x = np.array([[qi, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype = np.float32)
  parts.v = np.array([[pi, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype = np.float32)
  q = [parts.x[0, 0] - parts.x[1, 0]]
  p = [parts.mass[0]*(parts.v[0, 0] - parts.v[1, 0])]
  while (abs(q[-1]) <= abs(qi) and abs(p[-1]) <= abs(pi)):
    parts.f, e = pauli.forces(parts.x, parts.p)
    gorces, e = pauli.gorces(parts.x, parts.p)
    parts.x, parts.v = integ.step(parts.x, parts.v, gorces, parts.a)
    q.append(parts.x[0, 0] - parts.x[1, 0])
    p.append(parts.mass[0]*(parts.v[0, 0] - parts.v[1, 0]))
  reboto = (0<q[-1]*qi)
  return q, p, reboto

if (sys.argv[1] == "f"):
  Nq = 49
  if (3<=len(sys.argv)):
    Nq = int(sys.argv[2])
  po = np.linspace(-6, 0, Nq)
  fasesq = []
  fasesp = []
  for pi in po:
    print(pi)
    q, p, r = trayectoria(5, pi, pauli, integ, parts)
    fasesq.append(q)
    fasesp.append(p)
  po = np.linspace(0, 6, Nq)
  for pi in po:
    print(pi)
    q, p, r = trayectoria(-5, pi, pauli, integ, parts)
    fasesq.append(q)
    fasesp.append(p)
  if (len(sys.argv)==4):
    filename = sys.argv[3]
  else:
    filename = "datafases.txt"
  f = open(filename, "w")
  for i in range(len(fasesq)):
    for j in range(len(fasesq[i])-1):
      f.write(str(fasesq[i][j]) + " ")
    f.write(str(fasesq[i][-1])+"\n")
    for j in range(len(fasesp[i])-1):
      f.write(str(fasesp[i][j]) + " ")
    f.write(str(fasesp[i][-1])+"\n")
  f.close()
  for i in range(len(fasesq)):
    plt.plot(fasesq[i], fasesp[i], "b-")
  plt.xlabel(r"$\Delta q$")
  plt.ylabel(r"$\Delta p$")
  plt.axis([-5,5,-6, 6])
  plt.title("Espacio de fases")
  plt.show()



def trayectoria_rk(qi, pi, pauli, integ, parts):
  assert(pi*qi<=0)
  parts.x = np.array([[qi, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype = np.float32)
  parts.v = np.array([[pi, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype = np.float32)
  q = [parts.x[0, 0] - parts.x[1, 0]]
  p = [parts.mass[0]*(parts.v[0, 0] - parts.v[1, 0])]
  while (abs(q[-1])<=abs(qi) and abs(p[-1])<=6):
    parts.f, e = pauli.forces(parts.x, parts.p)
    gorces, e = pauli.gorces(parts.x, parts.p)
    x_temp, v_temp = integ.first_step(parts.x, parts.v, gorces, parts.a)
    parts.f, e = pauli.forces(x_temp, v_temp)
    gorces, e = pauli.gorces(x_temp, parts.mass[0]*v_temp)
    parts.x, parts.v = integ.last_step(parts.x, parts.v, gorces, parts.a)
    q.append(parts.x[0, 0] - parts.x[1, 0])
    p.append(parts.mass[0]*(parts.v[0, 0] - parts.v[1, 0]))
  reboto = (0<q[-1]*qi)
  return q, p, reboto

if (sys.argv[1] == "frk"):
  h = 0.0001
  integ = pexmd.integrator.RK2(h)
  Nq = 25
  if (3<=len(sys.argv)):
    Nq = int(sys.argv[2])
  po = np.linspace(-6, 0, Nq)
  fasesq = []
  fasesp = []
  energia = []
  for pi in po:
    q, p, r = trayectoria_rk(5, pi, pauli, integ, parts)
    fasesq.append(q)
    fasesp.append(p)
  po = np.linspace(0, 6, Nq)
  for pi in po:
    print(pi)
    q, p, r = trayectoria_rk(-5, pi, pauli, integ, parts)
    fasesq.append(q)
    fasesp.append(p)
  if (len(sys.argv)==4):
    filename = sys.argv[3]
  else:
    filename = "datafases.txt"
  f = open(filename, "w")
  for i in range(len(fasesq)):
    for j in range(len(fasesq[i])-1):
      f.write(str(fasesq[i][j]) + " ")
    f.write(str(fasesq[i][-1])+"\n")
    for j in range(len(fasesp[i])-1):
      f.write(str(fasesp[i][j]) + " ")
    f.write(str(fasesp[i][-1])+"\n")
  f.close()
  for i in range(len(fasesq)):
    plt.plot(fasesq[i], fasesp[i], "b-")
  plt.xlabel(r"$\Delta q$")
  plt.ylabel(r"$\Delta p$")
  plt.axis([-5,5,-6, 6])
  plt.title("Espacio de fases")
  plt.show()

def FixedPoint(parts, interact, dt, Niter = 5):
  Z_x = parts.x
  Z_v = parts.v
  for k in range(Niter):
    #forces, e = interact.forces(Z_x, parts.mass[0]*Z_v)
    #gorces, e = interact.gorces(Z_x, parts.mass[0]*Z_v)
    forces, gorces, e = interact.fgorces(Z_x, parts.mass[0]*Z_v)
    Z_x = parts.x + .5*dt*(Z_v - gorces)
    Z_v = parts.v + .5*dt*forces/parts.mass[0]
  return 2*Z_x - parts.x, 2*Z_v - parts.v

def trayectoria_fp(qi, pi, pauli, dt, parts):
  assert(pi*qi<=0)
  parts.x = np.array([[qi, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype = np.float32)
  parts.v = np.array([[pi, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype = np.float32)
  q = [parts.x[0, 0] - parts.x[1, 0]]
  p = [parts.mass[0]*(parts.v[0, 0] - parts.v[1, 0])]
  while (abs(q[-1])<=abs(qi)*1.05 and abs(p[-1])<=abs(pi)*1.05):
    parts.x, parts.v = FixedPoint(parts, pauli, dt)
    q.append(parts.x[0, 0] - parts.x[1, 0])
    p.append(parts.mass[0]*(parts.v[0, 0] - parts.v[1, 0]))
  reboto = (0<q[-1]*qi)
  return q, p, reboto

if (sys.argv[1] == "fp"):
  h = 0.001
  integ = pexmd.integrator.RK2(h)
  Nq = 25
  if (3<=len(sys.argv)):
    Nq = int(sys.argv[2])
  qinit = 5*qo
  pmax = 6*po
  pmin = 0
  pinits = np.linspace(-pmax, -pmin, Nq)
  pinits = np.array(pinits, dtype = np.float32)
  fasesq = []
  fasesp = []
  energia = []
  for pi in pinits:
    t = time.time()
    q, p, r = trayectoria_fp(qinit, pi, pauli, h, parts)
    fasesq.append(q)
    fasesp.append(p)
    t = time.time()-t
    print(pi, ":", t)
  pinits = np.linspace(pmin, pmax, Nq)
  pinits = np.array(pinits, dtype = np.float32)
  for pi in pinits:
    t = time.time()
    q, p, r = trayectoria_fp(-qinit, pi, pauli, h, parts)
    fasesq.append(q)
    fasesp.append(p)
    t = time.time()-t
    print(pi, ":", t)
  if (len(sys.argv)==4):
    filename = sys.argv[3]
  else:
    filename = "datafases.txt"
  f = open(filename, "w")
  for i in range(len(fasesq)):
    for j in range(len(fasesq[i])-1):
      f.write(str(fasesq[i][j]) + " ")
    f.write(str(fasesq[i][-1])+"\n")
    for j in range(len(fasesp[i])-1):
      f.write(str(fasesp[i][j]) + " ")
    f.write(str(fasesp[i][-1])+"\n")
  f.close()
  for i in range(len(fasesq)):
    plt.plot(fasesq[i], fasesp[i], "b-")
  plt.xlabel(r"$\Delta q$")
  plt.ylabel(r"$\Delta p$")
  plt.axis([-qinit, qinit, -pmax, pmax])
  plt.title("Espacio de fases")
  plt.savefig("fases_{0}_{1}.png".format(po, qo))
  plt.show()

if (sys.argv[1] == "efp"):
  h = 0.001
  if (len(sys.argv)==3):
    Nstep = int(Nstep*h/float(sys.argv[2]))
    h = float(sys.argv[2])
  q = np.zeros(Nstep+1)
  p = np.zeros(Nstep+1)
  pot = np.zeros(Nstep)
  Ecin = np.zeros(Nstep)
  rho = np.zeros(Nstep)
  q[0] = parts.x[0, 0] - parts.x[1, 0]
  p[0] = parts.mass[0]*(parts.v[0, 0] - parts.v[1, 0])
  t = time.time()
  for i in range(Nstep):
    parts.x, parts.v = FixedPoint(parts, pauli, h, 5)
    parts.f, pot[i] = pauli.forces(parts.x, parts.p)
    q[i+1] = parts.x[0, 0] - parts.x[1, 0]
    p[i+1] = parts.mass[0]*(parts.v[0, 0] - parts.v[1, 0])
    rho[i] = (q[i+1]-q[i])*(p[i+1]-p[i])
    Ecin[i] = (parts.v[0, 0]**2 + parts.v[1, 0]**2)/2
  t = time.time() - t
  print(t)
  e_std = np.std(Ecin+pot)
  e_mean = np.mean(Ecin+pot)
  print(e_mean, e_std, e_std/e_mean)
  plt.figure()
  plt.plot(q, p)
  plt.xlabel(r"$\Delta q$")
  plt.ylabel(r"$\Delta p$")
  plt.figure()
  plt.plot(rho)
  plt.xlabel('Paso')
  plt.ylabel(r"$\rho$")
  plt.figure()
  plt.plot(Ecin, "r-")
  plt.plot(pot, "b-")
  plt.plot(Ecin+pot, "k-")
  plt.xlabel("Paso")
  plt.title("Energia")
  plt.legend(["Cinetica", "Potencial", "Total"], loc=6)
  plt.show()


def barrido_D(D, tol=1E-2):
  filename = "barridoD_{0}.png".format(D)
  print(filename)
  return barrido(1, 1, filename, D, tol)

def barrido_PC(P,C, tol=1E-2):
  qo = np.sqrt(P*C)
  po = np.sqrt(P/C)
  filename = "barrido_{0}_{1}.png".format(qo*po, qo/po)
  return barrido(qo, po, filename, 10000, tol)

def barrido(qo, po, filename, D=10000, tol=1E-2, N=10):
  pauli = pexmd.interaction.Pauli(200, D, qo, po)
  parts = pexmd.particles.PointParticles(2)
  parts.x = np.array([[2.0, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype = np.float32)
  parts.v = np.array([[-2.0, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype = np.float32)
  parts.mass = 1
  h = 0.0001
  integ = pexmd.integrator.Euler(h)

  qinit = 5*qo
  pinf = 0    # Trayectoria que rebota
  psup = -3*po   # Trayectoria que no rebota
  q, p, r = trayectoria_fp(qinit, psup, pauli, integ.dt, parts)
  qs = []
  ps = []
  # Valores iniciales
  while(r):
    pinf = psup
    psup = 2*psup
    q, p, r = trayectoria_fp(qinit, psup, pauli, integ.dt, parts)
    print(pinf)
  pmax = abs(psup)
  qs.append(q)
  ps.append(p)
  # Comienzo busqueda
  i = 0
  while(abs(psup-pinf) > tol):
    t = time.time()
    pmed = (psup+pinf)/2
    q, p, r = trayectoria_fp(qinit, pmed, pauli, integ.dt, parts)
    qs.append(q)
    ps.append(p)
    if(r):
      pinf = pmed
    else:
      psup = pmed
    print("{0}) {1}: {2} segs".format(i, abs(psup-pinf), time.time()-t))
    i += 1
  plt.figure()
  for i in range(len(qs)):
    plt.plot(-np.array(qs[i]), -np.array(ps[i]), "b-")
    plt.plot(qs[i],ps[i], "b-")
  plt.xlabel(r'$\Delta q$')
  plt.ylabel(r'$\Delta p$')
  plt.axis([-qinit, qinit, -pmax, pmax])
  plt.savefig(filename)
  plt.close()#plt.show()
  return pmed

def avanzarN(parts, pot, dt, Nstep):
  for i in range(Nstep):
    parts.x, parts.v = FixedPoint(parts, pot, dt)
  return parts.x, parts.v

def vol_fases(qs, ps, Nskip, Nsamp, h=1E-3):
  n = len(qs)
  m = len(ps)
  res = np.zeros((n, m, Nsamp, 2))
  parts = pexmd.particles.PointParticles(2)
  pauli = pexmd.interaction.Pauli(200, 10000, 1, 1)
  parts.mass = 1
  for i in range(n):
    for j in range(m):
      parts.x = np.array([[qs[i], 0.0, 0.0], [0.0, 0.0, 0.0]], dtype = np.float32)
      parts.v = np.array([[-ps[j], 0.0, 0.0], [0.0, 0.0, 0.0]], dtype = np.float32)
      for k in range(Nsamp):
        avanzarN(parts, pauli, h, Nskip)
        res[i, j, k, :] = [parts.x[0, 0] - parts.x[1, 0], parts.mass[0]*(parts.v[0, 0] - parts.v[1, 0])]
  plt.figure()
  for i in range(n):
    for j in range(m):
      plt.plot(res[i, j, :, 0], res[i, j, :, 1], 'k.')
  plt.xlabel(r'$\Delta q$')
  plt.ylabel(r'$\Delta p$')
  plt.show()
  return res

def grilla():

def lado(x, y):
  n = len(x)
  i = int(n/2)
  while(0<i and 0<=x[i]):
    i -= 1
    while(i<n-1 and x[i]<=0):
      i += 1
      if (x[i] == 0 or i==n-1):
        return y[i]
      else:
        return y[i]-x[i]*(y[i+1]-y[i])/(x[i+1]-x[i])

def area(D, qo = 1, po = 1, N = 10):
  pauli = pexmd.interaction.Pauli(200, D, qo, po)
  parts = pexmd.particles.PointParticles(2)
  parts.x = np.array([[2.0, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype = np.float32)
  parts.v = np.array([[-2.0, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype = np.float32)
  parts.mass = 1
  h = 0.001
  integ = pexmd.integrator.Euler(h)

  pinf = 0    # Trayectoria que rebota; obtengo lq
  psup = -2   # Trayectoria que no rebota; obtengo lp
  q, p, r = trayectoria(5, psup, pauli, integ, parts)
  lq = []
  lp = []
  # Valores iniciales
  while(r):
    pinf = psup
    psup = 2*psup
    q, p, r = trayectoria(5, psup, pauli, integ, parts)
  lp.append(abs(lado(q, p)))
  if(pinf==0):
    lq.append(5)
  else:
    q, p, r = trayectoria(5, pinf, pauli, integ, parts)
    lq.append(abs(lado(p, q)))
  # Comienzo busqueda
  for i in range(N):
    pmed = (psup+pinf)/2
    q, p, r = trayectoria(5, pmed, pauli, integ, parts)
    if(r):
      pinf = pmed
      lq.append(abs(lado(p, q)))
      lp.append(lp[-1])
    else:
      psup = pmed
      lp.append(abs(lado(q, p)))
      lq.append(lq[-1])
  q, p, r = trayectoria(5, psup, pauli, integ, parts)
  plt.plot(q,p, "b-")
  print(p)
  print([-pi for pi in p])
  plt.plot(q,[-pi for pi in p], "b-")
  q, p, r = trayectoria(5, pinf, pauli, integ, parts)
  plt.plot(q,p, "r-")
  plt.plot([-qi for qi in q],p, "r-")
  plt.show()
  return lq, lp, lq[-1]*lp[-1]*np.pi

def areas_D(D, N=10):
  res = []
  i = 1
  for d in D:
    print(i, d)
    x = area(d, 1, 1, N)
    res.append(x[2])
    i += 1
  return res

def areas_qp(qp, D = 0.6, N = 10):
  res = []
  i = 1
  for qpi in qp:
    print(i, qo)
    x = area(D, qi, 1.0/qi, N)
    res.append(x[2])
    i += 1
  return res

if (sys.argv[1] == "vs"):
  Niter = 10000
  # Interaction
  po = 1
  qo = 1
  scut = 200
  D = 10000
  pauli = pexmd.interaction.Pauli(scut, D, qo, po)
  # Particles
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
  Npart = 1000
  m = 1
  parts = pexmd.particles.PointParticles(Npart)
  parts.x = particulas(Npart, 5)
  parts.mass = m
  pairs = np.array(list(it.combinations(range(Npart), 2)), dtype=np.int64)

  forces1, gorces1, e1 = pauli.fgorces(parts.x, parts.v, pairs)
  forces2, e2 = pauli.forces(parts.x, parts.v, pairs)
  gorces2, e3 = pauli.gorces(parts.x, parts.v, pairs)
  assert(np.std(forces1-forces2) < 1E-6)
  assert(np.std(gorces1-gorces2) < 1E-6)
  assert(e1==e2)
  assert(e2==e3)
  t1 = time.time()
  for i in range(Niter):
    forces1, gorces1, e1 = pauli.fgorces(parts.x, parts.v, pairs)
  t1 = time.time() - t1
  t2 = time.time()
  for i in range(Niter):
    forces2, e2 = pauli.forces(parts.x, parts.v, pairs)
    gorces2, e3 = pauli.gorces(parts.x, parts.v, pairs)
  t2 = time.time() - t2
  print(t1, t2, t2/t1, t1/t2)
