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
#m = 1
#parts = pexmd.particles.PointParticles(2)
#parts.x = np.array([[6.0, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype = np.float32)
#parts.v = np.array([[-4.0, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype = np.float32)
#parts.mass = m
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
  parts.v = np.array([[pi/parts.mass[0], 0.0, 0.0], [0.0, 0.0, 0.0]], dtype = np.float32)
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
  #plt.title("Espacio de fases")
  #plt.savefig("fases_{0}_{1}.png".format(po, qo))
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

# ------------------ Muestreo -------------------------- #

def muestreo(qinit, pmax, p_step, D, h=1E-3, rebs = False, pmin = 0):
  #pauli = pexmd.interaction.Pauli(200, D, 1, 1)
  pauli = pexmd.interaction.Pauli(200, D, 1.664, 120)
  h = 0.001
  #Nq = int(pmax/p_step) + 1
  pinits = np.arange(-pmax, -pmin, p_step)
  #pinits = np.linspace(-pmax, 0, Nq)
  pinits = np.array(pinits, dtype = np.float32)
  parts = pexmd.particles.PointParticles(2)
  parts.mass = 0.51
  fasesq = []
  fasesp = []
  rs = []
  energia = []
  for pi in pinits:
    t = time.time()
    q, p, r = trayectoria_fp(qinit, pi, pauli, h, parts)
    fasesq.append(q)
    fasesp.append(p)
    rs.append(r)
    t = time.time()-t
    print(pi, ":", t)
  fasesq = np.array(fasesq)
  fasesp = np.array(fasesp)
  if rebs:
    return fasesq, fasesp, rs
  else:
    return fasesq, fasesp

def muestreo_dif(qinit, pmax, p_step1, p_step2, D, h=1E-3):
  fasesq, fasesp, rs = muestreo(qinit, pmax, p_step1, D, h, True)
  i = len(rs) - sum(rs)
  print(rs)
  print(i, rs[i-1], rs[i], fasesp[i-1][0], fasesp[i][0])
  fasesq2, fasesp2 = muestreo(qinit, -fasesp[i-1][0], p_step2, D, h, False, -fasesp[i][0])
  fasesq = np.append(fasesq, fasesq2)
  fasesp = np.append(fasesp, fasesp2)
  return fasesq, fasesp



def plot_muestreo(fasesq, fasesp):
  for i in range(len(fasesq)):
    plt.plot(fasesq[i], fasesp[i], "b-")
    plt.plot(-np.array(fasesq[i]), -np.array(fasesp[i]), "b-")
  plt.xlabel(r"$\Delta q$")
  plt.ylabel(r"$\Delta p$")
  plt.axis([-fasesq[0][0], fasesq[0][0], -fasesp[0][0], fasesp[0][0]])
# ------------------ Barridos -------------------------- #

def barrido_D(D, tol=1E-2):
  filename = "barridoD_{0}.png".format(D)
  print(filename)
  return barrido(1, 1, filename, D, tol)

def barrido_PC(P,C, tol=1E-2):
  qo = np.sqrt(P*C)
  po = np.sqrt(P/C)
  filename = "barrido_{0}_{1}.png".format(qo*po, qo/po)
  return barrido(qo, po, filename, 10000, tol)

def barrido(qo, po, filename="None", D=10000, tol=1E-2):
  pauli = pexmd.interaction.Pauli(200, D, qo, po)
  parts = pexmd.particles.PointParticles(2)
  parts.x = np.array([[2.0, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype = np.float32)
  parts.v = np.array([[-2.0, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype = np.float32)
  parts.mass = 1
  h = 0.001
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
  if filename != "None":
    plt.figure()
    for i in range(len(qs)):
      plt.plot(-np.array(qs[i]), -np.array(ps[i]), "b-")
      plt.plot(qs[i],ps[i], "b-")
    plt.xlabel(r'$\Delta q$')
    plt.ylabel(r'$\Delta p$')
    plt.axis([-qinit, qinit, -pmax, pmax])
    plt.savefig(filename)
    plt.close()#plt.show()
  q, p, r = trayectoria_fp(qinit, pinf, pauli, integ.dt, parts)
  return np.array(q), np.array(p)


def cuadraturas(qs, ps):
  i0 = sum(ps<0)   # 0<ps[i] para todo i2<i
  qmin1 = min(qs[:i0])
  i1 = max((qs[:i0] == qmin1)*range(i0))
  qmin2 = min(qs[i0:])
  i2 = max((qs[i0:] == qmin2)*range(i0,len(qs)))
  trapecios = (qs[i1+1:i2] - qs[i1:i2-1])*(ps[i1+1:i2] + ps[i1:i2-1])  # No pongo el 1/2 porque sumo cada area 2 veces
  return abs(sum(trapecios))


def area(Ds, save=False, tol=1E-2):
  A = []
  if (type(tol)!=list and type(tol)!= np.array):
    tols = np.ones(len(Ds))*tol
  else:
    assert(len(tol)==len(Ds))
    tols = tol
  if save:
    names = "fases_D_{0}.png"
  else:
    names = "None"
  for i in range(len(Ds)):
    print("Barrido para D={0}".format(Ds[i]))
    q, p = barrido(qo, po, names.format(Ds[i]), Ds[i], tols[i])
    A.append(cuadraturas(q,p))
  return np.array(A)

# Calculo numerico del area de exclusion teorica
def area_teo(D, N):
  Hc = 0.5*(1+np.log(2*D))
  xs = np.arange(N)*(2*Hc-1)/N+1
  area = 4*(2*Hc-1)*np.sum(np.sqrt((xs-1-np.log(xs))/(2*Hc-xs)))/N
  return area

# ------------------ Volumen fases -------------------------- #

def avanzarN_fp(parts, pot, dt, Nstep):
  for i in range(Nstep):
    parts.x, parts.v = FixedPoint(parts, pot, dt)
  return parts.x, parts.v

def avanzarN_E(parts, pot, integ, Nstep):
  for i in range(Nstep):
    parts.f, gorces, e = pot.fgorces(parts.x, parts.v)
    parts.x, parts.v = integ.step(parts.x, parts.v, gorces, parts.a)
  return parts.x, parts.v

def vol_fases(qs, ps, Nskip, Nsamp, h=1E-3):
  n = len(qs)
  m = len(ps)
  res_e = np.zeros((n, m, Nsamp, 2))
  res_fp = np.zeros((n, m, Nsamp, 2))
  area_e = np.zeros(Nsamp)
  area_fp = np.zeros(Nsamp)
  parts = pexmd.particles.PointParticles(2)
  pauli = pexmd.interaction.Pauli(200, 10000, 1, 1)
  parts.mass = 1
  integ = pexmd.integrator.Euler(h)
  for i in range(n):
    for j in range(m):
      print(i*m+j+1, "/", n*m)
      parts.x = np.array([[qs[i], 0.0, 0.0], [0.0, 0.0, 0.0]], dtype = np.float32)
      parts.v = np.array([[-ps[j], 0.0, 0.0], [0.0, 0.0, 0.0]], dtype = np.float32)
      for k in range(Nsamp):
        avanzarN_fp(parts, pauli, h, Nskip)
        res_fp[i, j, k, :] = [parts.x[0, 0] - parts.x[1, 0], parts.mass[0]*(parts.v[0, 0] - parts.v[1, 0])]
      parts.x = np.array([[qs[i], 0.0, 0.0], [0.0, 0.0, 0.0]], dtype = np.float32)
      parts.v = np.array([[-ps[j], 0.0, 0.0], [0.0, 0.0, 0.0]], dtype = np.float32)
      for k in range(Nsamp):
        avanzarN_E(parts, pauli, integ, Nskip)
        res_e[i, j, k, :] = [parts.x[0, 0] - parts.x[1, 0], parts.mass[0]*(parts.v[0, 0] - parts.v[1, 0])]
  plt.figure()
  for k in range(Nsamp):
    matrix, area_e[k], reg_q, reg_p = grilla(res_e[:, :, k, 0].reshape((n*m)), res_e[:, :, k, 1].reshape((n*m)))
    matrix, area_fp[k], reg_q, reg_p = grilla(res_fp[:, :, k, 0].reshape((n*m)), res_fp[:, :, k, 1].reshape((n*m)))
    #plt.plot(reg_q, reg_p, "b.")

  for i in range(n):
    for j in range(m):
      plt.plot(res_e[i, j, :, 0], res_e[i, j, :, 1], 'ks')
      plt.plot(res_fp[i, j, :, 0], res_fp[i, j, :, 1], 'ro')
  plt.xlabel(r'$\Delta q$')
  plt.ylabel(r'$\Delta p$')
  plt.title("Espacio de fases")

  plt.figure()
  plt.plot(area_e, 'ro--')
  plt.plot(area_fp, 'bo--')
  plt.legend(['Euler', 'MPR'])
  plt.ylabel('Area')
  plt.show()
  return

def grilla(qs, ps):
  dq = 2E-4
  dp = 2E-4
  radio = 25
  n = len(qs)
  qo = min(qs)-(radio+1)*dq
  po = min(ps)-(radio+1)*dp
  qmax = max(qs)+(radio+1)*dq
  pmax = max(ps)+(radio+1)*dp
  N = int((qmax-qo)/dq) + 1
  M = int((pmax-po)/dp) + 1
  matrix = np.zeros((N,M))
  region_q = []
  region_p = []
  for k in range(n):
    i = int(np.round((qs[k]-qo)/dq))
    j = int(np.round((ps[k]-po)/dp))
    for l1 in range(-radio, radio+1):
      for l2 in range(-radio, radio+1):
        ii = i+l1
        jj = j+l2
        if (matrix[ii, jj] == 0 and l1**2+l2**2 <= radio**2):
          region_q.append(qo + ii*dq)
          region_p.append(po + jj*dp)
          matrix[ii, jj] = 1
  return matrix, sum(sum(matrix))*dq*dp, region_q, region_p

# --------------------- Para el VMD ------------------------- #

def data_vmd(filename="coso.lammpstrj", pi=2, Nsamp=500, Nskip = 5):
  parts = pexmd.particles.PointParticles(2)
  pauli = pexmd.interaction.Pauli(200, 10000, 1, 1)
  parts.x = np.array([[6, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype = np.float32)
  parts.v = np.array([[-pi, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype = np.float32)
  parts.mass = 1
  h = 1E-3
  header = "ITEM: TIMESTEP\n{0}\nITEM: NUMBER OF ATOMS\n2\nITEM: BOX BOUNDS pp pp pp\n0.0000000000000000e+00 1.2000000000000000e+00\n0.0000000000000000e+00 1.2000000000000000e+00\n0.0000000000000000e+00 1.2000000000000000e+00\nITEM: ATOMS type x y z vx vy vz\n"
  f = open(filename, "w")
  for i in range(Nsamp):
    print(i)
    parts.x, parts.v = avanzarN_fp(parts, pauli, h, Nskip)
    f.write(header.format(i*Nskip))
    f.write("1 {0} 0 0 {1} 0 0 \n".format(parts.x[0][0], parts.v[0][0]))
    f.write("1 {0} 0 0 {1} 0 0 \n".format(parts.x[1][0], parts.v[1][0]))
  f.close()

def make_movie(name, pi, Nstep=200, Nskip = 10):
  parts = pexmd.particles.PointParticles(2)
  pauli = pexmd.interaction.Pauli(200, 10000, 1, 1)
  parts.x = np.array([[5, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype = np.float32)
  parts.v = np.array([[-pi, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype = np.float32)
  parts.mass = 1
  h = 1E-3
  for i in range(Nstep):
    print(i)
    parts.x, parts.v = avanzarN_fp(parts, pauli, h, Nskip)
    plt.figure()
    plt.plot([parts.x[0,0]], [0], "ro")
    plt.plot([parts.x[1,0]], [0], "bo")
    plt.axis([-5, 5, -1, 1])
    plt.axis('off')
    if (i<10):
      plt.savefig("Videos/{0}_0{1}".format(name, i))
    else:
      plt.savefig("Videos/{0}_{1}".format(name, i))
    plt.close()



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

def policircular(qs, ps):
  qcm = sum(qs)/len(qs)
  pcm = sum(ps)/len(ps)
  Nangles = 1000
  thetas = np.linspace(0, 2*np.pi, Nangles)
  dtheta = thetas[1] - thetas[0]
  n = len(qs)
  rs = np.zeros(n)
  dists = np.zeros(n-1)
  for i in range(n):
    dists[0:i] = (qs[i] - qs[0:i])**2 + (ps[i] - ps[0:i])**2
    dists[i:n] = (qs[i] - qs[i+1:])**2 + (ps[i] - ps[i+1:])**2
    rs[i] = np.sqrt(min(dists))/2
  dists_cm = np.sqrt((qs-qcm)**2 + (ps-pcm)**2)
  q_exts = qs + rs*(qs-qcm)/dists_cm
  p_exts = ps + rs*(ps-pcm)/dists_cm
  radios = np.zeros(len(thetas))
  area = 0
  for i in range(len(thetas)):
    cos_diff = np.cos(thetas[i])*(qs-qcm)/dists_cm + np.sin(thetas[i])*(ps-pcm)/dists_cm
    radios[i] = max(((q_exts-qcm)*np.cos(thetas[i]) + (p_exts-pcm)*np.sin(thetas[i]))*cos_diff)
    area += 0.5*dtheta*radios[i]**2
  return area, radios


def policontorno(qs, ps):
  n = len(qs)
  qcm = sum(qs)/len(qs)
  pcm = sum(ps)/len(ps)
  Nangles = 1000
  thetas = np.linspace(0, 2*np.pi, Nangles)
  dtheta = thetas[1] - thetas[0]
  rs = bolas(qs, ps)
  dists_cm = np.sqrt((qs-qcm)**2 + (ps-pcm)**2)
  sin_alpha = (ps-pcm)/dists_cm
  cos_alpha = (qs-qcm)/dists_cm
  radios = np.zeros((Nangles, n))
  for i in range(n):
    sin_theta_alpha = np.sin(thetas)*cos_alpha[i] - sin_alpha[i]*np.cos(thetas)
    a = rs[i]
    c = dists_cm[i]
    for j in range(Nangles):
      if (c**2*sin_theta_alpha[j]**2 < a**2):
        radios[j, i] = c*np.sqrt(1-sin_theta_alpha[j]**2) + np.sqrt(a**2 - c**2*sin_theta_alpha[j]**2)
      else:
        #radios[j, i] = c*np.sqrt(1-sin_theta_alpha[j]**2)
        radios[j, i] = 0
  area = 0
  contorno = np.zeros(len(thetas))
  for i in range(len(thetas)):
    contorno[i] = max(radios[i, :])
  i = 0
  while (i<Nangles):
    j = 0
    k = 0
    while (i+j < Nangles and 0 < contorno[i+j]):
      j += 1
    if (i+j < Nangles):
      while (i+j+k < Nangles and contorno[i+j+k] == 0):
        k += 1
      for l in range(i+j, i+j+k):
        contorno[l] = (contorno[i+j+k]-contorno[i+j])*(thetas[l]-thetas[i+j])/(thetas[i+j+k]-thetas[i+j])+contorno[i+j]
    i = i+j+k
  for i in range(Nangles):
    for j in range(n):
      cmax = contorno[i]
      cmin = contorno[i]
      if ((cmax*np.cos(thetas[i])-qs[j])**2 + (cmax*np.sin(thetas[i])-ps[j])**2 < rs[j]**2):
        while ((cmax*np.cos(thetas[i])-qs[j])**2 + (cmax*np.sin(thetas[i])-ps[j])**2 < rs[j]**2):
          cmax *= 2
        for k in range(10):
          cmed = (cmax+cmin)/2
          if ((cmed*np.cos(thetas[i])-qs[j])**2 + (cmed*np.sin(thetas[i])-ps[j])**2 < rs[j]**2):
            cmin = cmed
          else:
            cmax = cmed
        contorno[i] = cmed
  a = 0
  p = 1
  promedio = ventana(a, Nangles)
  cont_aux = contorno
  for i in range(Nangles):
    cont_aux[i] = sum(promedio*contorno**p)
    promedio = np.roll(promedio, 1)
  contorno = cont_aux**(1/p)
  area = 0.5*dtheta*sum(contorno**2)
  return area, contorno

def ventana(a, N):
  res = np.zeros(N)
  res[0] = 1
  for i in range(a):
    res[i+1] = 1
    res[-i] = 1
  res = res/(2*a+1)
  return res

def vent_gauss(d, N):
  n = N/2
  res = np.exp(-(np.arange(N)-n)**2/d**2)
  res = res/sum(res)
  return res



def bolas(qs, ps):
  n = len(qs)
  rs = np.zeros(n)
  dists = np.zeros(n-1)
  for i in range(n):
    dists[0:i] = (qs[i] - qs[0:i])**2 + (ps[i] - ps[0:i])**2
    dists[i:n] = (qs[i] - qs[i+1:])**2 + (ps[i] - ps[i+1:])**2
    rs[i] = 0.5*np.sqrt(min(dists))
  return rs
