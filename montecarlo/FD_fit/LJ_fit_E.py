import numpy as np
import matplotlib.pylab as plt
#import scipy.optimize as sc
import itertools as it
import ctypes as ct
import sys
import glob
plt.ion()

#m = 1.043916 * 100
m = 1
h_bar = 6.582119
qo = 6
po = 2.067
D = 34.32*(h_bar/(po*qo))**3
rcut = np.sqrt(10)
pauli = [D, qo, po, np.sqrt(10)]
h = h_bar*2*np.pi
N = 10**3
Nsamp = 200

tipo = "rep"
rho = 'rho0'
Nbins = 100
Nreps = 8

nargs = len(sys.argv)
if (nargs >= 2):
  tipo = sys.argv[1]
if (nargs >= 3):
  rho = sys.argv[2]
if (nargs >= 4):
  Nbins = int(sys.argv[3])

if (rho[-1] == "0"):
  L = 1.357*10
if (rho[-1] == "1"):
  L = 0.983*10
if (rho[-1] == "2"):
  L = L/2
V = L**3

if (tipo=="f&v"):
  data = np.loadtxt("LUT_F32.txt")
  z = data[0,:]
  F32 = data[1,:]

  def mu(Y, T):
    if (30<=Y):
      return Ef*(1-np.pi**2/12*(T/Ef)**2)
    else:
      inf = 0
      sup = len(z)-1
      while (inf<sup-1):
        med = (inf+sup)//2
        if(F32[med]<Y):
          inf = med
        else:
          sup = med
      m = (F32[inf+1]-F32[inf])/(z[inf+1]-z[inf])
      return np.log(z[inf]+(Y-F32[inf])/m)*T

  pressure = ct.CDLL('../pressure.so')

  pressure_c = pressure.pressure_lj_PBC
  pressure_c.argtypes = [ct.c_voidp, ct.c_voidp, ct.c_longlong, ct.c_float,
                       ct.c_float, ct.c_float, ct.c_voidp, ct.c_float]
  pressure_c.restype = ct.c_float


  pairs = np.array(list(it.combinations(range(N), 2)), dtype=np.int64)
  def pressure(x):
    energ = 0
    force = np.zeros_like(x, dtype=np.float32)
    xp = x.ctypes.data_as(ct.c_voidp)
    pairsp = pairs.ctypes.data_as(ct.c_voidp)
    forcep = force.ctypes.data_as(ct.c_voidp)
    qF = pressure_c(xp, pairsp, len(pairs), 1.0, 1.0, rcut, forcep, L)
    return force, qF

  lennard = ct.CDLL('../../pexmd/interaction/lj.so')

  lennardforces_c = lennard.forces_PBC
  lennardforces_c.argtypes = [ct.c_voidp, ct.c_voidp, ct.c_longlong, ct.c_float,
                       ct.c_float, ct.c_float, ct.c_voidp, ct.c_float]
  lennardforces_c.restype = ct.c_float


  pairs = np.array(list(it.combinations(range(N), 2)), dtype=np.int64)
  def forces(x):
    energ = 0
    force = np.zeros_like(x, dtype=np.float32)
    xp = x.ctypes.data_as(ct.c_voidp)
    pairsp = pairs.ctypes.data_as(ct.c_voidp)
    forcep = force.ctypes.data_as(ct.c_voidp)
    qF = lennardforces_c(xp, pairsp, len(pairs), 1, 1, rcut, forcep, L)
    return force, qF

  def LJ_presion(q,p):
    P = 0
    dq = np.zeros(3)
    for i in range(1,N):
      for j in range(i):
        rij2 = 0
        dq = q[3*i:3*i+3]-q[3*j:3*j+3]
        for k in range(3):
          dq[k] = dq[k] + L*((dq[k] <- 0.5*L) - (0.5*L < dq[k]))
          rij2 += dq[k]**2
        rij6 = rij2**3
        P += np.sum(dq**2)*(2/rij6 - 1)/(rij2*rij6)
    P = 8*P/V
    P += np.sum(p**2)/m
    return P

  files = glob.glob("../LJ_fit/"+rho+"/distribucion_10_rep1_*")
  Ts = [f.split("_")[4][:-4] for f in files]
  n_temps = len(Ts)

  for k in range(n_temps):
    if Ts[k][2]== "E":
      Ts[k] = 10**float(Ts[k][3:])
    else:
      Ts[k] = float(Ts[k])
  idxs = np.argsort(Ts)
  Ts = np.array([Ts[i] for i in idxs])
  files = [files[i] for i in idxs]

  ns = np.zeros((n_temps, Nbins))
  bins = np.zeros((n_temps, Nbins+1))
  ns_q = np.zeros((n_temps, Nbins))
  bins_q = np.zeros((n_temps, Nbins+1))
  E = np.zeros((n_temps, Nbins))
  P = np.zeros(n_temps)
  qp = np.zeros(n_temps)
  pp = np.zeros(n_temps)
  qF = np.zeros(n_temps)
  mus = np.zeros(n_temps)
  for k in range(n_temps):
    data_q = np.array([], dtype=np.float32)
    data_p = np.array([], dtype=np.float32)
    for j in range(Nreps):
      data_aux = np.loadtxt("../LJ_fit/"+rho+"/distribucion_10_rep%d_%f.txt" %(j+1, Ts[k]), dtype=np.float32)
      data_q = np.concatenate([data_q, data_aux[:, 0]])
      data_p = np.concatenate([data_p, data_aux[:, 1]])
    print(k+1)
    for i in range(Nsamp*Nreps):
      fuerzas, qF_aux = pressure(data_q[3*N*i:3*N*(i+1)])
      pp[k] += np.sum(data_p[3*N*i:3*N*(i+1)]*data_p[3*N*i:3*N*(i+1)]/m)/(3*V*Nsamp*Nreps)
      qF[k] += qF_aux/(3*V*Nsamp*Nreps)
    P[k] = pp[k] + qF[k]
    data = np.zeros(len(data_p)//3)
    for i in range(len(data_p)//3):
      data[i] = np.sum(data_p[3*i:3*i+3]**2)/(2*m)
    ns[k, :], bins[k, :] = np.histogram(data, Nbins)
    ns_q[k, :], bins_q[k, :] = np.histogram(data_q, Nbins)
    ns[k, :] = ns[k, :]/(bins[k, 1]-bins[k, 0])
    ns_q[k, :] = ns_q[k, :]/(bins_q[k, 1]-bins_q[k, 0])

  ns = ns/(Nreps*Nsamp)
  ns_q = ns_q/(Nreps*3*Nsamp)
  N_max = np.max(ns)

  MB = lambda x, T: N*2*np.sqrt(x/np.pi)*np.exp(-x/T)/(T**1.5)
  FD = lambda x, mu, T: deg(x)/(np.exp((x-mu)/T)+1)

  seleccion = [0, 2, 4, 6, 8, 10, 12, 14, 15]
  i = 0
  plt.figure()
  for k in seleccion:
    i += 1
    E[k,:] = (bins[k,1:] + bins[k,:-1])/2
    if (n_temps > 1):
      plt.subplot(3, 3, i)
    plt.xlabel(r"$E$")
    plt.plot(E[k,:],  ns[k,:], "ko")
    plt.text((min(E[k,:])+ 5*max(E[k,:]))/6, 0.9*max(ns[k,:]), "T=%f" %(Ts[k]))
    rango = np.linspace(0, E[k,-1], 1000)
    exacto_MB = MB(rango, Ts[k])
    plt.plot(rango, exacto_MB, "r-")
  """
  plt.figure()
  for k in range(n_temps):
    q[k,:] = (bins_q[k,1:] + bins_q[k,:-1])/2
    if (n_temps > 1):
      plt.subplot(n_temps/3, 3, k+1)
    plt.xlabel(r"$q$")
    plt.plot(q[k,:],  ns_q[k,:], "ko")
    plt.axis([0, 60, 0, 100])
    plt.text((min(q[k,:])+ 5*max(q[k,:]))/6, 1.25*max(ns_q[k,:]), "T=%f" %(Ts[k]))
  """
  np.savetxt("../LJ_fit/"+rho+"/presiones_N=1000.txt", [Ts,P])
  plt.figure()
  plt.plot(Ts, P*V/N, "ro-")
  plt.plot(Ts, Ts, "k--")
  plt.legend([r'LJ', 'Boltzmann'], loc=9)
  plt.xlabel(r"$T$")
  plt.ylabel(r"$P/\rho$")
  plt.show()

if (tipo=="p"):
  rhos = ['rho0', 'rho1']
  V = np.array([1.357*10,0.983*10])**3
  i = 0
  plt.figure()
  for rho in rhos:
    data = np.loadtxt("../LJ_fit/"+rho+'/presiones_N=1000.txt')
    Ts = data[0,:]
    P = data[1,:]
    plt.plot(Ts, P*V[i]/N, "o--")
    i += 1
  plt.plot(Ts, Ts, "-")
  plt.xlabel(r"$T$")
  plt.ylabel(r"$P/\rho$")
  plt.legend([r"$\rho^* = 0.4$", r"$\rho^* = 0.95$", "Boltzmann"], loc=9)
  plt.show()
