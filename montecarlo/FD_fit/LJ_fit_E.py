import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
import scipy.ndimage
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
Nsamp = 1 #200

tipo = "rep"
rho = 'rho0'
Nbins = 100
Nreps = 10

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
  L = 0.7*10
V = L**3

pressure = ct.CDLL('../pressure.so')

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

if (tipo == "a"):

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

  P = np.zeros(n_temps)
  qp = np.zeros(n_temps)
  pp = np.zeros(n_temps)
  qF = np.zeros(n_temps)

  for k in range(n_temps):
    data_q = np.array([], dtype=np.float32)
    data_p = np.array([], dtype=np.float32)
    for j in range(Nreps):
      data_aux = np.loadtxt("../LJ_fit/"+rho+"/distribucion_10_rep%d_%f.txt" %(j+1, Ts[k]), dtype=np.float32)
      data_q = np.concatenate([data_q, data_aux[:, 0]])
      data_p = np.concatenate([data_p, data_aux[:, 1]])
    """
    for i in range(Nsamp*Nreps):
      fuerzas, qF_aux = pressure(data_q[3*N*i:3*N*(i+1)])
      pp[k] += np.sum(data_p[3*N*i:3*N*(i+1)]*data_p[3*N*i:3*N*(i+1)]/m)/(3*V*Nsamp*Nreps)
      qF[k] += qF_aux/(3*V*Nsamp*Nreps)
    P[k] = pp[k] + qF[k]
    """
    data = np.zeros(len(data_p)//3)
    for i in range(len(data_p)//3):
      data[i] = np.sum(data_p[3*i:3*i+3]**2)/(2*m)
    ns, bins = np.histogram(data, Nbins)
    ns_q, bins_q = np.histogram(data_q, Nbins)
    ns = ns/((bins[1]-bins[0])*Nreps*Nsamp)
    ns_q = ns_q/((bins_q[1]-bins_q[0])*3*Nreps*Nsamp)
    E = (bins[1:] + bins[:-1])/2
    filename = "../LJ_fit/"+rho+"/histograma_%d_T=%f.txt" %(Nbins, Ts[k])
    print(filename)
    np.savetxt(filename, [E, ns])

  #np.savetxt('../LJ_fit/'+rho+"/presiones_N=1000.txt", [Ts,P])


if (tipo == 'f&v'):

  MB = lambda x, T: N*2*np.sqrt(x/np.pi)*np.exp(-x/T)/(T**1.5)

  #seleccion = [0, 2, 4, 6, 8, 10, 12, 14, 15]
  seleccion = range(n_temps)
  i = 0
  plt.figure()
  for k in seleccion:
    data = np.loadtxt('../LJ_fit/'+rho+"/histograma_%d_T=%f.txt" %(Nbins, Ts[k]))
    E = data[0,:]
    ns = data[1,:]
    if (n_temps > 1):
      plt.subplot(3, 3, i+1)
    plt.xlabel(r"$E$")
    plt.plot(E, ns, "ko")
    plt.text((min(E)+ 5*max(E))/6, 0.9*max(ns), "T=%f" %(Ts[k]))
    #plt.text(3, 700, "T=%f" %(Ts[k]))
    rango = np.linspace(0, E[-1], 1000)
    exacto_MB = MB(rango, Ts[k])
    plt.plot(rango, exacto_MB, "r-")
    #plt.axis([0, 5, 0, 1000])
    i += 1

  plt.show()


if (tipo == "p"):
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


if(tipo == 'e'):

  deltas_c = pressure.delta_fases
  deltas_c.argtypes = [ct.c_voidp, ct.c_voidp, ct.c_voidp, ct.c_longlong,
                        ct.c_voidp, ct.c_voidp, ct.c_float]
  deltas_c.restype = ct.c_float

  pairs = np.array(list(it.combinations(range(N), 2)), dtype=np.int64)
  p = len(pairs)
  def deltas(x, p):
    dq = np.zeros(len(pairs), dtype=np.float32)
    dp = np.zeros(len(pairs), dtype=np.float32)
    xp = x.ctypes.data_as(ct.c_voidp)
    pp = p.ctypes.data_as(ct.c_voidp)
    pairsp = pairs.ctypes.data_as(ct.c_voidp)
    dq_p = dq.ctypes.data_as(ct.c_voidp)
    dp_p = dp.ctypes.data_as(ct.c_voidp)
    deltas_c(xp, pp, pairsp, len(pairs), dq_p, dp_p, L)
    return dq, dp

  dq = np.zeros(p*Nsamp*Nreps)
  dp = np.zeros(p*Nsamp*Nreps)
  for k in range(1):
    data_q = np.array([], dtype=np.float32)
    data_p = np.array([], dtype=np.float32)
    for j in range(Nreps):
      data_aux = np.loadtxt('../LJ_fit/'+rho+"/distribucion_10_rep%d_%f.txt" %(j+1, Ts[k]), dtype=np.float32)
      data_q = np.concatenate([data_q, data_aux[:, 0]])
      data_p = np.concatenate([data_p, data_aux[:, 1]])
    plt.figure()
    for i in range(Nsamp*Nreps):
      dq[p*i:p*(i+1)], dp[p*i:p*(i+1)] = deltas(data_q[3*N*i:3*N*(i+1)], data_p[3*N*i:3*N*(i+1)])
    counts,xbins,ybins,image = plt.hist2d(dq,dp,bins=100, norm=LogNorm(), cmap = plt.cm.rainbow)
    plt.colorbar()
    new_counts = scipy.ndimage.filters.gaussian_filter(counts, 1)
    CS = plt.contour(new_counts.transpose(),extent=[xbins[0],xbins[-1],ybins[0],ybins[-1]],
                  linewidths=3, colors = "black", levels = np.logspace(1, 3, 3))
    plt.clabel(CS, colors = "black", inline=True, fmt="%d", fontsize=20)
    plt.axis([0, xbins[-1], 0, ybins[-1]])
    plt.xlabel(r'$\Delta q$')
    plt.ylabel(r'$\Delta p$')
    plt.title(r'$T=%1.4fMeV$' %(Ts[k]))
  plt.show()


if (tipo == "gr"):

  gr_c = pressure.Gr
  gr_c.argtypes = [ct.c_voidp, ct.c_int, ct.c_float,
                       ct.c_float, ct.c_voidp, ct.c_int]
  gr_c.restype = ct.c_float


  pairs = np.array(list(it.combinations(range(N), 2)), dtype=np.int64)
  def Gr(x, dr):
    M = int(np.ceil(np.sqrt(3)*L/(2*dr)))
    gr = np.zeros(M, dtype=np.float32)
    xp = x.ctypes.data_as(ct.c_voidp)
    grp = gr.ctypes.data_as(ct.c_voidp)
    gr_c(xp, len(x)//3, dr, L, grp, M)
    return gr

  seleccion = range(n_temps)
  dr = 0.025

  for k in seleccion:
    data_q = np.array([], dtype=np.float32)
    data_p = np.array([], dtype=np.float32)
    for j in range(Nreps):
      data_aux = np.loadtxt('../LJ_fit/'+rho+"/distribucion_10_rep%d_%f.txt" %(j+1, Ts[k]), dtype=np.float32)
      data_q = np.concatenate([data_q, data_aux[:, 0]])
      data_p = np.concatenate([data_p, data_aux[:, 1]])
    plt.figure()
    gr = Gr(data_q[0:3*N], dr)/(Nsamp*Nreps)
    for i in range(1,Nsamp*Nreps):
      gr += Gr(data_q[3*N*i:3*N*(i+1)], dr)/(Nsamp*Nreps)
    plt.plot(dr*np.arange(len(gr)), gr, "b-")
    plt.plot([0, 0.5*L], [1, 1], "r-")
    plt.xlabel(r"$r$")
    plt.ylabel(r"$g(r)$")
    plt.title(r"$T=%1.4f\epsilon$    $\rho=%1.2f\sigma^{-3}$" %(Ts[k], N/V))
    plt.axis([0, 0.5*L, 0, max(gr)])
    plt.show()

if (tipo == "lin"):
  seleccion = [0, 2, 4, 6, 8, 10, 12, 14, 15]
  lindemann = np.zeros(n_temps)
  for k in range(n_temps):
    for j in range(Nreps):
      data_aux = np.loadtxt('../LJ_fit/'+rho+"/distribucion_10_rep%d_%f.txt" %(j+1, Ts[k]), dtype=np.float32)
      data_q = data_aux[:, 0]
      means = [np.mean(data_q[np.arange(200)*3*N+i]) for i in range(3*N)]
      vars = [np.linalg.norm(data_q[np.arange(200)*3*N+i]-means[i])**2 for i in range(3*N)]
      lindemann[k] += np.sqrt(np.mean(vars))*(N/V)**(1/3)/Nreps
  plt.plot(Ts, lindemann, "bo--")
  plt.show()
