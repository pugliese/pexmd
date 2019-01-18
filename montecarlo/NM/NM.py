import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
import scipy.ndimage
import itertools as it
import ctypes as ct
import sys
import glob
plt.ion()

m = 1.043916 * 100
h_bar =  6.582119
qo = 6
po = 2.067
D = 34.32*(h_bar/(po*qo))**3
scut = np.sqrt(10)
pauli = [D, qo, po, np.sqrt(10)]
h = h_bar*2*np.pi
N = 10**3
Nsamp = 1

tipo = "rep"
caso = 'layers/x1/'
Nbins = 100
Nreps = 1

nargs = len(sys.argv)
if (nargs >= 2):
  tipo = sys.argv[1]
if (nargs >= 3):
  caso = sys.argv[2]
if (nargs >= 4):
  Nbins = int(sys.argv[3])




files = glob.glob(caso+"energias_*")
rhos = [float(f.split("_")[1][:-4]) for f in files]
n_rhos = len(rhos)
rhos = np.sort(rhos)

if(tipo == 't'):
  E_kin = []
  E_nuc = []
  E_pauli = []
  E = []
  Es = np.zeros(n_rhos)

  for k in range(n_rhos):
    data = np.loadtxt(caso+"energias_%f.txt" %(rhos[k]), dtype=np.float32)
    E_kin.append(data[:,0])
    E_nuc.append(data[:,1])
    E_pauli.append(data[:,2])
    E.append(E_kin[-1] + E_nuc[-1] + E_pauli[-1])
    plt.subplot(n_rhos//2, 2, k+1)
    plt.plot(E[-1]/N, "b-")
    """
    plt.plot(E_kin[-1]/N, "r--")
    plt.plot(E_nuc[-1]/N, "k--")
    plt.plot(E_pauli[-1]/N, "g--")
    """
    plt.legend([r"$\rho=%f fm^{-3}$" %rhos[k]])
    Es[k] = np.mean(E[-1])/N
  plt.figure()
  plt.plot(rhos, Es, "o--")
  plt.show()

pressure = ct.CDLL('../pressure.so')

if(tipo == 'e'):

  deltas_c = pressure.delta_fases
  deltas_c.argtypes = [ct.c_voidp, ct.c_voidp, ct.c_voidp, ct.c_longlong,
                        ct.c_voidp, ct.c_voidp, ct.c_float]
  deltas_c.restype = ct.c_float

  n = N//4
  pairs = np.array(list(it.combinations(range(n), 2)), dtype=np.int64)
  p = len(pairs)
  def deltas(x, p, L):
    dq = np.zeros(len(pairs), dtype=np.float32)
    dp = np.zeros(len(pairs), dtype=np.float32)
    xp = x.ctypes.data_as(ct.c_voidp)
    pp = p.ctypes.data_as(ct.c_voidp)
    pairsp = pairs.ctypes.data_as(ct.c_voidp)
    dq_p = dq.ctypes.data_as(ct.c_voidp)
    dp_p = dp.ctypes.data_as(ct.c_voidp)
    deltas_c(xp, pp, pairsp, len(pairs), dq_p, dp_p, L)
    return dq, dp

  L = (N/rhos)**(1/3)

  dq = np.zeros(4*len(pairs))
  dp = np.zeros(4*len(pairs))

  for k in range(1):
    data_aux = np.loadtxt("distribucion_%f.txt" %(rhos[k]), dtype=np.float32)
    data_q = data_aux[:, 0]
    data_p = data_aux[:, 1]
    plt.figure()
    for i in range(4):
      dq[i*p:(i+1)*p], dp[i*p:(i+1)*p] = deltas(data_q[i*3*n:(i+1)*3*n], data_p[i*3*n:(i+1)*3*n], L[k])
    counts, xbins, ybins, image = plt.hist2d(dq, dp, bins=100, norm=LogNorm(), cmap = plt.cm.rainbow)
    plt.colorbar()
    new_counts = scipy.ndimage.filters.gaussian_filter(counts, 1)
    #plt.subplot(2, 2, i+1)
    CS = plt.contour(new_counts.transpose(),extent=[xbins[0],xbins[-1],ybins[0],ybins[-1]],
                  linewidths=3, colors = "black", levels = np.logspace(1, 3, 3))
    plt.clabel(CS, colors = "black", inline=True, fmt="%d", fontsize=20)
    plt.xlabel(r'$\Delta q$')
    plt.ylabel(r'$\Delta p$')
    plt.title(r'$\rho=%1.4f fm^{-3}$' %(rhos[k]))
  plt.show()


if (tipo == "gr"):

  gr_c = pressure.Gr
  gr_c.argtypes = [ct.c_voidp, ct.c_int, ct.c_float,
                       ct.c_float, ct.c_voidp, ct.c_int]
  gr_c.restype = ct.c_float


  pairs = np.array(list(it.combinations(range(N), 2)), dtype=np.int64)
  def Gr(x, dr, L):
    M = int(np.ceil(np.sqrt(3)*L/(2*dr)))
    gr = np.zeros(M, dtype=np.float32)
    xp = x.ctypes.data_as(ct.c_voidp)
    grp = gr.ctypes.data_as(ct.c_voidp)
    gr_c(xp, len(x)//3, dr, L, grp, M)
    return gr

  L = (N/rhos)**(1/3)
  dr = L/(2*400)
  gr = []
  i = 0
  for k in range(1):
    data_aux = np.loadtxt("distribucion_%f_10.txt" %(rhos[k]), dtype=np.float32)
    data_q = data_aux[:, 0]
    gr.append(Gr(data_q[0:3*N], dr[i], L[i])/(Nsamp*Nreps))
  for i in range(1):
    plt.subplot(3, 2, i+1)
    plt.plot(dr[i]*np.arange(len(gr[i])), gr[i], "b-")
    plt.plot([0, 0.5*L[i]], [1, 1], "r-")
    plt.xlabel(r"$r$")
    plt.ylabel(r"$g(r)$")
    plt.title(r"$T=0.5 MeV$ $\rho=%f fm^{-3}$" %(rhos[i]))
    plt.axis([0, 0.5*L[i], 0, max(gr[i])])
  plt.show()
