import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
import scipy.ndimage
import itertools as it
import ctypes as ct
import sys
import glob
plt.ion()

m = 1
rcut = 2.5
N = 4000

tipo = "rep"

nargs = len(sys.argv)
if (nargs >= 2):
  tipo = sys.argv[1]
if (nargs >= 3):
  rho = sys.argv[2]
if (nargs >= 4):
  Nbins = int(sys.argv[3])


pressure = ct.CDLL('../pressure.so')

files = glob.glob("termo_lennard_*")
rhos = [float(f.split("_")[2][:-10]) for f in files]
n_rhos = len(rhos)

idxs = np.argsort(rhos)
rhos = np.array([rhos[i] for i in idxs])
files = [files[i] for i in idxs]

Ts = np.loadtxt(files[0])[:,1]
n_temps = len(Ts)

if (tipo == "t"):
  Ecin = np.zeros((n_rhos, n_temps))
  Epot = np.zeros((n_rhos, n_temps))
  Etot = np.zeros((n_rhos, n_temps))
  i = 0
  plt.figure()
  for f in files:
    data = np.loadtxt(f)
    Ecin[i,:] = data[:,2]
    Epot[i,:] = data[:,3]
    Etot[i,:] = data[:,4]
    plt.plot(Ts, Etot[i,:])
    i+=1
  plt.show()


if (tipo == "gr"):

  gr_c = pressure.Gr
  gr_c.argtypes = [ct.c_voidp, ct.c_int, ct.c_float,
                       ct.c_float, ct.c_voidp, ct.c_int]
  gr_c.restype = ct.c_float


  def Gr(x, dr, L):
    M = int(np.ceil(np.sqrt(3)*L/(2*dr)))
    gr = np.zeros(M, dtype=np.float32)
    xp = x.ctypes.data_as(ct.c_voidp)
    grp = gr.ctypes.data_as(ct.c_voidp)
    gr_c(xp, len(x)//3, dr, L, grp, M)
    return gr

  def cargar_lammpstrj(filename, T_idx):
    f = open(filename, "r")
    for j in range(T_idx):
      for k in range(3):
        f.readline()
      N = int(f.readline())
      for k in range(5):
        f.readline()
      q = np.zeros(3*N, dtype=np.float32)
      p = np.zeros(3*N, dtype=np.float32)
      for i in range(N):
        data = f.readline()
        data = data.split(" ")
        for k in range(3):
          q[3*i+k] = float(data[2+k])
          p[3*i+k] = float(data[5+k])
    return q, p

  dr = 0.025
  T_idx = len(Ts)-1

  for rho in [0.9]:
    q, p = cargar_lammpstrj("configuracion_lennard_%1.2f" %rho + ".lammpstrj", T_idx)
    L = (len(q)/(3*rho))**(1./3)
    print(len(q), len(q)//3, L)

    gr = Gr(q, dr, L)
    plt.figure()
    plt.plot(dr*np.arange(len(gr)), gr, "b-")
    plt.plot([0, 0.5*L], [1, 1], "r-")
    plt.xlabel(r"$r$")
    plt.ylabel(r"$g(r)$")
    plt.title(r"$T=%1.4f\epsilon$    $\rho=%1.2f\sigma^{-3}$" %(Ts[T_idx], rho))
    plt.axis([0, 0.5*L, 0, max(gr)])
    plt.show()
