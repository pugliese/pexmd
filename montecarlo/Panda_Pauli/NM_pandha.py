import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
import scipy.ndimage
import itertools as it
import ctypes as ct
import sys
import glob
#import scipy.optimize as sc
plt.rcParams['axes.color_cycle'] = ["#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f","#e5c494"]
plt.ion()


tipo = "dist"
Nbins = 200
Nsamp = 50
N = 8000
m = 938

nargs = len(sys.argv)
if (nargs >= 2):
  tipo = sys.argv[1]


if (tipo == "termo"):
  if (nargs >= 3):
    folder = sys.argv[2]
  if (nargs >= 4):
    rhos = np.zeros(nargs-3)
    for i in range(nargs-3):
      rhos[i] = float(sys.argv[i+3])
  else:
    rhos = [0.16]

  leyenda = []
  for rho in rhos:
    data = np.loadtxt(folder+"/termo_%.3f.txt" %(rho))
    Ts = data[:,0]
    Ekin = data[:,1]
    Enuc = data[:,2]
    P = data[:,-4]
    Tv = data[:,-3]

    plt.figure(1)
    plt.plot(Ts, Ekin+Enuc, "b.-")
    plt.figure(2)
    plt.plot(Ts, Ekin, "r.-")
    leyenda.append(r"$\rho = %.3f$fm$^{-3}$")

  plt.figure(1)
  plt.xlabel("$T$ [MeV]")
  plt.ylabel("Energia [MeV]")
  plt.legend(leyenda, loc = "upper center")
  plt.figure(2)
  plt.xlabel("$T$ [MeV]")
  plt.ylabel("Energia cinética [MeV]")
  plt.legend(leyenda, loc = "upper center")
  plt.plot(Ts, 1.5*Ts, "k--")
  plt.show()


if (tipo == "Evsrho"):
  def almost_equal(a, b, tol=1E-3):
    res = np.zeros_like(a)
    for i in range(len(res)):
      res[i] = (abs(a[i]-b)<tol)
    return res

  if (nargs >= 3):
    folder = sys.argv[2]
  if (nargs >= 4):
    Ts = np.zeros(nargs-3)
    for i in range(nargs-3):
      Ts[i] = float(sys.argv[i+3])
  else:
    Ts = [0.5]

  files = glob.glob(folder+"/termo_*")
  rhos = np.sort(np.array([float(file[-9:-4]) for file in files]))
  rhos = np.delete(rhos, np.where(rhos==0.14)[0])
  rhos = np.delete(rhos, np.where(rhos==0.15)[0])
  rhos = np.delete(rhos, np.where(rhos==0.18)[0])
  n_rhos = len(rhos)
  n_temps = len(Ts)
  leyenda = []
  Ekin = np.zeros((n_rhos, n_temps))
  Enuc = np.zeros((n_rhos, n_temps))
  i = 0
  for rho in rhos:
    data = np.loadtxt(folder+"/termo_%.3f.txt" %(rho))
    idxs = [np.where(almost_equal(data[:,0], T))[0][0] for T in Ts]
    Ekin[i,:] = data[:,1][idxs]
    Enuc[i,:] = data[:,2][idxs]
    i += 1

  plt.figure(1)
  for j in range(n_temps):
    plt.plot(rhos, Ekin[:,j]+Enuc[:,j], "o-")
  plt.xlabel(r"$\rho$ [fm$^{-3}$]")
  plt.ylabel("Energia [MeV]")
  plt.legend([r"$T=%.3f$MeV" %T for T in Ts], loc = "upper center")



if (tipo == "sat"):
  if (nargs >= 3):
    folder = sys.argv[2]
  if (nargs >= 4):
    rhos = np.zeros(nargs-3)
    for i in range(nargs-3):
      rhos[i] = float(sys.argv[i+3])
  else:
    rhos = [0.16]

  leyenda = []
  data = np.loadtxt(folder+"/termo_%.3f.txt" %(rhos[0]))
  Ts = data[:,0]
  Ts = Ts[np.where(Ts<1)[0]]
  n_temps = len(Ts)
  n_rhos = len(rhos)
  Ekin = np.zeros((n_rhos, n_temps))
  Enuc = np.zeros((n_rhos, n_temps))
  for i in range(n_rhos):
    data = np.loadtxt(folder+"/termo_%.3f.txt" %(rhos[i]))
    Ekin[i,:] = data[:,1][np.where(Ts<1)[0]]
    Enuc[i,:] = data[:,2][np.where(Ts<1)[0]]
  E = Ekin + Enuc

  K = np.zeros(n_temps)
  E_sat = np.zeros(n_temps)
  rho_sat = np.zeros(n_temps)
  for j in range(n_temps):
    idx_sat = np.where(E[:,j]==min(E[:,j]))[0][0]
    E_sat[j] = E[idx_sat, j]
    rho_sat[j] = rhos[idx_sat]
    params = np.polyfit(rhos, E[:,j], 2)
    K[j] = params[0]*(18*rho_sat[j]**2)

  plt.figure()
  plt.plot(Ts, E_sat, "b.-")
  plt.xlabel("$T$ [MeV]")
  plt.ylabel("Energia de saturacion [MeV]")
  plt.figure()
  plt.plot(Ts, K, "b.-")
  plt.xlabel("$T$ [MeV]")
  plt.ylabel("Compresibilidad [MeV]")
  plt.figure()
  plt.plot(Ts, rho_sat, "b.-")
  plt.xlabel("$T$ [MeV]")
  plt.ylabel("Densidad de saturacion [fm$^{-3}$]")
  plt.show()
