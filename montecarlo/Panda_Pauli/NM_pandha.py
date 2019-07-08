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
    Epauli = data[:,3]
    P = data[:,-4]
    Tv = data[:,-3]

    plt.figure(1)
    plt.plot(Ts, Ekin+Enuc+Epauli, ".-")
    plt.figure(2)
    plt.plot(Ts, Tv, ".-")
    plt.figure(3)
    plt.plot(Ts, P, ".-")
    leyenda.append(r"$\rho = %.3f$fm$^{-3}$" %rho)

  plt.figure(1)
  plt.xlabel("$T$ [MeV]")
  plt.ylabel("Energia [MeV]")
  plt.legend(leyenda, loc = "upper left")
  plt.figure(2)
  plt.xlabel("$T$ [MeV]")
  plt.ylabel("Temperatura efectiva [MeV]")
  plt.legend(leyenda, loc = "upper left")
  plt.plot(Ts, Ts, "r-")
  plt.figure(3)
  plt.xlabel("$T$ [MeV]")
  plt.ylabel(r"Presion [MeV fm$^{-3}$]")
  plt.legend(leyenda, loc = "upper left")
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

if (tipo == "sc"):
  if (nargs >= 3):
    folder = sys.argv[2]
  if (nargs >= 4):
    rhos = np.zeros(nargs-3)
    for i in range(nargs-3):
      rhos[i] = float(sys.argv[i+3])
  else:
    rhos = [0.16]

  rhos = np.array(range(10, 25))/100
  E = []
  for rho in rhos:
    data = np.loadtxt("SC/data_%.3f.txt" %rho)
    E.append(data[0,0])
  E = np.array(E)
  plt.plot(rhos, E, "o-")
  params = np.polyfit(rhos, E, 2)
  print(np.polyval(params, 0.16), np.polyval([2*params[0], params[1]], 0.16), 9*(0.16**2)*params[0]*2)
  X = np.linspace(0.08, 0.25, 15001)
  plt.plot(X, np.polyval(params, X), "-")

  data = np.loadtxt("SC/data_0.160.txt")
  E_160_real = data[0, 0]
  E_160 = data[1:, 1:]
  As = data[1:,0]
  Bs = data[0,1:]

  Es = []
  for rho in rhos:
    data = np.loadtxt("SC/data_%.3f.txt" %rho)
    Es.append(data[1:,1:])
  Es = np.array(Es)
  plt.figure()
  plt.plot(rhos, E, "o-")
  leyendas = ["Pandharipande"]
  idx_0160 = np.where(rhos==0.16)[0]
  idxs = np.where(np.abs(Es[idx_0160,:,:]-E[idx_0160])<10)
  print(idxs)
  for i in np.unique(idxs[1]):
    for j in np.unique(idxs[2]):
      plt.plot(rhos, Es[:,i,j], ".--")
#      leyendas.append("a=%.2f | b=%.2f" %(As[i], Bs[j]))
  plt.legend(leyendas, loc = "upper left")

if (tipo == "sc2"):
  if (nargs >= 3):
    folder = sys.argv[2]
  else:
    rhos = [0.16]

  data = np.loadtxt(folder+"/Evsrho_temp.txt")
  rhos = data[:, 0]
  E = data[:,1]
  Ekin = data[:,2]
  idxs = np.arange(len(rhos)//2+1, 11, 1)
  idxs = np.arange(len(rhos))
  plt.figure()
  plt.plot(rhos, E, "o-")
  rhos = rhos[idxs]
  E = E[idxs]
  Ekin = Ekin[idxs]
  params = np.polyfit(rhos, E, 2)
  print("  Total :", params)
  X = np.linspace(min(rhos)*0.9, max(rhos)*1.1, 15001)
  plt.plot(X, np.polyval(params, X), "-")
  print("%.12f" %rhos[np.where(E==min(E))[0]], min(E))
  rho_min = -params[1]/(2*params[0])
  print(rho_min, np.polyval(params, rho_min))
  print("Derivadas", np.polyval([2*params[0], params[1]], 0.16), 9*(0.16**2)*params[0]*2)
  rhos = data[:, 0]
  E = data[:,1]
  Ekin = data[:,2]
  idxs = np.arange(0, len(rhos)//2-1, 1)
  rhos = rhos[idxs]
  E = E[idxs]
  Ekin = Ekin[idxs]
  params = np.polyfit(rhos, E, 2)
  print("  Total :", params)
  X = np.linspace(min(rhos)*0.9, max(rhos)*1.1, 15001)
  plt.plot(X, np.polyval(params, X), "-")
  print("%.12f" %rhos[np.where(E==min(E))[0]], min(E))
  rho_min = -params[1]/(2*params[0])
  print(rho_min, np.polyval(params, rho_min))
  print("Derivadas", np.polyval([2*params[0], params[1]], 0.16), 9*(0.16**2)*params[0]*2)

  """
  plt.figure()
  plt.plot(rhos, E-Ekin, "o-")
  params = np.polyfit(rhos, E-Ekin, 2)
  print("Cinetica:", params)
  plt.plot(X, np.polyval(params, X), "-")
  plt.figure()
  plt.plot(rhos, Ekin, "o-")
  params = np.polyfit(rhos, Ekin, 2)
  print("Potencial:", params)
  plt.plot(X, np.polyval(params, X), "-")
  """
