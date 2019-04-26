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
Nreps = 8
Nsamp = 100
N = 8000

nargs = len(sys.argv)
if (nargs >= 2):
  tipo = sys.argv[1]

if (tipo == "dist"):

  T = 10
  rho = 0.2
  D = 100
  if (nargs >= 3):
    T = float(sys.argv[2])
  if (nargs >= 4):
    rho = float(sys.argv[3])
  if (nargs >= 5):
    D = float(sys.argv[4])

  h_bar = 1
  m = 1
  deg = lambda E: 4*np.pi*np.sqrt(2*m**3*E)*(L/h)**3
  long_term = lambda T: (2*np.pi*h_bar**2/(m*T))**0.5
  Ef = h_bar**2/(2*m)*(6*np.pi**2*rho)**(2/3)

  data = np.loadtxt("FD_fit/LUT_F32.txt")
  z = data[0,:]
  F32 = data[1,:]

  F32_A = lambda z: (4/(3*np.pi**0.5))*(np.log(z)**1.5)*(1 + np.pi**2/(8*np.log(z)**2) + (7*np.pi**4/640)*np.log(z)**(-4))

  def dame_z(Y):
    F32_A = lambda z: (4/(3*np.pi**0.5))*(np.log(z)**1.5)*(1 + np.pi**2/(8*np.log(z)**2) + (7*np.pi**4/640)*np.log(z)**(-4))
    if (30<=Y):
      z_inf = 30.0
      z_sup = 2*z_inf
      while(F32_A(z_sup) < Y):
        z_sup *= 2
      z_med = (z_inf+z_sup)/2
      print(Y, z_sup)
      while (Y*1E-4 < np.abs(F32_A(z_med)-Y)):
        z_med = (z_inf+z_sup)/2
        if(F32_A(z_med)<Y):
          z_inf = z_med
        else:
          z_sup = z_med
      return z_med
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
      return z[inf]+(Y-F32[inf])/m

  def F_32(fug):
    if (z[-1]<=fug):
      return F32_A(fug)
    else:
      inf = 0
      sup = len(z)-1
      while (inf<sup-1):
        med = (inf+sup)//2
        if(z[med]<fug):
          inf = med
        else:
          sup = med
      m = (F32[inf+1]-F32[inf])/(z[inf+1]-z[inf])
      return F32[inf] + m*(fug-z[inf])


  data = np.loadtxt("prueba_pauli/dist_%.2f_%.2f_%.2f.txt" %(rho, D, T))
  ns, E_bins = np.histogram(data, Nbins)
  ns = ns/(Nsamp*(E_bins[1]-E_bins[0]))

  MB = lambda E, T: N*np.sqrt(4*E/(np.pi*T**3))*np.exp(-E/T)
  FD = lambda x, mu, T: deg(x)/(np.exp((x-mu)/T)+1)

  Es = (E_bins[1:]+E_bins[:-1])/2
  plt.figure()
  plt.plot(Es/T, ns, "k.")
  plt.plot(Es/T, MB(Es, T), "r-")
  plt.xlabel(r"$\varepsilon_{cin}$", fontsize=20)
  plt.ylabel(r"$f(\varepsilon_{cin})$", fontsize=20)
  plt.legend(["Pauli", "Boltzmann"])
  plt.title("$D^* = %.1f$\t$T^*=%.1f$" %(D, T))
  #plt.plot(Es, FD(Es, np.log(dame_z(long_term(Ts[k])**3*rho))*T, T), "b-")


if (tipo == "comp"):

  Teff = 10
  rho = 0.2
  D = 100
  if (nargs >= 3):
    Teff = float(sys.argv[2])
  if (nargs >= 4):
    rho = float(sys.argv[3])
  if (nargs >= 5):
    D = float(sys.argv[4])
  leyenda = []
  plt.figure()
  for d in range(-1,3):
    D = 10**(2+d)
    T = Teff*0.5**-d
    data = np.loadtxt("prueba_pauli/dist_%.2f_%.2f_%.2f.txt" %(rho, D, T))
    ns, E_bins = np.histogram(data, Nbins)
    ns = ns/(Nsamp*(E_bins[1]-E_bins[0]))

    MB = lambda E, T: N*np.sqrt(4*E/(np.pi*T**3))*np.exp(-E/T)
    FD = lambda x, mu, T: deg(x)/(np.exp((x-mu)/T)+1)

    Es = (E_bins[1:]+E_bins[:-1])/2
    plt.plot(Es/T, ns*T, ".")
    leyenda.append("$D^* = %.1f$\t$T^*=%.1f$" %(D, T))
  plt.plot(Es/T, MB(Es,T)*T, "r-")
  leyenda.append("Boltzmann")
  plt.legend(leyenda)
  plt.xlabel(r"$\varepsilon_{cin}/T$", fontsize=14)
  plt.ylabel(r"$f(\varepsilon_{cin})*T$", fontsize=14)


if (tipo == "Ekin"):
  rho = 0.2
  D = 100
  if (nargs >= 3):
    rho = float(sys.argv[2])
  if (nargs >= 4):
    D = float(sys.argv[3])
  leyenda = []
  plt.figure()
  for D in [10, 100, 1000]:
    data = np.loadtxt("prueba_pauli/termo_pauli_%.2f_%.2f.txt" %(rho, D))
    plt.plot(data[:,1], data[:,2]/data[:,1], '.--')
    leyenda.append("$D^*=%.1f$" %(D))
  plt.plot(data[:,1], 1.5*data[:,1]/data[:,1], "r-")
  plt.plot(data[:,1], data[40,2]/data[40,1]*data[:,1]/data[:,1], 'k--')
  leyenda.append("Boltzmann")
  plt.legend(leyenda, loc = 'upper center')
  plt.xlabel("T", fontsize=14)
  plt.ylabel(r"$E_{cin}$", fontsize=14)
