import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
import scipy.ndimage
import itertools as it
import ctypes as ct
import sys
import glob
import scipy.optimize as sc
plt.rcParams['axes.color_cycle'] = ["#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f","#e5c494"]
plt.ion()


tipo = "dist"
Nbins = 200
Nreps = 8
Nsamp = 50
N = 8000
params = 1

nargs = len(sys.argv)
if (nargs >= 2):
  tipo = sys.argv[1]


MB = lambda E: N*np.sqrt(4*E/(np.pi))*np.exp(-E)
FD = lambda E, a: N*np.sqrt(4*E/(np.pi))*np.exp(-E)/((a/dame_z(a))+a*np.exp(-E))

if (True):
  h_bar = 197.327
  m = 938
  deg = lambda E: 4*np.pi*np.sqrt(2*m**3*E)*(L/h)**3
  long_term = lambda T: (2*np.pi*h_bar**2/(m*T))**0.5

  data = np.loadtxt("../FD_fit/LUT_F32.txt")
  z32 = data[0,:]
  F32 = data[1,:]
  F32_A = lambda z: (4/(3*np.pi**0.5))*(np.log(z)**1.5)*(1 + np.pi**2/(8*np.log(z)**2) + (7*np.pi**4/640)*np.log(z)**(-4))

  data = np.loadtxt("../FD_fit/LUT_F52.txt")
  z52 = data[0,:]
  F52 = data[1,:]

  F52_A = lambda z: (8/(15*np.pi**0.5))*(np.log(z)**2.5)*(1 + 5*np.pi**2/(8*np.log(z)**2) - (5*3/16)*(7*np.pi**4/360)*np.log(z)**-4)

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
      sup = len(z32)-1
      while (inf<sup-1):
        med = (inf+sup)//2
        if(F32[med]<Y):
          inf = med
        else:
          sup = med
      m = (F32[inf+1]-F32[inf])/(z32[inf+1]-z32[inf])
      return z32[inf]+(Y-F32[inf])/m

  def F_32(fug):
    if (z32[-1]<=fug):
      return F32_A(fug)
    else:
      inf = 0
      sup = len(z32)-1
      while (inf<sup-1):
        med = (inf+sup)//2
        if(z32[med]<fug):
          inf = med
        else:
          sup = med
      m = (F32[inf+1]-F32[inf])/(z32[inf+1]-z32[inf])
      return F32[inf] + m*(fug-z32[inf])

  def F_52(fug):
    if (z52[-1]<=fug):
      return F52_A(fug)
    else:
      inf = 0
      sup = len(z52)-1
      while (inf<sup-1):
        med = (inf+sup)//2
        if(z52[med]<fug):
          inf = med
        else:
          sup = med
      m = (F52[inf+1]-F52[inf])/(z52[inf+1]-z52[inf])
      return F52[inf] + m*(fug-z52[inf])

  def P_rhoT(A):
    return F_52(dame_z(A))/A

if (tipo == "dist"):
  if (nargs >= 3):
    params = sys.argv[2]
  if (nargs >= 4):
    T = float(sys.argv[3])
  if (nargs >= 5):
    rho = float(sys.argv[4])
  if (nargs >= 6):
    Nsamp = int(sys.argv[5])

  f = open("params_%s/config_%.3f_%.3f.lammpstrj" %(params, rho, T), "r")
  Ekins = np.zeros(N*Nsamp)
  for i in range(Nsamp):
    for j in range(9):
      f.readline()
    for j in range(N):
      data = (f.readline()).split(" ")
      Ekins[i*N+j] = sum([float(data[-i])**2 for i in range(1,4)])/(2*m)
  ns, E_bins = np.histogram(Ekins, Nbins)
  ns = ns/(Nsamp*(E_bins[1]-E_bins[0]))
  Es = (E_bins[1:]+E_bins[:-1])/2

  params_aj, coso = sc.curve_fit(FD, Es/T, T*ns, long_term(T)**3*rho/4, bounds = ([0],[np.inf]))
  plt.figure()
  plt.plot(Es/T, T*ns, "k.")
  #plt.plot(Es/T, MB(Es/T), "r-")
  plt.plot(Es/T, FD(Es/T, long_term(T)**3*rho), "b-")
  plt.plot(Es/T, FD(Es/T, params_aj[0]), "k--")
  plt.xlabel(r"$\varepsilon_{cin}/T$", fontsize=20)
  plt.ylabel(r"$T f(\varepsilon_{cin})$", fontsize=20)
  #plt.legend(["Pauli", "Boltzmann", r"Fermi real $\lambda^3\rho=%.2f$" %(long_term(T)**3*rho), r"Fermi ajuste $\lambda^3\rho=%.2f$" %params_aj[0]])
  plt.legend(["Pauli", r"Fermi real $\lambda^3\rho=%.2f$" %(long_term(T)**3*rho), r"Fermi ajuste $\lambda^3\rho=%.2f$" %params_aj[0]])
  plt.title(r"params_%s | $T^*=%.1f$MeV | $\rho^*=%.2f$fm$^{-3}$" %(params, T, rho))


if (tipo == "termo"):
  if (nargs >= 3):
    params = sys.argv[2]
  if (nargs >= 4):
    rho = float(sys.argv[3])

  data = np.loadtxt("params_%s/termo_%.3f.txt" %(params, rho))
  Ts = data[:,1]
  Ekin = data[:,2]
  P = data[:,-3]
  Tv = data[:,-2]

  Ekin_FD = np.array([1.5*P_rhoT(long_term(T)**3*rho)*T for T in Ts])
  plt.figure()
  plt.plot(Ts, Ekin, "r.--")
  plt.plot(Ts, 1.5*P/rho, "b.--")
  plt.plot(Ts, Ekin_FD, "k-")
  plt.xlabel("$T$ [MeV]")
  plt.ylabel("Energia [MeV]")
  plt.title(r"params_%s | $\rho = %.3f$fm$^{-3}$" %(params, rho))
  plt.legend([r"$E_{cin}$", r"$1.5P/\rho$", "Fermi"], loc = "upper center")

pressure = ct.CDLL('../pressure.so')

if (tipo == "gr"):
  if (nargs >= 3):
    params = sys.argv[2]
  if (nargs >= 4):
    T = float(sys.argv[3])
  if (nargs >= 5):
    rho = float(sys.argv[4])
  if (nargs >= 6):
    Nsamp = int(sys.argv[5])

  gr_c = pressure.Gr
  gr_c.argtypes = [ct.c_voidp, ct.c_int, ct.c_float,
                       ct.c_float, ct.c_voidp, ct.c_int]
  gr_c.restype = ct.c_float

  L = (N/rho)**(1/3)
  #pairs = np.array(list(it.combinations(range(N), 2)), dtype=np.int64)
  def Gr(x, dr):
    M = int(np.ceil(np.sqrt(3)*L/(2*dr)))
    gr = np.zeros(M, dtype=np.float32)
    xp = x.ctypes.data_as(ct.c_voidp)
    grp = gr.ctypes.data_as(ct.c_voidp)
    gr_c(xp, len(x)//3, dr, L, grp, M)
    return gr

  def Gp(p, dp):
    max_dp = np.sqrt(3)*abs(max(p)-min(p))
    M = int(np.ceil(max_dp/dp))
    gp = np.zeros(M, dtype=np.float32)
    pp = p.ctypes.data_as(ct.c_voidp)
    gpp = gp.ctypes.data_as(ct.c_voidp)
    gr_c(pp, len(p)//3, dp, 5*max_dp, gpp, M)
    return gp

  f = open("params_%s/config_%.3f_%.3f.lammpstrj" %(params, rho, T), "r")
  qs = np.zeros(3*N, dtype = np.float32)
  ps = np.zeros(3*N*Nsamp, dtype = np.float32)
  dr = L/4000
  M = int(np.ceil(np.sqrt(3)*L/(2*dr)))
  gr_q = np.zeros(M, dtype = np.float32)
  for i in range(Nsamp):
    for j in range(9):
      f.readline()
    for j in range(N):
      data = (f.readline()).split(" ")
      qs[3*j:3*(j+1)] = np.array([float(data[i]) for i in range(2,5)], dtype = np.float32)
      ps[3*j+i*Nsamp:3*(j+1)+i*Nsamp] = np.array([float(data[i]) for i in range(2,5)], dtype = np.float32)
    gr_q += Gr(qs, dr)/Nsamp

  max_dp = np.sqrt(3)*abs(max(ps)-min(ps))
  dp = max_dp/1000
  M = int(np.ceil(max_dp/dp))
  gr_p = np.zeros(M, dtype = np.float32)
  for i in range(Nsamp):
    p_temp = np.array(ps[Nsamp*3*N*i:Nsamp*3*N*(i+1)], dtype=np.float32)
    gr_p += Gp(p_temp, dp)/Nsamp
  plt.figure()
  plt.plot(dr*(np.arange(len(gr_q))+0.5), gr_q, "-")
  plt.axis([0,L/2, 0, max(gr_q)*1.1])
  plt.figure()
  plt.plot(dp*(np.arange(M)+0.5), gr_p, "-")
  plt.show()
