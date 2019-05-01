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
Nsamp = 100
N = 8000

nargs = len(sys.argv)
if (nargs >= 2):
  tipo = sys.argv[1]


MB = lambda E, T: N*np.sqrt(4*E/(np.pi*T**3))*np.exp(-E/T)
FD = lambda E, a: N*np.sqrt(4*E/(np.pi*T**3))/(np.exp(E/T)*(a/dame_z(a))+a)

if (True):
  T = 10
  rho = 0.2
  D = 100

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


if (tipo == "dist"):

  if (nargs >= 3):
    T = float(sys.argv[2])
  if (nargs >= 4):
    rho = float(sys.argv[3])
  if (nargs >= 5):
    D = float(sys.argv[4])

  data = np.loadtxt("prueba_pauli/dist_%.2f_%.2f_%.2f.txt" %(rho, D, T))
  ns, E_bins = np.histogram(data, Nbins)
  ns = ns/(Nsamp*(E_bins[1]-E_bins[0]))

  Es = (E_bins[1:]+E_bins[:-1])/2

  params, coso = sc.curve_fit(FD, Es, ns, long_term(T)**3*rho, bounds = ([0],[np.inf]))
  plt.figure()
  plt.plot(Es/T, ns, "k.")
  plt.plot(Es/T, MB(Es, T), "r-")
  plt.plot(Es/T, FD(Es, params[0]), "b-")
  plt.xlabel(r"$\varepsilon_{cin}$", fontsize=20)
  plt.ylabel(r"$f(\varepsilon_{cin})$", fontsize=20)
  plt.legend(["Pauli", "Boltzmann", r"Fermi $\lambda^3\rho=%.2f$" %params[0]])
  plt.title(r"$D^* = %.1f$ | $T^*=%.1f$ | $\rho^*=%.2f$" %(D, T, rho))
  #plt.plot(Es, FD(Es, np.log(dame_z(long_term(Ts[k])**3*rho))*T, T), "b-")

if (tipo == "aj"):

  rho = 0.2
  D = 100
  if (nargs >= 3):
    rho = float(sys.argv[2])
  if (nargs >= 4):
    D = float(sys.argv[3])

  Ts = np.linspace(2.5, 30, 12)
  Ts_lindo = np.linspace(2.5,30,120)
  leyenda = []
  Ds = np.array([10.0, 50.0, 75.0, 100.0, 300.0, 500.0, 1000.0])
  Ds = np.array([10.0, 50.0, 100.0, 500.0, 1000.0])
  rhos = np.array([0.2, 0.4, 0.6, 0.8, 1.0])
  plt.figure()
  Asss = []
  for rho in rhos:
    Ass = []
    for D in Ds:
      As = []
      for T in Ts:
        data = np.loadtxt("prueba_pauli/dist_%.2f_%.2f_%.2f.txt" %(rho, D, T))
        ns, E_bins = np.histogram(data, Nbins)
        ns = ns/(Nsamp*(E_bins[1]-E_bins[0]))

        Es = (E_bins[1:]+E_bins[:-1])/2

        params, coso = sc.curve_fit(FD, Es, ns, long_term(T)**3*rho, bounds = ([0],[np.inf]))
        As.append(params[0])
      Ass.append(As)
      plt.loglog(Ts, As, "o--")
      leyenda.append(r"$D^*=%.2f$, $\rho^*=%.2f$" %(D,rho))
    Asss.append(Ass)
  #plt.plot(Ts_lindo, As[3]*(Ts[3]/Ts_lindo)**1.5, "r-")
  plt.xlabel("T")
  plt.ylabel(r"$\lambda^3\rho$")
  plt.legend(leyenda)
  Asss = np.array(Asss)

if (tipo == "ajvsrho"):
  Ts = np.linspace(2.5, 30, 12)
  Ts_lindo = np.linspace(2.5,30,120)
  Ds = np.array([10.0, 50.0, 75.0, 100.0, 300.0, 500.0, 1000.0])
  Ds = np.array([10.0, 50.0, 100.0, 500.0, 1000.0])
  rhos = np.array([0.2, 0.4, 0.6, 0.8, 1.0])
  markers = ["o--", "v--", "^--", "s--", "d--"]
  colors = ["k","r","b","g","y"]
  plt.figure()
  for j in range(len(rhos)):
    leyenda = []
    data = np.loadtxt("prueba_pauli/As_%.1f_vs_D_T.txt" %rhos[j])
    for i in range(len(Ds)):
      plt.loglog(Ds[i]**(1/3.5)/Ts, data[i,:]/(rhos[j]**0.8), colors[j]+markers[i])
      leyenda.append(r"$D^*=%.1f$" %Ds[i])
    plt.xlabel(r"$(D^*)^\frac{2}{7}/T^*$", fontsize=16)
    plt.ylabel(r"$A_{aj}/(\rho^*)^\frac{4}{5}$", fontsize=16)
    plt.title(r"$\rho^*$ = 0.2 | 0.4 | 0.6 | 0.8 | 1.0", fontsize=16)
    plt.legend(leyenda, loc = "lower right")
  plt.show()

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
  for d in range(-1,2):
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
  plt.title(r"$\rho^*=%.2f$" %rho)
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
  for D in [1, 10, 50, 100, 500, 1000, 5000]:
    data = np.loadtxt("prueba_pauli/termo_pauli_%.2f_%.2f.txt" %(rho, D))
    plt.plot(data[:,1], data[:,2]/data[:,1], '.--')
    leyenda.append("$D^*=%.1f$" %(D))
  plt.plot(data[:,1], 1.5*data[:,1]/data[:,1], "r-")
  plt.plot(data[:,1], 2*data[:,1]/data[:,1], 'k--')
  leyenda.append("Boltzmann")
  plt.legend(leyenda, loc = 'upper center')
  plt.xlabel("T", fontsize=14)
  plt.ylabel(r"$E_{cin}$", fontsize=14)

if (tipo == "P"):
  rho = 0.2
  D = 100
  if (nargs >= 3):
    rho = float(sys.argv[2])
  if (nargs >= 4):
    D = float(sys.argv[3])
  leyenda = []
  plt.figure()
  for D in [10, 50, 75, 100, 300, 500, 1000]:
    data = np.loadtxt("prueba_pauli/termo_pauli_%.2f_%.2f.txt" %(rho, D))
    plt.plot(data[:,1], data[:,-3], '.--')
    leyenda.append("$D^*=%.1f$" %(D))
  plt.plot(data[:,1], rho*data[:,1], "r-")
  leyenda.append("Boltzmann")
  plt.legend(leyenda, loc = 'upper center')
  plt.xlabel("$T$", fontsize=14)
  plt.ylabel("$P$", fontsize=14)
  plt.title(r"$\rho^*=%.2f$" %rho)


if (tipo == "TvsT"):
  rho = 0.2
  Do = 100
  if (nargs >= 3):
    rho = float(sys.argv[2])
  if (nargs >= 4):
    Do = float(sys.argv[3])
  leyenda = []
  data_Do = np.loadtxt("prueba_pauli/termo_pauli_%.2f_%.2f.txt" %(rho, Do))
  Ts_Do = data_Do[90:-8,1]
  Es = data_Do[90:-8,2]/Ts_Do
  Ds = np.array([10, 30, 50, 75, 100, 300, 500, 750, 1000, 5000], dtype = np.float32)
  if (rho!=0.2):
    Ds = np.array([10, 50, 100, 500, 1000], dtype = np.float32)
  Ts = np.zeros((len(Ds), len(Ts_Do)))
  ms = np.zeros_like(Ds)
  bs = np.zeros_like(Ds)
  alfas = np.zeros_like(Ds)
  j = 0
  leyenda = []
  plt.figure()
  for D in Ds:
    data_D = np.loadtxt("prueba_pauli/termo_pauli_%.2f_%.2f.txt" %(rho, D))
    i = 0
    for Ekin in Es:
      deltas_Ekin = np.abs(data_D[:,2]/data_D[:,1]-Ekin)
      idx = np.where(deltas_Ekin==np.min(deltas_Ekin))[0][0]
      Ts[j, i] = data_D[idx,1]
      i += 1
    coefs = np.polyfit(np.log(Ts_Do), np.log(Ts[j,:]), 1)
    ms[j], bs[j] = coefs[0], coefs[1]
    #alfas[j] = np.log(D/100.0)/np.log(ms[j])
    plt.loglog(Ts_Do, Ts[j,:], ".--")
    leyenda.append("$D^*=%.2f$" %D)
    j += 1
  plt.legend(leyenda, loc = 'upper left')
  plt.xlabel("T", fontsize=14)
  plt.ylabel(r"$T_{equiv}$", fontsize=14)
  plt.figure()
  #plt.plot(np.log(Ds/Do), bs, "bo-")
  #plt.xlabel("log($D/D_0$)")
  plt.plot(np.log(2*Ds), np.exp(bs), "bo-")
  plt.xlabel("log(2D)")
  plt.ylabel("$f(D; D_0)$")
  plt.title("$D_0 = %.2f$" %Do)
  plt.grid()
  plt.figure()
  for i in range(len(Es)):
    plt.plot(np.log(2*Ds), Ts[:, i], ".--")
  plt.xlabel("log(2D)", fontsize=14)
  plt.ylabel("T", fontsize=14)
  plt.show()


  """
if (tipo == "TvsT"):
  rho = 0.2
  D = 100
  if (nargs >= 3):
    rho = float(sys.argv[2])
  if (nargs >= 4):
    D = float(sys.argv[3])
  leyenda = []
  data_100 = np.loadtxt("prueba_pauli/termo_pauli_%.2f_%.2f.txt" %(rho, 100))
  m_1000 = np.zeros(40)
  b_1000 = np.zeros(40)
  m_10 = np.zeros(40)
  b_10 = np.zeros(40)
  alfa_10 = np.zeros(40)
  alfa_1000 = np.zeros(40)
  for i in range(40):
    Ts_100 = data_100[60:-(i+1),1]
    Es = data_100[60:-(i+1),2]/Ts_100
    Ts_10 = []
    data_10 = np.loadtxt("prueba_pauli/termo_pauli_%.2f_%.2f.txt" %(rho, 10))
    Ts_1000 = []
    data_1000 = np.loadtxt("prueba_pauli/termo_pauli_%.2f_%.2f.txt" %(rho, 1000))
    for Ekin in Es:
      deltas_Ekin = np.abs(data_10[:,2]/data_10[:,1]-Ekin)
      idx = np.where(deltas_Ekin==np.min(deltas_Ekin))[0][0]
      Ts_10.append(data_10[idx,1])
      j_10 = idx
      deltas_Ekin = np.abs(data_1000[:,2]/data_1000[:,1]-Ekin)
      idx = np.where(deltas_Ekin==np.min(deltas_Ekin))[0][0]
      Ts_1000.append(data_1000[idx,1])
      j_1000 = idx
    Ts_10 = np.array(Ts_10)
    Ts_1000 = np.array(Ts_1000)
    m_10[i], b_10[i] = np.polyfit(Ts_100, Ts_10, 1)
    m_1000[i], b_1000[i] = np.polyfit(Ts_100, Ts_1000, 1)
    alfa_10[i] = np.log(0.1)/np.log(m_10[i])
    alfa_1000[i] = np.log(10)/np.log(m_1000[i])
  plt.figure()
  plt.plot(alfa_10, "r.")
  plt.plot(alfa_1000, "b.")
  """
  """
  plt.figure()
  plt.plot(Ts_100, Ts_10, "r.")
  plt.plot(Ts_100, b_10+Ts_100*m_10, "r-")
  plt.plot(Ts_100, Ts_1000, "b.")
  plt.plot(Ts_100, b_1000+m_1000*Ts_100, "b-")
  print(r"D* = 10: ", m_10, b_10, "--->", np.log(0.1)/np.log(m_10))
  print(r"D*=1000: ", m_1000, b_1000, "--->", np.log(10)/np.log(m_1000))
  plt.legend(["$D^*=10$", "Fit $D^*=10$","$D^*=1000$", "Fit $D^*=1000$"], loc = 'upper center')
  plt.xlabel("T", fontsize=14)
  plt.ylabel(r"$T_{equiv}$", fontsize=14)
  """
