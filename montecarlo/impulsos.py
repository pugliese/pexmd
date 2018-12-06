import numpy as np
import matplotlib.pylab as plt
import scipy.optimize as sc
import sys

Tinit =  0.003
T_step = 0.0005
formato = "%f"
tipo = "n"
Nbins = 50
masa = "mx100_Lx1.25_"

if (len(sys.argv) >= 2):
  orden = int(sys.argv[1])
  Tinit =  3*10**orden
  T_step = 0.5*10**orden

if (len(sys.argv)==3):
  tipo = sys.argv[2]

ns = np.zeros((6, Nbins))
bins = np.zeros((6, Nbins+1))
p = np.zeros((6, Nbins))
for k in range(6):
  To = Tinit - T_step*k
  filename = "impulsos/impulsos_"+masa+formato%(To)+".txt"
  data = np.loadtxt(filename)
  ns[k, :], bins[k, :] = np.histogram(abs(data), Nbins)

N_max = np.max(ns)

if (tipo == "n" or tipo == "a"):
  plt.figure()
  for k in range(6):
    To = Tinit - T_step*k
    filename = "impulsos/impulsos_"+masa+formato%(To)+".txt"
    data = np.loadtxt(filename)
    plt.subplot(3, 2, k+1)
    plt.xlabel(r"$p$")
    plt.hist(abs(data), 50)
    plt.text((min(abs(data))+ 5*max(abs(data)))/6, 0.9*N_max, "T="+formato %(To))
    plt.axis([min(abs(data)), max(abs(data)), 0, N_max])
  plt.show()
if (tipo == "r" or tipo == "a"):
  plt.figure()
  for k in range(6):
    To = Tinit - T_step*k
    filename = "impulsos/impulsos_"+masa+formato%(To)+".txt"
    data = np.loadtxt(filename)
    data = data/np.sqrt(To)
    plt.subplot(3, 2, k+1)
    plt.xlabel(r"$p/\sqrt{T}$")
    plt.hist(abs(data), 50)
    plt.text((min(abs(data))+ 5*max(abs(data)))/6, 0.9*N_max, "T="+formato %(To))
    plt.axis([min(abs(data)), max(abs(data)), 0, N_max])
  plt.show()

if (tipo == "f"):
  FD = lambda x, A, muT, Tm: A/(np.exp((0.5*x**2-muT)/Tm)+1)
  FB = lambda x, A, Tm: A*np.exp(-0.5*x**2/Tm)
  As = np.zeros(6)
  muTs = np.zeros(6)
  Tms_B = np.zeros(6)
  As_B = np.zeros(6)
  Tms = np.zeros(6)
  for k in range(6):
    p[k,:] = (bins[k,1:] + bins[k,:-1])/2
    params, coso = sc.curve_fit(FD, p[k,:], ns[k,:], [N_max, 0, 100000], bounds = ([0,0,0],[np.inf,np.inf,np.inf]))
    As[k] = params[0]
    muTs[k] = params[1]
    Tms[k] = params[2]
    params, coso = sc.curve_fit(FB, p[k,:], ns[k,:], [N_max, 1000], bounds = ([0,0],[np.inf,np.inf]))
    As_B[k] = params[0]
    Tms_B[k] = params[1]

  plt.figure()
  for k in range(6):
    To = Tinit - T_step*k
    filename = "impulsos/impulsos_"+masa+formato%(To)+".txt"
    data = np.loadtxt(filename)
    plt.subplot(3, 2, k+1)
    plt.xlabel(r"$p$")
    plt.plot(p[k,:],  ns[k,:], "ko")
    plt.plot(np.linspace(p[k,0], p[k,-1], 1000),  FD(np.linspace(p[k,0], p[k,-1], 1000), As[k], muTs[k], Tms[k]), "b-")
    plt.plot(np.linspace(p[k,0], p[k,-1], 1000),  FB(np.linspace(p[k,0], p[k,-1], 1000), As_B[k], Tms_B[k]), "r--")
    plt.text((min(abs(data))+ 5*max(abs(data)))/6, 0.9*N_max, "T="+formato %(To))
    plt.axis([min(abs(data)), max(abs(data)), 0, N_max])
  print(As, muTs, Tms)
  plt.show()
