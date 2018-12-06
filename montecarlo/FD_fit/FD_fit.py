import numpy as np
import matplotlib.pylab as plt
#import scipy.optimize as sc
import sys


n_temps = 10
log10_Tsup =  0
log10_Tinf =  -4
log10_step = (log10_Tsup - log10_Tinf)/(n_temps - 1)

tipo = "v"
Nbins = 50

rho = 'rho0'

Ts = 10**np.linspace(log10_Tinf, log10_Tsup, n_temps)


ns = np.zeros((n_temps, Nbins))
bins = np.zeros((n_temps, Nbins+1))
p = np.zeros((n_temps, Nbins))
for k in range(n_temps):
  filename = "distribucion_"+rho+"_%f.txt"%(Ts[k])
  data = np.loadtxt(filename)
  ns[k, :], bins[k, :] = np.histogram(abs(data), Nbins)

N_max = np.max(ns)

if (tipo == "v"):
  plt.figure()
  for k in range(n_temps):
    filename = "distribucion_"+rho+"_%f.txt"%(Ts[k])
    data = np.loadtxt(filename)
    if (n_temps > 1):
        plt.subplot(n_temps/2, 2, k+1)
    plt.xlabel(r"$p$")
    plt.hist(abs(data), Nbins)
    plt.text((min(abs(data))+ 5*max(abs(data)))/6, 0.9*max(ns[k,:]), "T=%f" %(Ts[k]))
    #plt.axis([min(abs(data)), max(abs(data)), 0, N_max])
  plt.show()

if (tipo == "f" or tipo == 'f&v'):
  FD = lambda x, A, muT, Tm: A/(np.exp((0.5*x**2-muT)/Tm)+1)
  As = np.zeros(n_temps)
  muTs = np.zeros(n_temps)
  Tms = np.zeros(n_temps)
  for k in range(n_temps):
    p[k,:] = (bins[k,1:] + bins[k,:-1])/2
    params, coso = sc.curve_fit(FD, p[k,:], ns[k,:], [N_max, 0, Ts[k]], bounds = ([0,0,0],[np.inf,np.inf,np.inf]))
    As[k] = params[0]
    muTs[k] = params[1]
    Tms[k] = params[2]

  if(tipo == 'f&v'):
    plt.figure()
    for k in range(n_temps):
      Ts[k] = 10**(log10_Tsup - log10_step*k)
      filename = "distribucion_"+rho+"_%f.txt"%(Ts[k])
      data = np.loadtxt(filename)
      if (n_temps > 1):
        plt.subplot(n_temps/2, 2, k+1)
      plt.xlabel(r"$p$")
      plt.plot(p[k,:],  ns[k,:], "ko")
      plt.plot(np.linspace(p[k,0], p[k,-1], 1000),  FD(np.linspace(p[k,0], p[k,-1], 1000), As[k], muTs[k], Tms[k]), "b-")
      plt.plot(np.linspace(p[k,0], p[k,-1], 1000),  FB(np.linspace(p[k,0], p[k,-1], 1000), As_B[k], Tms_B[k]), "r--")
      plt.text((min(abs(data))+ 5*max(abs(data)))/6, 0.9*max(ns[k,:]), "T=%f" %(Ts[k]))
      #plt.axis([min(abs(data)), max(abs(data)), 0, N_max])
      plt.show()
    np.savetxt('data.txt', Ts, Tms, muTs, As)
