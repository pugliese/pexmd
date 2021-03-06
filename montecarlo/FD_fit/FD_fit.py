import numpy as np
import matplotlib.pylab as plt
import scipy.optimize as sc
import sys
import glob

m = 1.043916 * 100

tipo = "v"
rho = 'rho0'
Nbins = 100

nargs = len(sys.argv)
if (nargs >= 2):
  tipo = sys.argv[1]
if (nargs >= 3):
  rho = sys.argv[2]
if (nargs >= 4):
  Nbins = int(sys.argv[3])

files = glob.glob("distribucion_"+rho+"_ideal_*")
Ts = [f.split("_")[3][:-4] for f in files]
n_temps = len(Ts)

for k in range(n_temps):
  if Ts[k][2]== "E":
    Ts[k] = 10**float(Ts[k][3:])
  else:
    Ts[k] = float(Ts[k])
idxs = np.argsort(Ts)
Ts = np.array([Ts[i] for i in idxs])
files = [files[i] for i in idxs]

ns = np.zeros((n_temps, Nbins))
bins = np.zeros((n_temps, Nbins+1))
p = np.zeros((n_temps, Nbins))
for k in range(n_temps):
  data_aux = np.loadtxt(files[k])
  data = np.zeros(len(data_aux)/3)
  for i in range(int(len(data_aux)/3)):
    data[i] = np.sum(data_aux[3*i:3*i+3]**2)

  ns[k, :], bins[k, :] = np.histogram(data, Nbins)

N_max = np.max(ns)
N_maxs = np.max(ns, axis=1)

if (tipo == "v"):
  plt.figure()
  for k in range(n_temps):
    if (n_temps > 1):
      plt.subplot(n_temps/3, 3, k+1)
    plt.xlabel(r"$p$")
    plt.hist(data, Nbins)
    plt.text((min(data)+ 5*max(data))/6, 0.9*max(ns[k,:]), "T=%f" %(Ts[k]))
  plt.show()

if (tipo == "f" or tipo == 'f&v'):
  MB = lambda x, A, mum, Tm: A*np.exp(-0.5*x**2/Tm)
  FD = lambda x, A, mum, Tm: A/(np.exp((0.5*x**2-mum)/Tm)+1)
  As = np.zeros(n_temps)
  mums = np.zeros(n_temps)
  Tms = np.zeros(n_temps)
  for k in range(n_temps):
    p[k,:] = (bins[k,1:] + bins[k,:-1])/2
    params, coso = sc.curve_fit(FD, p[k,:], ns[k,:], [10, 0, Ts[k]*m], bounds = ([0,-np.inf,0],[8*N_max,np.inf,np.inf]))
    As[k] = params[0]
    mums[k] = params[1]
    Tms[k] = params[2]

  if(tipo == 'f&v'):
    plt.figure()
    for k in range(0,n_temps-1,5):
      if (n_temps > 1):
        plt.subplot(n_temps/(5*2), 2, k/5+1)
      plt.xlabel(r"$p$")
      plt.plot(p[k,:],  ns[k,:], "ko")
      plt.plot(np.linspace(p[k,0], p[k,-1], 1000),  FD(np.linspace(p[k,0], p[k,-1], 1000), As[k], mums[k], Tms[k]), "b-")
      plt.text((min(data)+ 5*max(data))/6, 0.9*max(ns[k,:]), "T=%f" %(Ts[k]))
    plt.figure()
    plt.plot(Ts, Tms/m, "bo--")
    plt.plot(Ts, Ts, "r-")
    plt.xlabel(r"$T$")
    plt.ylabel(r"$T_{FD}$")
    plt.figure()
    #plt.semilogx(Ts, mums/m, "bo--")
    plt.plot(Ts, mums/m, "bo--")
    plt.xlabel(r"$T$")
    plt.ylabel(r"$\mu_{FD}$")
    print(Tms/(Ts*m))
    plt.show()
    np.savetxt('data_'+rho+'.txt', [Ts, Tms, mums, As])
