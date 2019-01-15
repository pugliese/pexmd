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
N = 18**3
Nsamp = 200

tipo = "rep"
rho = 'rho0'
Nbins = 100
Nreps = 1

nargs = len(sys.argv)
if (nargs >= 2):
  tipo = sys.argv[1]
if (nargs >= 3):
  rho = sys.argv[2]
if (nargs >= 4):
  Nbins = int(sys.argv[3])

L = 12*N**(1.0/3)
if (rho[-1] == "0"):
  L = L/3
if (rho[-1] == "1"):
  L = L/2
V = L**3

files = glob.glob("energias_*")
rhos = [float(f.split("_")[1][:-4]) for f in files]
n_rhos = len(rhos)
rhos = np.sort(rhos)
E = []
E_kin = []
Es = np.zeros(n_rhos)

for k in range(n_rhos):
  data = np.loadtxt("energias_%f.txt" %(rhos[k]))
  E.append(data[:,0])
  E_kin.append(data[:,1])
  plt.subplot(n_rhos//2, 2, k+1)
  plt.plot(E[-1]/N, "b-")
  plt.plot(E_kin[-1]/N, "r--")
  plt.legend([r"$\rho=%f fm^{-3}$" %rhos[k]])
  Es[k] = np.mean(E[-1])/N
plt.figure()
plt.plot(rhos, Es, "o--")
plt.show()
