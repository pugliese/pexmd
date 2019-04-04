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


tipo = "rep"
rho = 'rho0'
Nbins = 100
Nreps = 8
Nsamp = 200

nargs = len(sys.argv)
if (nargs >= 2):
  tipo = sys.argv[1]
if (nargs >= 3):
  rho = sys.argv[2]
if (nargs >= 4):
  Nbins = int(sys.argv[3])


N = 10**3
if(rho[0]== "M"):
  L = 2*1.644*N**(1.0/3)
  m = 938.27203 * 100
  h_bar = 197.327
  qo = 1.644
  po = 120
  D = 207*(h_bar/(po*qo))**3
else:
  L = 12*N**(1.0/3)
  m = 1.043916 * 100
  h_bar =  6.582119
  qo = 6
  po = 2.067 # En MeV 1E-22 s/fm
  D = 34.32*(h_bar/(po*qo))**3
  factor_p = 29.98 # Para convertir a MeV/c

if (rho[-1] == "0"):
  L = L/3
if (rho[-1] == "1"):
  L = L/2
V = L**3

scut = np.sqrt(10)
pauli = [D, qo, po, np.sqrt(10)]
h = h_bar*2*np.pi

files = glob.glob(rho+"/distribucion_10_rep1_*")
#files = glob.glob(rho+"/energia_*")
Ts = [f.split("_")[3][:-4] for f in files]
#Ts = [f.split("_")[1][:-4] for f in files]
n_temps = len(Ts)

for k in range(n_temps):
  if Ts[k][2]== "E":
    Ts[k] = 10**float(Ts[k][3:])
  else:
    Ts[k] = float(Ts[k])
idxs = np.argsort(Ts)
Ts = np.array([Ts[i] for i in idxs])
files = [files[i] for i in idxs]

pressure = ct.CDLL('../pressure.so')


if(tipo == 't'):
  E_kin = []
  E = []
  acept = []

  for k in range(n_temps):
    data = np.loadtxt(rho+"/energia_%f.txt" %(Ts[k]), dtype=np.float32)
    E_kin.append(data[:,0])
    E.append(data[:,1])
    acept.append(data[:,2])
    plt.subplot(n_temps//2, 2, k+1)
    plt.plot(E[-1]/N, "b-")
    plt.legend([r"$T=%f fm^{-3}$" %Ts[k]])
  plt.figure()
  for k in range(n_temps):
    plt.subplot(n_temps//2, 2, k+1)
    plt.plot(acept[k]/(20*np.arange(1,len(acept[k])+1)), "b-")
    plt.legend([r"$T=%f fm^{-3}$" %Ts[k]])
  plt.show()

if (tipo == "h" or tipo == "h&p"):

  pressure_c = pressure.pressure_pauli_PBC
  pressure_c.argtypes = [ct.c_voidp, ct.c_voidp, ct.c_voidp, ct.c_longlong, ct.c_float,
                       ct.c_float, ct.c_float, ct.c_float, ct.c_voidp, ct.c_voidp, ct.c_float]
  pressure_c.restype = ct.c_float

  pairs = np.array(list(it.combinations(range(N), 2)), dtype=np.int64)
  def pressure(x, p):
    energ = 0
    force = np.zeros_like(x, dtype=np.float32)
    gorce = np.zeros_like(x, dtype=np.float32)
    xp = x.ctypes.data_as(ct.c_voidp)
    pp = p.ctypes.data_as(ct.c_voidp)
    pairsp = pairs.ctypes.data_as(ct.c_voidp)
    forcep = force.ctypes.data_as(ct.c_voidp)
    gorcep = gorce.ctypes.data_as(ct.c_voidp)
    qF = pressure_c(xp, pp, pairsp, len(pairs), D, qo, po, scut, forcep, gorcep, L)
    return force, gorce, qF

  P = np.zeros(n_temps)
  Tv = np.zeros(n_temps)
  qp = np.zeros(n_temps)
  pp = np.zeros(n_temps)
  qF = np.zeros(n_temps)

  for k in range(n_temps):
    filename = rho+"/histograma_%d_T=%f.txt" %(Nbins, Ts[k])
    print(filename)
    data_q = np.array([], dtype=np.float32)
    data_p = np.array([], dtype=np.float32)
    for j in range(Nreps):
      data_aux = np.loadtxt(rho+"/distribucion_10_rep%d_%f.txt" %(j+1, Ts[k]), dtype=np.float32)
      data_q = np.concatenate([data_q, data_aux[:, 0]])
      data_p = np.concatenate([data_p, data_aux[:, 1]])
    print("Presion")
    if(tipo == "h&p"):
      for i in range(Nsamp*Nreps):
        fuerzas, guerzas, qF_aux = pressure(data_q[3*N*i:3*N*(i+1)], data_p[3*N*i:3*N*(i+1)])
        qp[k] += np.sum(data_p[3*N*i:3*N*(i+1)]*(data_p[3*N*i:3*N*(i+1)]/m-guerzas))/(3*V*Nsamp*Nreps)
        qF[k] += qF_aux/(3*V*Nsamp*Nreps)
      P[k] = qF[k] + qp[k]
      Tv[k] = qp[k]*V/N

    data = np.zeros(len(data_p)//3)
    for i in range(len(data_p)//3):
      data[i] = np.sum(data_p[3*i:3*i+3]**2)/(2*m)
    ns, bins = np.histogram(data, Nbins)
    ns_q, bins_q = np.histogram(data_q, Nbins)
    ns = ns/((bins[1]-bins[0])*Nreps*Nsamp)
    ns_q = ns_q/((bins_q[1]-bins_q[0])*3*Nreps*Nsamp)
    E = (bins[1:] + bins[:-1])/2
    np.savetxt(filename, [E, ns])

  if(tipo == "h&p"):
    np.savetxt(rho+"/presiones_N=1000.txt", [Ts,P])
    np.savetxt(rho+"/temperaturas_N=1000.txt", [Ts,Tv])


deg = lambda E: 4*np.pi*np.sqrt(2*m**3*E)*(L/h)**3
long_term = lambda T: (2*np.pi*h_bar**2/(m*T))**0.5
Ef = h_bar**2/(2*m)*(6*np.pi**2*N/V)**(2/3)

data = np.loadtxt("LUT_F32.txt")
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


if (tipo == "f&vg"):
  MB = lambda x, T: N*2*np.sqrt(x/np.pi)*np.exp(-x/T)/(T**1.5)
  FD = lambda x, mu, T: deg(x)/(np.exp((x-mu)/T)+1)

  if (rho[0]!="M"):
    Ts_quiero = Ts[8:]
    seleccion = [list(Ts).index(T) for T in Ts_quiero]
  if (rho[0] == "M"):
    seleccion = range(len(Ts))

  mus = np.zeros(len(seleccion))
  mus_aj = np.zeros(len(seleccion))
  Ts_aj = np.zeros(len(seleccion))
  i = 0
  for k in seleccion:
    plt.figure()
    data = np.loadtxt(rho+"/histograma_%d_T=%f.txt" %(Nbins, Ts[k]))
    E = data[0,:]
    ns = data[1,:]
    plt.plot(E,  ns/1000, "ko")
    rango = np.linspace(0, E[-1], 10000)
    exacto_MB = MB(rango, Ts[k])/1000
    plt.plot(rango, exacto_MB, "r-")
    print(dame_z(long_term(Ts[k])**3*N/V), long_term(Ts[k])**3*N/V)
    mus[i] = np.log(dame_z(long_term(Ts[k])**3*N/V))*Ts[k]
    exacto_FD = FD(rango, mus[i], Ts[k])/1000

    params, coso = sc.curve_fit(FD, E, ns, [0, Ts[k]], bounds = ([-np.inf,0],[Ef,np.inf]))

    plt.plot(rango, exacto_FD, "k--")
    plt.plot(rango, FD(rango, params[0], params[1])/1000, "b-")
    plt.xlabel(r"$E$", fontsize=20)
    plt.ylabel(r"$f(E)$", fontsize=20)
    plt.xticks([])
    plt.yticks([])
    plt.title(r"Parametros Dorso - $\rho^* = 3.375$", fontsize=20)
    if (rho[0] == "M"):
      plt.title(r"Parametros Maruyama - $\rho^* = 3.375$", fontsize=20)
    plt.axis([0, max(E), 0, 1.05*max(exacto_MB)])
    plt.legend(["Data", "Boltzmann", "Fermi", "Ajuste FD"], fontsize=20)
    plt.text((min(E)+ 4*max(E))/6, 0.4*max(exacto_MB), "T=%1.2fMeV" %(Ts[k]), fontsize=16)
    plt.savefig(rho+"/gif/histo_%d.png" %(len(seleccion)-i-1))
    plt.close()
    i+=1

if (tipo == "f&v"):
  MB = lambda x, T: N*2*np.sqrt(x/np.pi)*np.exp(-x/T)/(T**1.5)
  FD = lambda x, mu, T: deg(x)/(np.exp((x-mu)/T)+1)

  if (rho[0]!="M"):
    Ts_quiero = [0.01, 0.03, 0.05, 0.1, 0.3, 0.5, 1, 2, 3]
    seleccion = [list(Ts).index(T) for T in Ts_quiero]
  if (rho[0] == "M"):
    seleccion = range(9)
    if (rho[-1] == "0"):
      seleccion = [0,1,3,5,7,9,11,12,13]

  mus = np.zeros(len(seleccion))
  mus_aj = np.zeros(len(seleccion))
  Ts_aj = np.zeros(len(seleccion))
  i = 0
  plt.figure()
  plt.title(r"$\rho^* = %f$" %((qo/L)**3))
  for k in seleccion:
    data = np.loadtxt(rho+"/histograma_%d_T=%f.txt" %(Nbins, Ts[k]))
    E = data[0,:]
    ns = data[1,:]
    if (n_temps > 1):
      plt.subplot(3, 3, i+1)
      plt.subplots_adjust(left = 0.05, right = 0.98, bottom = 0.05, top = 0.98, wspace = 0.1, hspace = 0.125)
    plt.plot(E,  ns/1000, "ko")
    rango = np.linspace(0, E[-1], 10000)
    exacto_MB = MB(rango, Ts[k])/1000
    plt.plot(rango, exacto_MB, "r-")
    print(dame_z(long_term(Ts[k])**3*N/V), long_term(Ts[k])**3*N/V)
    mus[i] = np.log(dame_z(long_term(Ts[k])**3*N/V))*Ts[k]
    exacto_FD = FD(rango, mus[i], Ts[k])/1000

    params, coso = sc.curve_fit(FD, E, ns, [0, Ts[k]], bounds = ([-np.inf,0],[Ef,np.inf]))
    mus_aj[i] = params[0]
    Ts_aj[i] = params[1]

    plt.plot(rango, exacto_FD, "k--")
    plt.plot(rango, FD(rango, params[0], params[1])/1000, "b-")
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.axis([0, max(E), 0, 1.05*max(exacto_MB)])
    plt.text((min(E)+ 4*max(E))/6, 0.4*max(exacto_MB), "T=%1.3fMeV" %(Ts[k]), fontsize=16)
    if (i==7):
      plt.xlabel(r"$\varepsilon$ [$MeV$]", fontsize=18)
    if (i==3):
      plt.ylabel(r"$f(E)$ [$GeV^{-1}$]", fontsize=18)
    if (i==4):
      plt.legend(["Data", "Boltzmann", "Fermi", "Ajuste FD"], fontsize=20)
    if (i==1):
      plt.text(0.3*max(E), 0.865*max(exacto_MB), r"$\rho^*$ = %1.3f"%(qo**3*N/V), fontsize=26, bbox=dict(facecolor='grey', edgecolor='black', boxstyle='round,pad=1', alpha=0.4))
    if (i==0):
      if (rho[0]=="M"):
        plt.text(0.3*max(E), 0.75*max(exacto_MB), "Parametros\nMaruyama", fontsize=26, bbox=dict(facecolor='grey', edgecolor='black', boxstyle='round,pad=1', alpha=0.4))
      else:
        plt.text(0.3*max(E), 0.75*max(exacto_MB), "Parametros\nDorso", fontsize=26, bbox=dict(facecolor='grey', edgecolor='black', boxstyle='round,pad=1', alpha=0.4))
    i+=1
  plt.show()
  mus = np.array([np.log(dame_z(long_term(Ts[k])**3*N/V))*T for T in Ts])


if (tipo == "fp"):
  FD = lambda x, mu, T: deg(x)/(np.exp((x-mu)/T)+1)

  pre_rho = rho[:-4]
  if (rho[0]!="M"):
    Ts_quiero = [0.01, 0.03, 0.05, 0.1, 0.3, 0.5, 1, 2, 3]
  if (rho[0]=="M"):
    Ts_quiero = [0.5, 1, 2, 3, 4, 5, 10, 15, 20]
  seleccion = [list(Ts).index(T) for T in Ts_quiero]

  files = glob.glob(pre_rho+"rho2"+"/distribucion_10_rep1_*")
  Ts = [float(f.split("_")[3][:-4]) for f in files]
  Ts = np.sort(Ts)
  n_temps = len(Ts)

  mus = np.zeros((3, n_temps))
  mus_aj = np.zeros((3, n_temps))
  Ts_aj = np.zeros((3, n_temps))
  densidad = np.zeros((3, n_temps))

  L = 10*qo*2*np.array([1.0/3,1.0/2,1.0])
  V = L**3

  Ef = h_bar**2/(2*m)*(6*np.pi**2*N/V)**(2/3)
  for j in range(3):
    rho = pre_rho + "rho" + str(j)
    deg = lambda E: 4*np.pi*np.sqrt(2*m**3*E)*(L[j]/h)**3
    for k in range(n_temps):
      data = np.loadtxt(rho+"/histograma_%d_T=%f.txt" %(Nbins, Ts[k]))
      E = data[0,:]
      ns = data[1,:]
      mus[j,k] = np.log(dame_z(long_term(Ts[k])**3*N/V[j]))*Ts[k]
      params, coso = sc.curve_fit(FD, E, ns, [0, Ts[k]], bounds = ([-np.inf,0],[Ef[j],np.inf]))
      mus_aj[j,k] = params[0]
      Ts_aj[j,k] = params[1]
      densidad[j,k] = F_32(np.exp(mus_aj[j,k]/Ts_aj[j,k]))/long_term(Ts_aj[j,k])**3

  for j in range(3):
    plt.figure()
    plt.plot(Ts, mus_aj[j,:], "o-", linewidth=2)
    plt.plot(Ts, mus[j,:],"s-", linewidth=2)
    plt.legend(["Ajustada", "Fermi"], loc = 3)
    if (rho[0] == "M"):
      Tmax = 20
    else:
      Tmax = 3
    mu_min = min([min(mus[j,:]), min(mus_aj[j,:])])
    plt.plot([0,0.46*Tmax], [Ef[j], Ef[j]], "k--")
    plt.plot([0.54*Tmax,1.05*Tmax], [Ef[j], Ef[j]], "k--")
    plt.text(Tmax*0.475, Ef[j]+mu_min/50, r"$E_F$", fontsize=16)
    plt.axis([0,Tmax*1.05,1.1*mu_min, -mu_min/10])
    plt.xlabel("$T$ [$MeV$]", fontsize=14)
    plt.ylabel(r"$\mu$ [$MeV$]", fontsize=14)
    if (rho[0]=="M"):
      plt.title(r"Parametros Maruyama - $\rho^*=%1.3f$" %(qo**3*N/V[j]))
    if (rho[0]!="M"):
      plt.title(r"Parametros Dorso - $\rho^*=%1.3f$" %(qo**3*N/V[j]))

  plt.figure()
  plt.plot(Ts, densidad[0,:]*V[0]/N, "o--", linewidth=2)
  plt.plot(Ts, densidad[1,:]*V[1]/N, "s--", linewidth=2)
  plt.plot(Ts, densidad[2,:]*V[2]/N, "^--", linewidth=2)
  plt.plot([-0.1,Ts[-1]*1.05], [1,1], "k-")
  plt.xlabel(r"$T$ [$MeV$]", fontsize=14)
  plt.ylabel(r"$\rho_{aj}/\rho$", fontsize=14)
  plt.legend([r"$\rho^*$ = %1.3f" %(qo**3*N/v) for v in V])
  if (rho[0]=="M"):
    plt.title("Parametros Maruyama")
    plt.axis([-0.1,Ts[-1]*1.05, 0.95, 1.2])
  if (rho[0]!="M"):
    plt.title("Parametros Dorso")
    plt.axis([-0.1,Ts[-1]*1.05, 0.9, 1.8])

  plt.figure()
  plt.plot(Ts, Ts_aj[0,:], "o-", linewidth=2)
  plt.plot(Ts, Ts_aj[1,:], "s-", linewidth=2)
  plt.plot(Ts, Ts_aj[2,:], "^-", linewidth=2)
  plt.plot([0,max(Ts)*1.05],[0,max(Ts)*1.05],"r-")
  plt.xlabel("$T$ [$MeV$]", fontsize=14)
  plt.ylabel("$T_{aj}$ [$MeV$]", fontsize=14)
  plt.axis([0,max(Ts)*1.05, 0, max(Ts)*1.05])
  plt.legend([r"$\rho^*$ = %1.3f" %(qo**3*N/v) for v in V], loc=2)
  if (rho[0]=="M"):
    plt.title("Parametros Maruyama")
  if (rho[0]!="M"):
    plt.title("Parametros Dorso")
  plt.show()

if (tipo == "p"):

  deg = lambda E: 4*np.pi*np.sqrt(2*m**3*E)*(L/h)**3
  long_term = lambda T: (2*np.pi*h_bar**2/(m*T))**0.5
  Ef = h_bar**2/(2*m)*(6*np.pi**2*N/V)**(2/3)

  data = np.loadtxt("LUT_F32.txt")
  z32 = data[0,:]
  F32 = data[1,:]

  data = np.loadtxt("LUT_F52.txt")
  z52 = data[0,:]
  F52 = data[1,:]

  F52_A = lambda z: (8/(15*np.pi**0.5))*(np.log(z)**2.5)*(1 + 5*np.pi**2/(8*np.log(z)**2) - (5*3/16)*(7*np.pi**4/360)*np.log(z)**-4)

  def Pres(T,V):
    Y = long_term(T)**3*N/V
    z = dame_z(Y)
    if (F52[-1]<=Y):
      return F52_A(z)
    else:
      inf = 0
      sup = len(z52)-1
      while (inf<sup-1):
        med = (inf+sup)//2
        if(z52[med]<z):
          inf = med
        else:
          sup = med
      pend = (F52[inf+1]-F52[inf])/(z52[inf+1]-z52[inf])
      return (pend*(z-z52[inf])+F52[inf])*T/long_term(T)**3

  rhos = ['rho0', 'rho1', 'rho2']
  forma = ['s--', '^--', 'o--']
  V = (12*10*np.array([1/3,1/2,1]))**3
  parametros = ["", "Maruyama/"]
  P_degs = np.zeros((2,3))
  P_degs_FD = np.zeros((2,3))
  for eso in range(2):
    i = 0
    rho = parametros[eso]+"rho0"
    plt.figure()
    if (rho[0] == "M"):
      rhos = ['Maruyama/rho0','Maruyama/rho1','Maruyama/rho2']
      V = (2*1.644*10*np.array([1.0/3.0, 1.0/2.0, 1.0/1.0]))**3
    rango_T = np.linspace(min(Ts), max(Ts), 1000)
    T_min = 0.1
    if (rho[0]=="M"):
      T_min = 0.5
    for rho in rhos:
      data = np.loadtxt(rho+'/presiones_N=1000.txt')
      Ts = data[0,:]
      P = data[1,:]
      P_degs[eso,i] = P[list(Ts).index(T_min)]
      plt.plot(Ts, P*V[i]/N, forma[i])
      i += 1
    P_degs_FD[eso,:] = np.array([Pres(T_min, v) for v in V])
    plt.axis([0, 1.05*max(Ts), 0, 1.05*max(Ts)])
    plt.plot([0, 1.05*max(Ts)], [0, 1.05*max(Ts)], "k-")
    plt.xlabel(r"$T$ [MeV]", fontsize=14)
    plt.ylabel(r"$P/\rho$ [MeV]", fontsize=14)
    plt.grid()
    leyenda = [r"$\rho^*$ = %1.3f "%(qo**3*N/V[k]) for k in range(3)]+["Gas ideal"]
    plt.legend(leyenda, loc=2, fontsize=14)
    if(rho[0]=="M"):
      plt.title("Parametros Maruyama")
    else:
      plt.title("Parametros Dorso")

    plt.figure()
    i = 0
    for rho in rhos:
      data = np.loadtxt(rho+'/temperaturas_N=1000.txt')
      Ts = data[0,:]
      Tv = data[1,:]
      plt.plot(Ts, Tv, forma[i])
      i += 1
    plt.axis([0, 1.05*max(Ts), 0, 1.05*max(Ts)])
    plt.plot([0, 1.05*max(Ts)], [0, 1.05*max(Ts)], "k-")
    plt.xlabel(r"$T$ [MeV]", fontsize=14)
    plt.ylabel(r"$T_V$ [MeV]", fontsize=14)
    plt.grid()
    leyenda = [r"$\rho^*$ = %1.3f "%(qo**3*N/V[k]) for k in range(3)]+["Gas ideal"]
    plt.legend(leyenda, loc=2, fontsize=14)
    if(rho[0]=="M"):
      plt.title("Parametros Maruyama")
    else:
      plt.title("Parametros Dorso")
    plt.show()

  plt.figure()
  plt.plot(P_degs[0,:]/P_degs_FD[0,:], "s--")
  plt.plot(P_degs[1,:]/P_degs_FD[1,:], "o--")
  plt.xticks([0,1,2], [3.375, 1, 0.125])
  plt.xlabel(r"$\rho^*$")
  plt.ylabel(r"$P/P_{FD}$")
  plt.axis([-0.5, 2.5, 1, 3.5])
  plt.grid()
  plt.title("Presión de degeneración")
  plt.legend(["Dorso", "Maruyama"])



deltas_c = pressure.delta_fases
deltas_c.argtypes = [ct.c_voidp, ct.c_voidp, ct.c_voidp, ct.c_longlong,
                    ct.c_voidp, ct.c_voidp, ct.c_float]
deltas_c.restype = ct.c_float

pairs = np.array(list(it.combinations(range(N), 2)), dtype=np.int64)
p = len(pairs)
def deltas(x, p):
  dq = np.zeros(len(pairs), dtype=np.float32)
  dp = np.zeros(len(pairs), dtype=np.float32)
  xp = x.ctypes.data_as(ct.c_voidp)
  pp = p.ctypes.data_as(ct.c_voidp)
  pairsp = pairs.ctypes.data_as(ct.c_voidp)
  dq_p = dq.ctypes.data_as(ct.c_voidp)
  dp_p = dp.ctypes.data_as(ct.c_voidp)
  deltas_c(xp, pp, pairsp, len(pairs), dq_p, dp_p, np.float32(L))
  return dq, dp


def comparar_gas_ideal(T, r = "0", zoom=False):
  n = len(pairs)
  dq = np.zeros(Nreps*n, dtype=np.float32)
  dp = np.zeros(Nreps*n, dtype=np.float32)
  for i in range(Nreps):
    sim_q = np.random.random(3*1000)*L
    sim_p = np.random.normal(0,1,3*1000)*np.sqrt(m*T)*factor_p
    sim_q = np.array(sim_q, dtype=np.float32)
    sim_p = np.array(sim_p, dtype=np.float32)
    dq[n*i:n*(i+1)], dp[n*i:n*(i+1)] = deltas(sim_q, sim_p)
  plt.figure()
  counts,xbins,ybins,image = plt.hist2d(dq, dp, bins=200, norm=LogNorm(), cmap = plt.cm.rainbow)
  new_counts = scipy.ndimage.filters.gaussian_filter(counts, 1)
  plt.colorbar()
  CS = plt.contour(new_counts.transpose(),extent=[xbins[0],xbins[-1],ybins[0],ybins[-1]],
                linewidths=3, colors = "black", levels = np.logspace(1, 2, 2))
  plt.clabel(CS, colors = "black", inline=True, fmt="%d", fontsize=20)
  plt.xlabel(r'$\Delta q$ [fm]', fontsize=14)
  plt.ylabel(r'$\Delta p$ [MeV/c]', fontsize=14)
  plt.title(r'Gas ideal - $\rho^*=%1.3f$ - $T=%1.4fMeV$' %(qo**3*N/V, T), fontsize=16)
  plt.axis([0,np.sqrt(3)*L/2, 0, 4000])
  plt.grid()
  if (zoom):
    plt.axis([0, 4*qo, 0, 4*po*factor_p])
  # DORSO
  rho = "rho"+r
  data = np.loadtxt(rho+"/histograma_2D_%f.txt" %T)
  xbins = data[0,1:]
  ybins = data[1:,0]
  counts = data[1:-1,1:-1]
  Nbins_2D = len(xbins)-1
  dq = []
  dp = []
  for x in range(len(xbins)-1):
    for y in range(len(ybins)-1):
      for i in range(int(counts[x,y])):
        dq.append((xbins[x+1]+xbins[x])*0.5)
        if (rho[0]!="M"):
          dp.append((ybins[y+1]+ybins[y])*0.5*factor_p)
        else:
          dp.append((ybins[y+1]+ybins[y])*0.5)
  plt.figure()
  counts,xbins,ybins,image = plt.hist2d(dq, dp, bins=Nbins_2D, norm=LogNorm(), cmap = plt.cm.rainbow)
  new_counts = scipy.ndimage.filters.gaussian_filter(counts, 1)
  plt.colorbar()
  CS = plt.contour(new_counts.transpose(),extent=[xbins[0],xbins[-1],ybins[0],ybins[-1]],
                linewidths=3, colors = "black", levels = np.logspace(1, 2, 2))
  plt.clabel(CS, colors = "black", inline=True, fmt="%d", fontsize=20)
  plt.xlabel(r'$\Delta q$ [fm]', fontsize=14)
  plt.ylabel(r'$\Delta p$ [MeV/c]', fontsize=14)
  plt.title(r'Dorso - $\rho^*=%1.3f$ - $T=%1.4fMeV$' %(qo**3*N/V, T), fontsize=16)
  plt.axis([0,np.sqrt(3)*L/2, 0, 4000])
  plt.grid()
  if (zoom):
    plt.axis([0, 4*qo, 0, 4*po*factor_p])
  plt.show()
  return dq, dp

if(tipo == 'ea'):
  Nsamp = 1

  files = glob.glob(rho+"/distribucion_10_rep1_*")
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

  Nbins_2D = 200

  dq = np.zeros(p*Nsamp*Nreps)
  dp = np.zeros(p*Nsamp*Nreps)
  for k in range(len(Ts)):
    print("Analizando T=%fMeV"%Ts[k])
    data_q = np.array([], dtype=np.float32)
    data_p = np.array([], dtype=np.float32)
    for j in range(Nreps):
      data_aux = np.loadtxt(rho+"/distribucion_10_rep%d_%f.txt" %(j+1, Ts[k]), dtype=np.float32)
      data_q = np.concatenate([data_q, data_aux[:, 0]])
      data_p = np.concatenate([data_p, data_aux[:, 1]])
    for i in range(Nsamp*Nreps):
      dq[p*i:p*(i+1)], dp[p*i:p*(i+1)] = deltas(data_q[3*N*i:3*N*(i+1)], data_p[3*N*i:3*N*(i+1)])
    counts,xbins,ybins,image = plt.hist2d(dq, dp, bins=Nbins_2D, norm=LogNorm(), cmap = plt.cm.rainbow)
    plt.close()
    data = np.zeros((Nbins_2D+2,Nbins_2D+2))
    data[0,1:] = xbins
    data[1:,0] = ybins
    data[1:-1,1:-1] = counts
    np.savetxt(rho+"/histograma_2D_%f.txt" %Ts[k], data)


if (tipo == "ev"):
  files = glob.glob(rho+"/distribucion_10_rep1_*")
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

  if (rho[0]=="M"):
    seleccion = [0, 3, -1]
    if (rho[-1]=="0"):
      seleccion = [0, 5, -1]
  else:
    seleccion = [0, 4, 8]
    if (rho[-1]=="0"):
      seleccion = [8, 12, 16]

  for k in seleccion:
    data = np.loadtxt(rho+"/histograma_2D_%f.txt" %Ts[k])
    xbins = data[0,1:]
    ybins = data[1:,0]
    counts = data[1:-1,1:-1]
    Nbins_2D = len(xbins)-1
    dq = []
    dp = []
    for x in range(len(xbins)-1):
      for y in range(len(ybins)-1):
        for i in range(int(counts[x,y])):
          dq.append((xbins[x+1]+xbins[x])*0.5)
          if (rho[0]!="M"):
            dp.append((ybins[y+1]+ybins[y])*0.5*factor_p)
          else:
            dp.append((ybins[y+1]+ybins[y])*0.5)
    plt.figure()
    counts,xbins,ybins,image = plt.hist2d(dq, dp, bins=Nbins_2D, norm=LogNorm(), cmap = plt.cm.rainbow)
    new_counts = scipy.ndimage.filters.gaussian_filter(counts, 1)
    plt.colorbar()
    CS = plt.contour(new_counts.transpose(),extent=[xbins[0],xbins[-1],ybins[0],ybins[-1]],
                  linewidths=3, colors = "black", levels = np.logspace(1, 2, 2))
    plt.clabel(CS, colors = "black", inline=True, fmt="%d", fontsize=20)
    plt.xlabel(r'$\Delta q$ [fm]', fontsize=14)
    plt.ylabel(r'$\Delta p$ [MeV/c]', fontsize=14)
    if (rho[0]=="M"):
      plt.title(r'Maruyama - $\rho^*=%1.3f$ - $T=%1.4fMeV$' %(qo**3*N/V, Ts[k]), fontsize=16)
    else:
      plt.title(r'Dorso - $\rho^*=%1.3f$ - $T=%1.4fMeV$' %(qo**3*N/V, Ts[k]), fontsize=16)
    #plt.axis([0, 4*qo, 0, 4*po*factor_p])
  plt.show()


if (tipo == "eg"):
  files = glob.glob(rho+"/distribucion_10_rep1_*")
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

  def contorno_inferior(counts):
    n = len(counts[:,0])
    c = []
    for q in range(n):
      c.append(min(np.where(counts[q,:]!=0))[0])
    c = np.array(c)
    return c

  def contorno_izquierdo(counts):
    n = len(counts[0,:])
    c = []
    for p in range(n):
      zeros = np.where(counts[:,p]!=0)[0]
      if (len(zeros)==0):
        c.append(n+1)
      else:
        c.append(min(zeros))
    c = np.array(c)
    return c

  As = np.zeros_like(Ts)
  plt.figure()
  formas=["o","v","s"]
  for i in range(3):
    rho = rho[:-1]+str(i)
    for k in range(n_temps):
      data = np.loadtxt(rho+"/histograma_2D_%f.txt" %Ts[k])
      xbins = data[0,1:]
      ybins = data[1:,0]
      counts = data[1:-1,1:-1]
      Nbins_2D = len(xbins)-1

      c_inf = contorno_inferior(counts)
      c_izq = contorno_izquierdo(counts)
      m_q = min(list(c_inf).index(min(c_inf)), int(5*qo/(xbins[1]-xbins[0])))
      m_p = min(list(c_izq).index(min(c_izq)), int(5*po/(ybins[1]-ybins[0])))
      As[k] = 4*np.sum(counts[:m_q,:m_p]==0)*(xbins[1]-xbins[0])*(ybins[1]-ybins[0])/(qo*po)
    plt.plot(Ts, As, formas[i]+"-")
  plt.xlabel(r"$T$ [MeV]", fontsize = 14)
  plt.ylabel(r"$A^*$", fontsize = 14)
  plt.legend([r"$\rho^*=3.375$", r"$\rho^*=1.000$", r"$\rho^*=0.125$"], loc=3)
  if (rho[0]=="M"):
    plt.plot([0,Ts[-1]], [39.486138343113048, 39.486138343113048], "k--")
    plt.axis([0,5,15,45])
    plt.title("Maruyama", fontsize = 14)
  else:
    plt.plot([0,Ts[-1]], [26.351538860175374, 26.351538860175374], "k--")
    plt.axis([0,0.5,10,31])
    plt.title("Dorso", fontsize = 14)
  plt.show()




if (tipo[:2] == "gr"):

  gr_c = pressure.Gr
  gr_c.argtypes = [ct.c_voidp, ct.c_int, ct.c_float,
                       ct.c_float, ct.c_voidp, ct.c_int]
  gr_c.restype = ct.c_float

  pairs = np.array(list(it.combinations(range(N), 2)), dtype=np.int64)
  def Gr(x, dr):
    M = int(np.ceil(np.sqrt(3)*L/(2*dr)))
    gr = np.zeros(M, dtype=np.float32)
    xp = x.ctypes.data_as(ct.c_voidp)
    grp = gr.ctypes.data_as(ct.c_voidp)
    gr_c(xp, len(x)//3, dr, L, grp, M)
    return gr

  dr = L/(2*400)
  gr = []
  seleccion = range(8)
  seleccion = [0, int(0.3*len(Ts)), int(0.6*len(Ts)), len(Ts)-1]
  if (rho[0]=="M"):
    Ts_quiero = [0.5, 3, 20]
  if (rho[0]!="M"):
    Ts_quiero = [0.01, 0.1, 1]
  if (len(tipo)==3):
    Ts_quiero = [1E-4, 1E-3, 5E-3]
  seleccion = [list(Ts).index(T) for T in Ts_quiero]
  i = 0
  for k in seleccion:
    data_q = np.array([], dtype=np.float32)
    data_p = np.array([], dtype=np.float32)
    for j in range(Nreps):
      data_aux = np.loadtxt(rho+"/distribucion_10_rep%d_%f.txt" %(j+1, Ts[k]), dtype=np.float32)
      data_q = np.concatenate([data_q, data_aux[:, 0]])
      data_p = np.concatenate([data_p, data_aux[:, 1]])
    gr.append(Gr(data_q[0:3*N], dr)/(Nsamp*Nreps))
    for j in range(1,Nsamp*Nreps):
      q_aux = np.array(data_q[3*N*j:3*N*(j+1)], dtype=np.float32)
      gr[i] += Gr(q_aux, dr)/(Nsamp*Nreps)
    i += 1
  for i in range(len(seleccion)):
    plt.figure()
    plt.plot(dr*np.arange(len(gr[i])), gr[i], "b-")
    plt.plot([0, 0.5*L], [1, 1], "r-")
    plt.xlabel(r"$r$ [fm]", fontsize=16)
    plt.ylabel(r"$g(r)$", fontsize=16)
    plt.grid()
    plt.tick_params(axis='both', which='major', labelsize=14)
    if (rho[0]!="M"):
      plt.title(r"Parametros Dorso - $T=%1.2f MeV$ - $\rho^*=%1.3f $"%(Ts[seleccion[i]], qo**3*N/V), fontsize=16)
      if (len(tipo)==3):
        plt.title(r"Parametros Dorso - $T=%1.4f MeV$ - $\rho^*=%1.3f $"%(Ts[seleccion[i]], qo**3*N/V), fontsize=16)
    if (rho[0]=="M"):
      plt.title(r"Parametros Maruyama - $T=%1.1f MeV$ - $\rho^*=%1.3f $"%(Ts[seleccion[i]], qo**3*N/V), fontsize=16)
    plt.axis([0, 0.5*L, 0, max(gr[i])])
  plt.show()


if (tipo == "evt"):
  colores = ["s--", "v--", "o--"]
  E = np.zeros((len(Ts),3))
  plt.figure()
  for j in range(3):
    rhoj = rho[:-1]+str(j)
    for i in range(len(Ts)):
      data = np.loadtxt(rhoj+"/histograma_100_T=%1.6f.txt" %(Ts[i]))
      E[i,j] = np.sum(data[1,:]*data[0,:])*(data[0,1]-data[0,0])
    plt.plot(Ts, (2/3)*E[:,j]/1000, colores[j])
  plt.plot(Ts, Ts, "r-")
  plt.legend([r"$\rho^*=%1.3f$" %(r) for r in [3.375, 1, 0.125]], loc=4)
  plt.title("Parametros "+(rho[0]=="M")*"Maruyama"+(1-(rho[0]=="M"))*"Dorso")
  plt.grid()
  plt.xlabel(r"$T$ [MeV]", fontsize=16)
  plt.ylabel(r"$\frac{2}{3}\frac{E_{cin}}{N}$ [MeV]", fontsize=16)
  plt.axis([0, max(Ts)*1.05, 0, max(Ts)*1.05])
  plt.show()


if (tipo == "evp"):
  colores = ["s--", "v--", "o--"]
  E = np.zeros((len(Ts),3))
  P = np.zeros((len(Ts),3))
  V = N*(np.array([1.0/3, 1.0/2, 1])*qo*2)**3
  plt.figure()
  for j in range(3):
    rhoj = rho[:-1]+str(j)
    data = np.loadtxt(rhoj+"/presiones_N=1000.txt")
    P_aux = data[1,:]
    P[:,j] = P_aux[[list(data[0,:]).index(T) for T in Ts]]
    for i in range(len(Ts)):
      data = np.loadtxt(rhoj+"/histograma_100_T=%1.6f.txt" %(Ts[i]))
      E[i,j] = np.sum(data[1,:]*data[0,:])*(data[0,1]-data[0,0])
    plt.semilogx(Ts, E[:,j]/(1.5*P[:,j]*V[j]), colores[j])
  plt.plot(Ts, np.zeros_like(Ts)+1, "r-")
  plt.legend([r"$\rho^*=%1.3f$" %(r) for r in [3.375, 1, 0.125]], loc=4)
  plt.title("Parametros "+(rho[0]=="M")*"Maruyama"+(1-(rho[0]=="M"))*"Dorso")
  plt.grid()
  plt.xlabel(r"$T$ [MeV]", fontsize=16)
  plt.ylabel(r"$\frac{2E_{cin}}{3PV}$", fontsize=20)
  plt.axis([0, max(Ts)*1.05, 0.75, 1.05])
  plt.show()
