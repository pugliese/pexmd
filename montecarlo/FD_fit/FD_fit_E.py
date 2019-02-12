import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
import scipy.ndimage
import itertools as it
import ctypes as ct
import sys
import glob
plt.ion()


tipo = "rep"
rho = 'rho0'
Nbins = 100
Nreps = 1
Nsamp = 2

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
  po = 2.067
  D = 34.32*(h_bar/(po*qo))**3

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
    while (Y*1E-2 < np.abs(F32_A(z_med)-Y)):
      z_med = (z_inf+z_sup)/2
      if(F32_A(z_med)<Y):
        z_inf = z_med
      else:
        z_sup = z_med
    print("Fin")
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

if (tipo == "f&v"):

  MB = lambda x, T: N*2*np.sqrt(x/np.pi)*np.exp(-x/T)/(T**1.5)
  FD = lambda x, mu, T: deg(x)/(np.exp((x-mu)/T)+1)

  #seleccion = range(n_temps)

  if (rho[-1] == "0"):
    seleccion = [8,10,12,13,14,15,16, 20, -1]
  if (rho[-1] == "1"):
    seleccion = [0,1,2,3,4,5,6, 8, 10]
  if (rho[-1] == "2"):
    seleccion = range(9)
  if (rho[0] == "M"):
    seleccion = range(9)

  mus = np.zeros(len(seleccion))
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
      #plt.subplots_adjust(left  = 0.125, right = 0.9, bottom = 0.1, top = 0.9, wspace = 0.2, hspace = 0.2)
    plt.plot(E,  ns/1000, "ko")
    #plt.text(1, 1700, "T=%f" %(Ts[k]))
    rango = np.linspace(0, E[-1], 10000)
    exacto_MB = MB(rango, Ts[k])/1000
    plt.plot(rango, exacto_MB, "r-")
    print(dame_z(long_term(Ts[k])**3*N/V), long_term(Ts[k])**3*N/V)
    mus[i] = np.log(dame_z(long_term(Ts[k])**3*N/V))*Ts[k]
    exacto_FD = FD(rango, mus[i], Ts[k])/1000
    plt.plot(rango, exacto_FD, "k--")
    #plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.tick_params(axis='both', which='major', labelsize=14)

    plt.axis([0, max(E), 0, 1.05*max(exacto_MB)])
    plt.text((min(E)+ 4*max(E))/6, 0.5*max(exacto_MB), "T=%1.3fMeV" %(Ts[k]), fontsize=16)
    #plt.axis([0, 2, 0, 2500])
    if (i==7):
      plt.xlabel(r"$E$ [$MeV$]", fontsize=18)
    if (i==3):
      plt.ylabel(r"$f(E)$ [$GeV^{-1}$]", fontsize=18)
    if (i==4):
      plt.legend(["Data", "Boltzmann", "Fermi"], fontsize=20)
    i+=1
  plt.show()
  mus = np.array([np.log(dame_z(long_term(Ts[k])**3*N/V))*T for T in Ts])

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
  forma = ['sr--', '^b--', 'og--']
  V = (12*10*np.array([1/3,1/2,1]))**3
  i = 0
  plt.figure()
  if (rho[0] == "M"):
    rhos = ['Maruyama/rho0','Maruyama/rho1','Maruyama/rho2']
    V = (2*1.644*10*np.array([1.0/3.0, 1.0/2.0, 1.0/1.0]))**3
  rango_T = np.linspace(min(Ts), max(Ts), 1000)
  for rho in rhos:
    data = np.loadtxt(rho+'/presiones_N=1000.txt')
    Ts = data[0,:]
    P = data[1,:]
    print(P)
    plt.plot(Ts, P*V[i]/N, forma[i])
    #plt.plot(rango_T, [Pres(T,V[i])*V[i]/N for T in rango_T], forma[i][1:-1])
    i += 1
  plt.plot(Ts, Ts, "k-")
  plt.xlabel(r"$T$")
  plt.ylabel(r"$P/\rho$")
  leyenda = [r"$\rho = %f fm^{-3}$"%(N/V[k]) for k in range(3)]
  plt.legend(leyenda, loc=2)
  #plt.legend([, r"FD $\rho = 0.015625fm^{-3}$", r"$\rho = 0.00462963fm^{-3}$", r"FD $\rho = 0.00462963fm^{-3}$", r"$\rho = 0.0005787fm^{-3}$", r"FD $\rho = 0.0005787fm^{-3}$", "Boltzmann"], loc=9)
  plt.show()


if(tipo == 'e'):

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
    deltas_c(xp, pp, pairsp, len(pairs), dq_p, dp_p, L)
    return dq, dp

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
    seleccion = [0, 5, -1]
  else:
    seleccion = [0, 4, 8]

  dq = np.zeros(p*Nsamp*Nreps)
  dp = np.zeros(p*Nsamp*Nreps)
  for k in seleccion:
    data_q = np.array([], dtype=np.float32)
    data_p = np.array([], dtype=np.float32)
    for j in range(Nreps):
      data_aux = np.loadtxt(rho+"/distribucion_10_rep%d_%f.txt" %(j+1, Ts[k]), dtype=np.float32)
      data_q = np.concatenate([data_q, data_aux[:, 0]])
      data_p = np.concatenate([data_p, data_aux[:, 1]])
    plt.figure()
    for i in range(Nsamp*Nreps):
      dq[p*i:p*(i+1)], dp[p*i:p*(i+1)] = deltas(data_q[3*N*i:3*N*(i+1)], data_p[3*N*i:3*N*(i+1)])
    """
    heatmap, xedges, yedges = np.histogram2d(dq, dp, bins=100)
    dQ, dP = np.meshgrid((xedges[1:]+xedges[:-1])/2, (yedges[1:]+yedges[:-1])/2)
    plt.hist2d(dq, dp, bins=100)
    plt.contour(dP, dQ, heatmap, levels = (np.arange(0,5)*np.max(heatmap)/6+10))
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    plt.clf()
    plt.imshow(heatmap.T, extent=extent, origin='lower')
    plt.cm.rainbow
    """
    counts,xbins,ybins,image = plt.hist2d(dq,dp,bins=100, norm=LogNorm(), cmap = plt.cm.rainbow)
    plt.colorbar()
    new_counts = scipy.ndimage.filters.gaussian_filter(counts, 1)
    CS = plt.contour(new_counts.transpose(),extent=[xbins[0],xbins[-1],ybins[0],ybins[-1]],
                  linewidths=3, colors = "black", levels = np.logspace(1, 3, 3))
    plt.clabel(CS, colors = "black", inline=True, fmt="%d", fontsize=20)
    plt.xlabel(r'$\Delta q$')
    plt.ylabel(r'$\Delta p$')
    plt.title(r'$T=%1.4fMeV$' %(Ts[k]))
  plt.show()


if (tipo == "gr"):

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
      gr[i] += Gr(data_q[3*N*j:3*N*(j+1)], dr)/(Nsamp*Nreps)
    i += 1
  for i in range(len(seleccion)):
    plt.subplot(4, 2, i+1)
    plt.plot(dr*np.arange(len(gr[i])), gr[i], "b-")
    plt.plot([0, 0.5*L], [1, 1], "r-")
    plt.xlabel(r"$r$")
    plt.ylabel(r"$g(r)$")
    plt.title(r"$T=%f MeV$ $\rho=%f fm^{-3}$"%(Ts[seleccion[i]], N/V))
    plt.axis([0, 0.5*L, 0, max(gr[i])])
  plt.show()


"""
if (tipo == "v"):
  plt.figure()
  for k in range(n_temps):
    E[k,:] = (bins[k,1:] + bins[k,:-1])/2
    if (n_temps > 1):
      plt.subplot(n_temps//2, 2, k+1)
    plt.xlabel(r"$E$")
    plt.plot(E[k,:],  ns[k,:], "ko")
    plt.text((min(E[k,:])+ 5*max(E[k,:]))/6, 0.9*max(ns[k,:]), "T=%f" %(Ts[k]))
  plt.figure()
  for k in range(n_temps):
    q[k,:] = (bins_q[k,1:] + bins_q[k,:-1])/2
    if (n_temps > 1):
      plt.subplot(n_temps/3, 3, k+1)
    plt.xlabel(r"$q$")
    plt.plot(q[k,:],  ns_q[k,:], "ko")
    plt.axis([0, 60, 0, 100])
    plt.text((min(q[k,:])+ 5*max(q[k,:]))/6, 1.25*max(ns_q[k,:]), "T=%f" %(Ts[k]))
  plt.show()

if (tipo == "f" or tipo == 'f&v'):
  #MB = lambda x, A, mu, T: deg(x)*np.exp(-x/T)*(N/V)*(2*np.pi*h_bar**2/(m*T))**1.5
  MB = lambda x, mu, T: mu*2*np.sqrt(x/np.pi)*np.exp(-x/T)/(T**1.5)
  FD = lambda x, mu, T: deg(x)/(np.exp((x-mu)/T)+1)
  As = np.zeros(n_temps)
  mus = np.zeros(n_temps)
  Ts_aj = np.zeros(n_temps)
  mus_MB = np.zeros(n_temps)
  Ts_aj_MB = np.zeros(n_temps)
  Z = np.zeros(n_temps)
  for k in range(n_temps):
    Z_k = lambda z: F32(z) - long_term(Ts[k])**3*N/V
    #Z[k] = sc.root(Z_k, 1).x
    Z[k] = raiz(Z_k, 0.5, 2)
    E[k,:] = (bins[k,1:] + bins[k,:-1])/2
    params, coso = sc.curve_fit(FD, E[k,:], ns[k,:], [0.1, Ts[k]], bounds = ([-np.inf,0],[np.inf,np.inf]))
    mus[k] = params[0]
    Ts_aj[k] = params[1]
    params, coso = sc.curve_fit(MB, E[k,:], ns[k,:], [100, Ts[k]], bounds = ([-np.inf,0],[np.inf,np.inf]))
    mus_MB[k] = params[0]
    Ts_aj_MB[k] = params[1]

  if(tipo == 'f&v'):
    plt.figure()
    for k in range(n_temps):
      if (n_temps > 1):
        plt.subplot(n_temps//3, 3, k+1)
      plt.xlabel(r"$E$")
      plt.plot(E[k,:],  ns[k,:], "ko")
      rango = np.linspace(0, E[k,-1], 1000)
      exacto_FD = FD(rango, min(Ef,np.log(Z[k])*Ts[k]), Ts[k])
      plt.plot(rango, exacto_FD, "k--")
      exacto_MB = MB(rango, 13**3, Ts[k])
      plt.plot(rango, exacto_MB, "r-")
      plt.text((min(E[k,:])+ 2*max(E[k,:]))/3, 0.9*max(ns[k,:]), "T=%fMeV" %(Ts[k]))
    plt.figure()
    plt.plot(Ts, Ts_aj_MB, "bo--")
    plt.plot(Ts, Ts, "r-")
    plt.xlabel(r"$T$")
    plt.ylabel(r"$T_{aj}$")
    plt.figure()
    #plt.semilogx(Ts, mus, "bo--")
    plt.plot(Ts, mus, "bo--")
    plt.xlabel(r"$T$")
    plt.ylabel(r"$\mu_{FD}$")
    print(Ts_aj/(Ts*m))
    plt.show()
    np.savetxt('data_'+rho+'.txt', [Ts, Ts_aj, mus, As])
"""
