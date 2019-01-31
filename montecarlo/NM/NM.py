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

Vo = 25.93
r1 = 1.757
r2 = 1.771
p1 = 6.2
p2 = 3.0
d = 3.35
a = 5.0/6.0
rcut = 6

h = h_bar*2*np.pi
N = 10**3
Nsamp = 1
Nreps = 1

tipo = "rep"
caso = 'layers/x1/'
Nbins = 100

nargs = len(sys.argv)
if (nargs >= 2):
  tipo = sys.argv[1]
if (nargs >= 3):
  caso = sys.argv[2]


files = glob.glob(caso+"energias_*")
rhos = [float(f.split("_")[1][:-4]) for f in files]
n_rhos = len(rhos)
rhos = np.sort(rhos)

if(tipo == 't'):
  E_kin = []
  E_nuc = []
  E_pauli = []
  E = []
  Es = np.zeros(n_rhos)

  for k in range(n_rhos):
    data = np.loadtxt(caso+"energias_%f.txt" %(rhos[k]), dtype=np.float32)
    E_kin.append(data[:,0])
    E_nuc.append(data[:,1])
    E_pauli.append(data[:,2])
    E.append(E_kin[-1] + E_nuc[-1] + E_pauli[-1])
    """
    plt.subplot(n_rhos//2, 2, k+1)
    plt.plot(E[-1]/N, "b-")
    plt.plot(E_kin[-1]/N, "r--")
    plt.plot(E_nuc[-1]/N, "k--")
    plt.plot(E_pauli[-1]/N, "g--")
    plt.legend([r"$\rho=%f fm^{-3}$" %rhos[k]])
    """
    Es[k] = np.mean(E[-1])/N
  plt.figure()
  plt.plot(rhos, Es, "o-")
  plt.plot([0.04, 0.04], [-50, -10], "k--")
  plt.plot([0.06, 0.06], [-50, -10], "k--")
  plt.plot([0.095, 0.095], [-50, -10], "k--")
  plt.plot([0.1625, 0.1625], [-50, -10], "k--")
  plt.text(0.01, -40, "ñoqui")
  plt.text(0.045, -47, "s\np\na\ng\nh\ne\nt\nt\ni")
  plt.text(0.063, -20, "lasagna")
  plt.text(0.115, -20, "agujero")
  plt.xlabel(r"$\rho [fm^{-3}]$")
  plt.ylabel(r"$E [MeV]$")
  plt.show()

pressure = ct.CDLL('../pressure.so')

if(tipo == 'e'):

  seleccion = range(n_rhos)

  if (nargs >= 4):
    seleccion = []
    for i in range(3,nargs):
      seleccion.append(int(sys.argv[i]))

  deltasPBC_c = pressure.delta_fases
  deltasPBC_c.argtypes = [ct.c_voidp, ct.c_voidp, ct.c_voidp, ct.c_longlong,
                        ct.c_voidp, ct.c_voidp, ct.c_float]
  deltasPBC_c.restype = ct.c_float

  deltas_c = pressure.delta_fases_sin_PBC
  deltas_c.argtypes = [ct.c_voidp, ct.c_voidp, ct.c_voidp, ct.c_longlong,
                        ct.c_voidp, ct.c_voidp]
  deltas_c.restype = ct.c_float

  g = 4
  n = N//g
  pairs = np.array(list(it.combinations(range(n), 2)), dtype=np.int64)
  p = len(pairs)
  def deltas(x, p):
    dq = np.zeros(len(pairs), dtype=np.float32)
    dp = np.zeros(len(pairs), dtype=np.float32)
    xp = x.ctypes.data_as(ct.c_voidp)
    pp = p.ctypes.data_as(ct.c_voidp)
    pairsp = pairs.ctypes.data_as(ct.c_voidp)
    dq_p = dq.ctypes.data_as(ct.c_voidp)
    dp_p = dp.ctypes.data_as(ct.c_voidp)
    deltas_c(xp, pp, pairsp, len(pairs), dq_p, dp_p)
    return dq, dp

  L = (N/rhos)**(1/3)

  dq = np.ones(g*len(pairs))
  dp = np.ones(g*len(pairs))

  for k in seleccion:
    data_aux = np.loadtxt(caso+"distribucion_%f.txt" %(rhos[k]), dtype=np.float32)
    data_q = data_aux[:, 0]
    data_p = data_aux[:, 1]
    plt.figure()
    for i in range(g):
      q_temp = np.array(data_q[i*3*n:(i+1)*3*n], dtype=np.float32)
      p_temp = np.array(data_p[i*3*n:(i+1)*3*n], dtype=np.float32)
      dq[i*p:(i+1)*p], dp[i*p:(i+1)*p] = deltas(q_temp, p_temp)
    counts, xbins, ybins, image = plt.hist2d(dq, dp, bins=100, norm=LogNorm(), cmap = plt.cm.rainbow)
    plt.colorbar()
    new_counts = scipy.ndimage.filters.gaussian_filter(counts, 1)
    CS = plt.contour(new_counts.transpose(),extent=[xbins[0],xbins[-1],ybins[0],ybins[-1]],
                  linewidths=3, colors = "black", levels = np.logspace(0, 2, 3))
    plt.clabel(CS, colors = "black", inline=True, fmt="%d", fontsize=20)
    plt.xlabel(r'$\Delta q$')
    plt.ylabel(r'$\Delta p$')
    plt.axis([0, xbins[-1], 0, ybins[-1]])
    plt.title(r'$\rho=%1.4f fm^{-3}$' %(rhos[k]))
  plt.show()


if (tipo == "gr"):

  seleccion = range(n_rhos)

  if (nargs >= 4):
    seleccion = []
    for i in range(3,nargs):
      seleccion.append(int(sys.argv[i]))

  gr_c = pressure.Gr
  gr_c.argtypes = [ct.c_voidp, ct.c_int, ct.c_float,
                       ct.c_float, ct.c_voidp, ct.c_int]
  gr_c.restype = ct.c_float


  pairs = np.array(list(it.combinations(range(N), 2)), dtype=np.int64)
  def Gr(x, dr, L):
    M = int(np.ceil(np.sqrt(3)*L/(2*dr)))
    gr = np.zeros(M, dtype=np.float32)
    xp = x.ctypes.data_as(ct.c_voidp)
    grp = gr.ctypes.data_as(ct.c_voidp)
    gr_c(xp, len(x)//3, dr, L, grp, M)
    return gr

  L = (N/rhos)**(1/3)
  dr = L/(2*400)
  gr = []
  i = 0
  for k in seleccion:
    data_aux = np.loadtxt(caso+"distribucion_%f.txt" %(rhos[k]), dtype=np.float32)
    data_q = data_aux[:, 0]
    gr.append(Gr(data_q[0:3*N], dr[k], L[k])/(Nsamp*Nreps))
  i = 0
  for k in seleccion:
    plt.subplot(len(seleccion)//2, 2, i+1)
    plt.plot(dr[i]*np.arange(len(gr[i])), gr[i], "b-")
    plt.plot([0, 0.5*L[i]], [1, 1], "r-")
    plt.xlabel(r"$r$")
    plt.ylabel(r"$g(r)$")
    plt.title(r"$T=0.5 MeV$ $\rho=%f fm^{-3}$" %(rhos[k]))
    plt.axis([0, 0.5*L[k], 0, max(gr[i])])
    i += 1
  plt.show()


if (tipo == "vp"):
  T = 0.5
  data = np.loadtxt(caso+'presiones_N=1000.txt')
  rhos = data[0,:]
  P = data[1,:]
  plt.figure()
  #plt.plot(rhos, T*rhos, "k-")
  plt.plot(rhos, P, "bo--")
  plt.xlabel(r"$\rho$  [$fm^{-3}$]")
  plt.ylabel(r"$P$   [$MeV$ $fm^{-3}$]")
  plt.plot([0.04, 0.04], [-50, 100], "k-")
  plt.plot([0.06, 0.06], [-50, 100], "k-")
  plt.plot([0.095, 0.095], [-50, 100], "k-")
  plt.plot([0.1625, 0.1625], [-50, 100], "k-")
  plt.text(0.01, 10, "ñoqui")
  plt.text(0.045, 17, "s\np\na\ng\nh\ne\nt\nt\ni")
  plt.text(0.063, 20, "lasagna")
  plt.text(0.115, 20, "agujero")
  plt.grid()
  plt.axis([0, 0.25, -50, 100])
  plt.show()

if (tipo == "ap"):
  pressure_pauli_c = pressure.pressure_pauli_PBC
  pressure_pauli_c.argtypes = [ct.c_voidp, ct.c_voidp, ct.c_voidp, ct.c_longlong, ct.c_float,
                       ct.c_float, ct.c_float, ct.c_float, ct.c_voidp, ct.c_voidp, ct.c_float]
  pressure_pauli_c.restype = ct.c_float

  n = N//4

  pairs1 = np.array(list(it.combinations(range(n), 2)), dtype=np.int64)
  pairs2 = np.array(list(it.combinations(range(n,2*n), 2)), dtype=np.int64)
  pairs3 = np.array(list(it.combinations(range(2*n,3*n), 2)), dtype=np.int64)
  pairs4 = np.array(list(it.combinations(range(3*n,4*n), 2)), dtype=np.int64)
  pairs_pauli = np.concatenate([pairs1, pairs2, pairs3, pairs4])
  def pressure_pauli(x, p, L):
    energ = 0
    force = np.zeros_like(x, dtype=np.float32)
    gorce = np.zeros_like(x, dtype=np.float32)
    xp = x.ctypes.data_as(ct.c_voidp)
    pp = p.ctypes.data_as(ct.c_voidp)
    pairsp = pairs_pauli.ctypes.data_as(ct.c_voidp)
    forcep = force.ctypes.data_as(ct.c_voidp)
    gorcep = gorce.ctypes.data_as(ct.c_voidp)
    qF = pressure_pauli_c(xp, pp, pairsp, len(pairs_pauli), D, qo, po, scut, forcep, gorcep, L)
    return force, gorce, qF

  pressure_qcnm_c = pressure.pressure_qcnm_layers
  pressure_qcnm_c.argtypes = [ct.c_voidp, ct.c_voidp, ct.c_longlong, ct.c_float, ct.c_float, ct.c_float,
              ct.c_float, ct.c_float, ct.c_float, ct.c_float, ct.c_float, ct.c_voidp, ct.c_float, ct.c_int]
  pressure_qcnm_c.restype = ct.c_float

  pairs_nuc = np.array(list(it.combinations(range(N), 2)), dtype=np.int64)
  def pressure_nuc(x, L):
    energ = 0
    force = np.zeros_like(x, dtype=np.float32)
    xp = x.ctypes.data_as(ct.c_voidp)
    pairsp = pairs_nuc.ctypes.data_as(ct.c_voidp)
    forcep = force.ctypes.data_as(ct.c_voidp)
    qF = pressure_qcnm_c(xp,  pairsp, len(pairs_nuc), Vo, r1, p1, r2, p2, d, a, rcut, forcep, L, 1)
    return force, qF

  V = N/rhos
  L = V**(1.0/3)
  qF = np.zeros_like(rhos)
  qp = np.zeros_like(rhos)
  P = np.zeros_like(rhos)
  for k in range(n_rhos):
    data = np.loadtxt(caso+'distribucion_%f.txt' %(rhos[k]), dtype=np.float32)
    data_q = np.array(data[:, 0], dtype=np.float32)
    data_p = np.array(data[:, 1], dtype=np.float32)

    fuerzas, guerzas, qF_aux = pressure_pauli(data_q, data_p, L[k])
    qp[k] += np.sum(data_p*(data_p/m-guerzas))/(3*V[k])
    qF[k] += qF_aux/(3*V[k])
    fuerzas, qF_aux = pressure_nuc(data_q, L[k])
    qF[k] += qF_aux/(3*V[k])

    P[k] = qF[k] + qp[k]
  np.savetxt(caso+"presiones_N=1000.txt", [rhos, P])


if (tipo=="vmd"):
  for k in range(n_rhos):
    new_filename = "sistema_{0}.lammpstrj".format(rhos[k])
    data = np.loadtxt(caso+"distribucion_%f.txt" %(rhos[k]), dtype=np.float32)
    L = (1000/rhos[k])**(1.0/3.0)
    header1 = "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n{0}\n".format(1000)
    header2 = "ITEM: BOX BOUNDS pp pp pp\n{0} {1}\n{2} {3}\n{4} {5}\n".format(0, L, 0, L, 0, L)
    header3 = "ITEM: ATOMS type x y z\n"
    f = open(new_filename, "w")
    f.write(header1 + header2 + header3)
    for i in range(len(data[:,0])//3):
      f.write("{0} {1} {2} {3}\n".format(i//250, data[3*i, 0], data[3*i+1, 0], data[3*i+2, 0]))
    f.close()
