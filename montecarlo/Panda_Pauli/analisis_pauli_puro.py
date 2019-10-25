import numpy as np
import matplotlib.pylab as plt
import matplotlib.animation as animation
import matplotlib.gridspec as gridspec
import sys
from GA import *
from mediador import *
plt.ion()

if(sys.argv[1]=="pauli"):
	formato = "Pauli_puro/B3/termo_B3_fija_%.3f.txt"
	rhos = np.array([0.04, 0.08, 0.12, 0.16, 0.20, 0.24])

	n_rhos = len(rhos)

	Ek = [[] for rho in rhos]
	Ep = [[] for rho in rhos]
	for r,rho in zip(range(len(rhos)),rhos):
		data = np.loadtxt(formato %rho)
		Ek[r] = data[:,1]
		Ep[r] = data[:,3]
	Ts = data[:,0]
	Ek = np.array(Ek)
	Ep = np.array(Ep)

	Ek_fermi = kinetic_energy(rhos/4, Ts)
	plt.figure()
	plt.plot(rhos, Ek[:, -1], "ko--")
	plt.plot(rhos, 0.6*Ef_rho*(rhos/4)**(2/3), "k-")
	plt.xlabel(r"$\rho$ [fm$^{-3}$]", fontsize=18)
	plt.ylabel("$E_k$ [MeV]", fontsize=18)
	plt.axis([0.03, 0.25, 0, 30])
	plt.xticks(rhos)
	plt.yticks(np.concatenate([Ek[:,-1], 0.6*Ef_rho*(rhos/4)**(2/3)]))
	plt.grid()

	plt.figure()
	colores = ["y", "r", "b", "k", "c", "g"]
	for i in range(n_rhos):
		plt.plot(Ts[:-20], Ek[i,:-20], colores[i]+"o--")
	for i in range(n_rhos):
		plt.plot(Ts[:-20], Ek_fermi[i,:-20], colores[i]+"-", linewidth=5)
		plt.text(-0.12, Ek_fermi[i,-20], r"$%.2f$fm$^{-3}$" %rhos[i], color=colores[i], fontsize=14)
	plt.axis([-0.15, 4, 5, 35])
	plt.xlabel("$T$ [MeV]", fontsize=18)
	plt.ylabel("$E_k$ [MeV]", fontsize=18)
	plt.grid()
	plt.show()
elif (sys.argv[1]=="pauli_red"):
	redes = sys.argv[2:]
	rhos = np.array([0.150, 0.155, 0.158, 0.160, 0.162, 0.165, 0.167, 0.170])
	formato1 = "Pauli_puro/B1/termo_maruyama_dqx0.0_%.3f.txt"
	formato2 = "Pauli_puro/B2/termo_B2_fija_%.3f.txt"
	formato3 = "Pauli_puro/B3/termo_B3_fija_%.3f.txt"
	n_rhos = len(rhos)
	Ek = np.zeros((len(redes),n_rhos))
	Ep = np.zeros((len(redes),n_rhos))
	leyenda = []
	for i,red in zip(range(len(redes)), redes):
		formato = formato1*(red=="B1") + formato2*(red=="B2") + formato3*(red=="B3")
		for r,rho in zip(range(n_rhos),rhos):
			data = np.loadtxt(formato %rho)
			Ek[i,r] = data[-1,1]
			Ep[i,r] = data[-1,3]
		To = data[-1,0]
		Ek = np.array(Ek)
		Ep = np.array(Ep)
		plt.figure(1)
		plt.plot(rhos, Ek[i,:], "o--")
		plt.figure(2)
		plt.plot(rhos, Ek[i,:] + Ep[i,:], "o--")
		leyenda.append(red)
	plt.figure(1)
	plt.plot(rhos, 0.6*Ef_rho*(rhos/4)**(2/3), "k-")
	plt.xlabel(r"$\rho$ [fm$^{-3}$]", fontsize=18)
	plt.ylabel("$E_k$ [MeV]", fontsize=18)
	plt.xlim(0.150, 0.170)
	plt.legend(leyenda+["Fermi"], loc="lower right")
	plt.grid()
	plt.figure(2)
	plt.xlabel(r"$\rho$ [fm$^{-3}$]", fontsize=18)
	plt.ylabel("$E_k + E_p$ [MeV]", fontsize=18)
	plt.xlim(0.150, 0.170)
	plt.legend(leyenda, loc="lower right")
	plt.grid()
	plt.show()
elif (sys.argv[1]=="calorica"):
	folder = sys.argv[2]
	rhos = np.array([0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20])
	n_rhos = len(rhos)

	Ek = [[] for rho in rhos]
	Ep = [[] for rho in rhos]
	En = [[] for rho in rhos]
	plt.figure(1)
	leyenda = []
	for r,rho in zip(range(len(rhos)),rhos):
		data = np.loadtxt("%s/termo_QCNM_%.3f.txt" %(folder, rho))
		Ts =  data[:,0]
		Ek[r] = data[:,1]
		En[r] = data[:,2]
		Ep[r] = data[:,3]
		E = Ek[r] + En[r] + Ep[r]
		plt.plot(Ts, E, "-", lw=3)
		leyenda.append(r"$\rho=%.2f$fm$^{-3}$" %rho)
		#plt.text(np.max(Ts), np.max(E), r"$%.2f$fm$^{-3}$" %rho, fontsize=10)
	Ts = data[:,0]
	Ek = np.array(Ek)
	En = np.array(En)
	Ep = np.array(Ep)
	plt.xlabel(r"$T$ [MeV]", fontsize=18)
	plt.ylabel("$E$ [MeV]", fontsize=18)
	plt.xlim(0, 4)
	plt.legend(leyenda, loc="upper left")
	plt.grid()
	plt.show()
elif (sys.argv[1]=="Evsrho"):
	folder = sys.argv[2]
	rhos = np.array([0.04, 0.05, 0.06, 0.08, 0.10, 0.12, 0.14, 0.15, 0.16, 0.17, 0.18, 0.20])
	n_rhos = len(rhos)
	Ts = np.array([0.1, 0.5, 1, 2, 4])
	n_ts = len(Ts)

	E = np.zeros((n_rhos, n_ts))
	Ek = np.zeros((n_rhos, n_ts))
	En = np.zeros((n_rhos, n_ts))
	Ep = np.zeros((n_rhos, n_ts))
	plt.figure(1)
	for r,rho in zip(range(len(rhos)),rhos):
		#data = np.loadtxt("%s/termo_QCNM_%.3f.txt" %(folder, rho))
		data = np.loadtxt("%s/termo_%.3f.txt" %(folder, rho))
		idxs = [np.where(np.abs(T - data[:,0]) < 0.001)[0][0] for T in Ts]
		Ek[r,:] = data[:,1][idxs]
		En[r,:] = data[:,2][idxs]
		Ep[r,:] = data[:,3][idxs]
		E[r,:] = (Ek + En + Ep)[r,:]
	leyenda = []
	for t in range(n_ts):
		plt.plot(rhos, E[:,t], "o-")
		leyenda.append(r"$T=%.2f$MeV" %Ts[t])
	plt.xlabel(r"$\rho$ [fm$^{-3}$]", fontsize=18)
	plt.ylabel("$E$ [MeV]", fontsize=18)
	plt.legend(leyenda, loc="upper center")
	plt.grid()
	plt.show()
