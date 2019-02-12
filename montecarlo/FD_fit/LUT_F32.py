import numpy as np
import matplotlib.pylab as plt


C = np.array([0]+[(-1)**k/(k+1)**1.5 for k in range(20)])

F32_S = lambda z: np.sum(C*z**np.arange(21))

def Wynn(c, z, p):
  E = np.zeros((p+1, p+2))
  E[0,0] = c[0]
  for i in range(1,p+1):
    E[i, 0] = E[i-1, 0] + c[i]*z**i
  for m in range(p):
    for n in range(p-m):
      E[n, m+1] = E[n+1, m-1] + 1/(E[n+1, m]-E[n, m])
  return E[0,p]

F32_W = lambda z: Wynn(C, z, 20)

F32_A = lambda z: (4/(3*np.pi**0.5))*(np.log(z)**1.5)*(1 + np.pi**2/(8*np.log(z)**2) + (7*np.pi**4/640)*np.log(z)**-4)

F32_APB = lambda z: (4/(3*np.pi**0.5))*(np.log(z)**1.5)

error_rel = lambda a,b: 2*np.abs((a-b)/(a+b))

z1 = np.linspace(0, 0.7, 701)
z2 = np.linspace(0.7001, 20, 2000)
z3 = np.logspace(1.305, 5, 1000)

F1 = [F32_S(z) for z in z1]
F2 = [F32_W(z) for z in z2]
F3 = F32_A(z3)

z = np.concatenate([z1, z2, z3])
F = np.concatenate([F1, F2, F3])

np.savetxt("LUT_F32.txt",[z, F])


""" Ejemplo con la f1(z) = log(z+1)
c_log = [0]+[(-1)**(n+1)/n for n in range(1,21)]
zs = np.linspace(1, 20)
log_2 = [Wynn(c_log, z, 2) for z in zs]
log_6 = [Wynn(c_log, z, 6) for z in zs]
log_10 = [Wynn(c_log, z, 10) for z in zs]
log_20 = [Wynn(c_log, z, 20) for z in zs]
exacto = np.log(zs+1)
plt.figure()
plt.loglog(zs, error_rel(log_2, exacto), "c-")
plt.loglog(zs, error_rel(log_6, exacto), "g-")
plt.loglog(zs, error_rel(log_10, exacto), "r-")
plt.loglog(zs, error_rel(log_20, exacto), "b-")
plt.axis([1, 27, 1E-16, 1])
plt.text(21, error_rel(log_2[-1], exacto[-1])*0.6, "p=1")
plt.text(21, error_rel(log_6[-1], exacto[-1])*0.7, "p=3")
plt.text(21, error_rel(log_10[-1], exacto[-1])*0.7, "p=5")
plt.text(21, error_rel(log_20[-1], exacto[-1])*0.7, "p=10")
plt.title("Error del algoritmo de Wynn al computar log(z+1)")
plt.xlabel("z")
plt.ylabel("Error relativo")
plt.show()

plt.figure()
plt.plot(zs, log_2, "c--")
plt.plot(zs, log_6, "g--")
plt.plot(zs, log_10, "r--")
plt.plot(zs, exacto, "k-")
plt.axis([1, 22, 0.6, 3.2])
plt.text(20.1, log_2[-1], "p=1")
plt.text(20.1, log_6[-1], "p=3")
plt.text(20.1, log_10[-1]-0.02, "p=5")
plt.text(20.1, exacto[-1]+0.02, "Exacto")
plt.title("Computos de log(z+1)")
plt.xlabel("z")
plt.ylabel("log(z+1)")
plt.show()
"""

""" Plot de la F32
z1 = np.linspace(0, 0.65, 651)
z2 = np.linspace(0.6501, 20, 2000)
z3 = np.logspace(1.305, 4, 2000)
zs = np.concatenate([z1,z2,z3])
F1 = [F32_S(z) for z in z1]
F2 = [F32_W(z) for z in z2]
F3 = F32_A(z3)
F = np.concatenate([F1,F2,F3])
plt.loglog(zs, F, "b-")
plt.plot([0.65, 0.65], [1E-3, 1E2], "k--")
plt.plot([20, 20], [1E-3, 1E2], "k--")
plt.axis([1E-3, 1E4, 1E-3, 1E2])
plt.text(0.005, .1, "Sumas\nparciales", fontsize=16)
plt.text(2, .15, "Wynn", fontsize=16)
plt.text(100, .15, "Sommerfeld", fontsize=16)
plt.xlabel("z", fontsize=14)
plt.ylabel(r"$f_{3/2}$(z)", fontsize=14)
plt.show()
"""
