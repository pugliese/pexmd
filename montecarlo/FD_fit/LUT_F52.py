import numpy as np
import matplotlib.pylab as plt


C = np.array([0]+[(-1)**k/(k+1)**2.5 for k in range(20)])

F52_S = lambda z: np.sum(C*z**np.arange(21))

def Wynn(c, z, p):
  E = np.zeros((p+1, p+2))
  E[0,0] = c[0]
  for i in range(1,p+1):
    E[i, 0] = E[i-1, 0] + c[i]*z**i
  for m in range(p):
    for n in range(p-m):
      E[n, m+1] = E[n+1, m-1] + 1/(E[n+1, m]-E[n, m])
  return E[0,p]


F52_W = lambda z: Wynn(C, z, 20)

F52_A = lambda z: (8/(15*np.pi**0.5))*(np.log(z)**2.5)*(1 + 5*np.pi**2/(8*np.log(z)**2) - (5*3/16)*(7*np.pi**4/360)*np.log(z)**-4)

F52_APB = lambda z: (8/(15*np.pi**0.5))*(np.log(z)**2.5)

error_rel = lambda a,b: 2*np.abs(a-b)/(a+b)

z1 = np.linspace(0, 0.75, 751)
z2 = np.linspace(0.751, 20, 2000)
z3 = np.logspace(1.305, 5, 1000)

F1 = [F52_S(z) for z in z1]
F2 = [F52_W(z) for z in z2]
F3 = F52_A(z3)

z = np.concatenate([z1, z2, z3])
F = np.concatenate([F1, F2, F3])

np.savetxt("LUT_F52.txt",[z, F])
