import numpy as np
import matplotlib.pyplot as plt

V0 = 34e6           #[eV]
dx = 17e-15         #[m]
c = 3e8
m = 3.727e9/c**2    #[eV]

hbar = 6.582e-16    #[eVs]

n = 1000
eps = 1e-5          #Small tolerance to make sure E =! 0 & V0

E_arr = np.linspace(0+eps,V0-eps,n)

T = lambda E: 1./(1. + V0**2/(4.*E*(V0-E))*(np.sinh(dx/hbar*np.sqrt(2.*m*(V0-E))))**2)

plt.plot(E_arr/V0,T(E_arr))
plt.xlabel('$E/V_0$'); plt.ylabel('T')
plt.title('Transmission as a function of energy')
plt.show()

E_alfa_1 = 4.08e6       #232Th, [eV]
E_alfa_2 = 9.85e6       #218Th, [eV]
E_alfa_3 = 4.8e6        #226Ra, [eV]
E_alfa_4 = 9.93e6       #226Ra with v = 0.073c, [eV]

print T(E_alfa_1), T(E_alfa_2), T(E_alfa_3), T(E_alfa_4)
print 1./T(E_alfa_3)
print T(E_alfa_1)/T(E_alfa_2)
