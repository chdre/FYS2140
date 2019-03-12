import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sparse
import scipy.sparse.linalg
import scipy as sp

x0 = 5.0            #[fm]
R = 7.3             #[fm]
deltax = 17          #[fm], length of the potential
k = 1.38            #[fm]
a = 1.0             #[fm]
V0 = 34.0           #[MeV]
c = 3e5             #[fm/as]
m = 3.727e3/c**2    #[MeV]
hbar = 6.582e-4     #[MeVas]
Zalpha = 2
ZD = 86
ke = 1.44           #[MeVfm]

#Constants
c1 = -1j*hbar/(2.*m)
c2 = 1j/hbar

#Setup
T = 0.01           #Total run time
dt = 1e-5
t = 0
n = 1000       #Length of the run
dx = 1e-3

#Functions
psi_init = 1./(2*np.pi*a**2)**(1./4)*np.exp(1j*k*x0)
psi0 = lambda x: 1./(2*np.pi*a**2)**(1./4)*np.exp(-(x-x0)**2/(4.*a**2))*np.exp(1j*k*x)

#Arrays
x = np.linspace(x0,2*R+deltax,n)    #Positions

#Calculating second derivative
psi = psi0(x)

#Matrix calculations
I = sparse.identity(n)
data = np.ones((3,n))
data[1] = -2*data[1]
diags = [-1,0,1]
D = c1*sparse.spdiags(data, diags, n, n)/dx**2

V_data = np.zeros(n)

box_potential = True

if box_potential == True:
    for i in range(n):
        if R <= x[i] <= R+deltax:
            V_data[i] = V0
        else:
            V_data[i] = 0
else:
    for i in range(n):
        if 0 <= x[i] < R:
            V_data[i] = 0
        elif R <= x[i]:
            V_data[i] = Zalpha*ZD*ke/x[i]

V_diags = [0]
V = c2*sparse.spdiags(V_data, V_diags, n, n)

plt.ion()
fig = plt.figure()
ax = fig.add_subplot(111)
line, = ax.plot(x, abs(psi)**2, label='$|\Psi(x,t)|^2$')

ax.plot(x,V_data/np.max(V_data), label='$V(x)/V_0$')
fig.suptitle("Solution of Schrodinger's equation")
ax.set_xlabel('$x$ [fm]')
ax.set_ylabel('$|\Psi(x, t)|^2$ [1/fm] and $V(x)/V_0$ [MeV]')
ax.legend(loc='best')
plt.draw()

while t < T:
    A = (I - 0.5*dt*(D + V))
    b = (I + 0.5*dt*(D + V))*psi

    psi = sparse.linalg.spsolve(A,b)

    t += dt

    line.set_ydata(abs(psi)**2)
    plt.draw()
    plt.pause(0.01)

plt.ioff()
plt.show()













#
