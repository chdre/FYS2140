import numpy as np
import matplotlib.pyplot as plt

x0 = 5          #[fm]
a = 1           #[fm]
k = 1.38        #[fm]^(-1)

A = (1./(2*np.pi*a**2))**(1./4)

def psi(x):
    abs_square = A**2*np.exp(-(x-x0)**2/(2.*a**2))
    return abs_square

n = 100
x = np.linspace(0,2*x0,n)

plt.plot(x,psi(x))
plt.xlabel('x [fm]'); plt.ylabel('$|\\Psi(x,0)|^2$')
plt.title('Probability density')
plt.show()

hbar = 6.582e-16
print (k*1e15)**2*hbar**2/(2*3.727e9/(3e8)**2)
