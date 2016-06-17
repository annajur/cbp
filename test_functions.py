import numpy as np
from math import exp

import matplotlib.pyplot as plt

# density:

M_bar = 2.0*10**10 #M_sun
x_0 = 1.49 #kpc
y_0 = 0.58 #kpc
z_0 = 0.40 #kpc
q = 0.6
rho_0 = 1 #!!!!!!!!!!!!!!!!!!!

def r_1(x, y, z):
    return((((x/x_0)**2 + (y/y_0)**2)**2 + (z/z_0)**4)**0.25)

def r_2(x, y, z):
    return(((q**2*(x**2 + y**2) + z**2)/z_0**2)**0.5)

def rho(r_1, r_2):
    lrho = []
    for i in range(len(r_1)):
        lrho.append(rho_0*(exp(-r_1[i]**2/2.) + r_2[i]**(-1.85)*exp(-r_2[i])))
    return(lrho)     

x = np.linspace(0.1, 10.1, 21)
y = np.linspace(0.1, 10.1, 21)
z = np.zeros(len(x))
'''
for i in range(10):
    plt.plot(x, rho(r_1(x,y,z+i/10), r_2(x,y,z+i/10)))
'''
# Gegenbauer polynomials

# Legendre polynomials:
'''
def legendrepoly(n, x):
    if n == 0:
        return([1 for _ in range(len(x))])
    elif n == 1:
        return(x)
    else:
        return((2*n + 1)/(n + 1)*x*legendrepoly(n-1,x) - n/(n+1)*legendrepoly(n-2,x))
#Px[n+1] = (2*n + 1)/(n + 1)*x*P[n] - n/(n+1)*P[n-1]


print(legendrepoly(0,x))
print(legendrepoly(1,x))
print(legendrepoly(2,x))

for i in range(5):
    plt.plot(x, legendrepoly(i,x))

'''

plt.show()