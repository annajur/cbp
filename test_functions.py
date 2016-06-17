import numpy as np
from math import exp

import matplotlib.pyplot as plt

# Gegenbauer polynomials

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
    #return(rho_0*(exp(-r_1**2/2.) + r_2**(-1.85)*exp(-r_2)))

x = np.linspace(0.1, 10., 100)
y = np.linspace(0.1, 10., 100)
z = np.zeros(len(x))

r1 = r_1(x,y,z)
r2 = r_2(x,y,z)

plt.plot(rho(r1, r2), x)
plt.show()