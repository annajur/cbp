# cython: boundscheck=False
# cython: nonecheck=False
# cython: cdivision=True
# cython: wraparound=False

import numpy as np
cimport numpy as np

from scipy import integrate
from scipy.optimize import fsolve

cdef extern from "math.h":
    double sqrt(double)

cdef double cy_FracInt(double x, double y, double z, double a2, double b2, double c2, double tau, double n):
    return (1 - x*x/(a2 + tau) - y*y/(a2*b2 + tau) - z*z/(a2*c2 + tau))**n/sqrt((a2 + tau)*(a2*b2 + tau)*(a2*c2 + tau))

cdef double cy_integrand_x(double tau, double x, double y, double z, double a2, double b2, double c2, double n):
    return -x/(a2 + tau) * cy_FracInt(x,y,z,a2,b2,c2,tau,n)

cdef double cy_integrand_y(double tau, double x, double y, double z, double a2, double b2, double c2, double n):
    return -y/(a2*b2 + tau) * cy_FracInt(x,y,z,a2,b2,c2,tau,n)

cdef double cy_integrand_z(double tau, double x, double y, double z, double a2, double b2, double c2, double n):
    return -z/(a2*c2 + tau) * cy_FracInt(x,y,z,a2,b2,c2,tau,n)

cpdef double func(double tau, double x, double y, double z, double a2, double b2, double c2):
    return x*x/(a2+tau) + y*y/(a2*b2+tau) + z*z/(a2*c2+tau) - 1.

cpdef double lowerlim(double x, double y, double z, double a2, double b2, double c2):
    """Lower limit of the integrals"""
    if sqrt(x*x/a2 + y*y/(a2*b2) + z*z/(a2*c2)) >= 1:
        return fsolve(func,0,args=(x,y,z,a2,b2,c2))[0]
    else:
        return 0

cpdef cy_forceInt(double x, double y, double z, double a2, double b2, double c2, int i, double n):
    if i == 0:
        return integrate.quad(cy_integrand_x, lowerlim(x,y,z,a2,b2,c2), np.inf, args=(x,y,z,a2,b2,c2,n))[0]
    elif i == 1:
        return integrate.quad(cy_integrand_y, lowerlim(x,y,z,a2,b2,c2), np.inf, args=(x,y,z,a2,b2,c2,n))[0]
    elif i == 2:
        return integrate.quad(cy_integrand_z, lowerlim(x,y,z,a2,b2,c2), np.inf, args=(x,y,z,a2,b2,c2,n))[0]
