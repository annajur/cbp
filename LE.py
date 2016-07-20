import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from astropy import units
from galpy.orbit import Orbit
from FerrersPotential import FerrersPotential as FP

def allorbits(x,y):
    xo = [x[2*i] for i in range(int(len(x)/2))]
    xo1 = [x[2*i+1] for i in range(int(len(x)/2))]
    yo = [y[2*i] for i in range(int(len(x)/2))]
    yo1 = [y[2*i+1] for i in range(int(len(x)/2))]
    return [xo,yo],[xo1,yo1]

x = []
y = []
def evolveorbit(icon, ti, tau, pot):
    global x
    global y
    o = Orbit(vxvv=icon) # [R,vR,vT,z,vz,phi]
    tf = ti+tau
    ts = np.linspace(ti,tf,100)
    o.integrate(ts, pot, method = 'leapfrog')
    x.append(o.x(ts[0]))
    y.append(o.y(ts[0]))
    return [o.R(tf),o.vR(tf),o.vT(tf),o.z(tf),o.vz(tf),o.phi(tf)]

def dvector(o,d):
    return np.array(d)-np.array(o)

def alpha(delta):
    size = np.linalg.norm(delta)
    return size, delta/size

def LEs(time, size, initsize):
    return np.log(size/initsize)

def lyapunov(o,tau, potential, Tm):
    global x
    global y
    x,y = [],[]
    time,LE = [],[]
    continuing = True
    i = 1
    w = [1.,0.,0.,0.,0.,0.]
    initsize = 1e-5
    while continuing:
        icdo = list(np.array(o)+initsize*np.array(w))
        newo = evolveorbit(o, tau*i, tau, potential)
        newdo = evolveorbit(icdo, tau*i, tau, potential)
        wj = dvector(o=newo,d=newdo)
        size, veps0 = alpha(wj)
        LE.append(LEs(tau*i, size, initsize))
        time.append(tau*i)
        if i*tau > Tm:
            break
        i += 1
        o = newo 
        w = veps0
    A = np.array([sum(LE[:i]) / time[i] for i in range(len(time))])
    return A, time