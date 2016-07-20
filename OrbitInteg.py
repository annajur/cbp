import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from astropy import units
from galpy.orbit import Orbit
from FerrersPotential import FerrersPotential as FP

def rot(omega, t):
    temp = [[np.cos(t*omega), -np.sin(t*omega)], [np.sin(t*omega), np.cos(t*omega)]]
    return np.array(temp)

def inrotframe(orbit, ts, potential):
    x, y = [], []
    for item in ts:
        x.append(orbit.x(item))
        y.append(orbit.y(item))
    xy = np.zeros([len(x),2])
    xy[:,0] = x
    xy[:,1] = y
    omega = potential.OmegaP()
    xrot, yrot = np.zeros(len(ts)), np.zeros(len(ts))
    for i in range(len(ts)):
        xrot[i],yrot[i] = np.dot(xy[i],rot(omega, ts[i]))
    return xrot, yrot

pmw = FP(a = 8*units.kpc, b = 0.35, c = 0.2375, normalize = True, omegab = 10.*units.km/units.s/units.kpc)

o = Orbit(vxvv = [1.,0.25,0.6,0.,0,0])
time = np.linspace(0,400,2000)
o.integrate(time, pmw, method = 'leapfrog')

o.plot(d1 = 'x', d2 = 'y')
plt.show()