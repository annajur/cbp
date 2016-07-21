#!usr/bin/python

import sys
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from astropy import units
from galpy.orbit import Orbit
from FerrersPotential import FerrersPotential as FP
matplotlib.rcParams['figure.figsize'] = (10, 10)


if len(sys.argv) != 3:
    print ('Usage: %s DATA, specify both input file containing initial conditions and log file' % (os.path.basename(sys.argv[0])))
    sys.exit(1)
input_filename = sys.argv[1]
input_file = open(input_filename, 'r')
dataline = input_file.read().split('\n')

for line in dataline:
    if line != '':
        initc_string = line.split('\t')
input_file.close()

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

output_filename = sys.argv[2]
output_file = open(output_filename, 'w')
for i in range(len(initc_string)):
    output_file.write(initc_string[i]+', ')

phi = [i*np.pi/12 for i in range(12)]

initc = [float(initc_string[i]) for i in range(len(initc_string))]

print(initc)
print('baf')

vxvvs = []
# ====== setting vR ======
for i in range(4):
    vR = -0.7+0.05*i
    output_file.write(str(vR)+', ')
    tmp = initc[:]
    tmp[1] = vR
    vxvvs.append(tmp)
output_file.close()

pmw = FP(a = 8*units.kpc, b = 0.35, c = 0.2375, normalize = True, omegab = 10.*units.km/units.s/units.kpc)
ts = np.linspace(0,50,1500)

for i in range(len(phi)):
    for j in range(len(vxvvs)):
        vxvvs[j][5] = phi[i]
        orbit = Orbit(vxvv = vxvvs[j])
        try:
            orbit.integrate(ts, pmw)
            #orbit.plot(d1='x',d2='y', color = 'dodgerblue', overplot = False)
            ts = np.linspace(0,50,1000)
            xr, yr = inrotframe(orbit,ts,pmw)
            plt.plot(xr,yr, c = 'crimson')
            plt.xlabel(r'$x/R_0$', fontsize = 17)
            plt.ylabel(r'$y/R_0$', fontsize = 17)
            xr, yr = [],[]
            if i<10:
                name = str(i)+'pi_12_0'+str(j)+'.png'
            else:
                name = str(i)+'pi_12_'+str(j)+'.png'
            plt.savefig(name)
            plt.close()
        except (RuntimeWarning,FutureWarning,RuntimeError,AttributeError,DeprecationWarning,ZeroDivisionError):
            print('I am sorry, but something went wrong :(')

# create submit file, compile script, check where are files saved and if figures are saved as stdout, 
# and what I can do about output directory