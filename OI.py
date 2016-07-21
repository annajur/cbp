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

#print(input_filename)

input_file = open(input_filename, 'r')
dataline = input_file.read().split('\n')
#print(dataline)
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

pmw = FP(a = 8*units.kpc, b = 0.35, c = 0.2375, normalize = True, omegab = 10.*units.km/units.s/units.kpc)

output_filename = sys.argv[2]
output_file = open(output_filename, 'w')

phi = [i*np.pi/12 for i in range(12)]

initc = [float(initc_string[i]) for i in range(len(initc_string))]

print('baf')

# ====== setting vR ======
for i in range(1,3):
    initc[1] = -0.7+0.05*i
    vxvvs.append(initc)

# ====== setting vT ======
#for i in range(1,41):
#    initc[2] = -0.7+0.05*i
#    vxvvs.append(initc)

print(vxvvs)

for i in range(len(phi)):
    initc.append(phi[i])
    for j in range(1,41):
        try:
            orbits[i].integrate(ts, pmw)
            #orbits[i].plot(d1='x',d2='y', color = 'dodgerblue', overplot = False)
            ts = np.linspace(0,50,1000)
            xr, yr = inrotframe(vxvvs[i],ts,pmw)
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


'''
R = 1
svxvv = '[R,0.01,0.01,0.,0.,np.pi/2]'

output_file = open('/home/annaj/cbp_usrp/pretty_pictures/in_rot/4/setR/log','w')
output_file.write(svxvv+'\n'+'R = ')

for i in range(1,41):
    R = 0.025*i
    initc.append([R,0.01,0.01,0.,0.,0.]) #np.pi/2])
    output_file.write(str(R)+', ')
ts = np.linspace(0,100,2400)  
output_file.close()


'''
'''
for i in range(1,21):
    vR =  -0.2 + 0.02*i
    initc.append([1.35,vR,0.1,0.,0.,np.pi/4])
    output_file.write(str(vR)+', ')
ts = np.linspace(0,100,800)
output_file.close()
'''
'''
orbits = []
for i in range(len(initc)):
    orbits.append(Orbit(vxvv = initc[i]))


output_file.close()