import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from astropy import units
from galpy.orbit import Orbit
from FerrersPotential import FerrersPotential as FP
matplotlib.rcParams['figure.figsize'] = (10, 10)

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


initc = []
R = 1
svxvv = '[0.2,0.2,vT,0.,0.,0.]'

output_file = open('/home/annaj/cbp_usrp/pretty_pictures/in_rot/6/setvT/log','w')
output_file.write(svxvv+'\n'+'R = ')

'''
for i in range(1,41):
    R = 0.025*i
    initc.append([R,0.01,0.01,0.,0.,np.pi/2])
    output_file.write(str(R)+', ')
ts = np.linspace(0,100,2000)  
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

for i in range(1,41):
    vT = -0.7+0.05*i
    initc.append([0.2,0.2,vT,0.,0.,0.])
    output_file.write(str(vT)+', ')
ts = np.linspace(0,50,1500)
output_file.close()

orbits = []
for i in range(len(initc)):
    orbits.append(Orbit(vxvv = initc[i]))

for i in range(len(orbits)):
    try:
        orbits[i].integrate(ts, pmw)
        #orbits[i].plot(d1='x',d2='y', color = 'dodgerblue', overplot = False)
        ts = np.linspace(0,50,1500)
        xr, yr = inrotframe(orbits[i],ts,pmw)
        plt.plot(xr,yr, c = 'crimson')
        plt.xlabel(r'$x/R_0$', fontsize = 17)
        plt.ylabel(r'$y/R_0$', fontsize = 17)
        xr, yr = [],[]
        if i<10:
            name = '/home/annaj/cbp_usrp/pretty_pictures/in_rot/6/setvT/r0'+str(i)+'.png'
        else:
            name = '/home/annaj/cbp_usrp/pretty_pictures/in_rot/6/setvT/r'+str(i)+'.png'
        plt.savefig(name)
        plt.close()
    except (RuntimeWarning,FutureWarning,RuntimeError,AttributeError,DeprecationWarning,ZeroDivisionError):
        print('I am sorry, but something went wrong :(')

'''
--log--

out_rot:

[R,-0.1,0.1,0.,0.,0.]       R = 1 + 0.05*i
[R,-0.1,0.1,0.,0.,np.pi/4]  R = 1 + 0.05*i
[R,-0.1,0.1,0.,0.,np.pi/2]  R = 1 + 0.05*i

[R,-0.1,0.2,0.,0.,0.]
[R,-0.1,0.3,0.,0.,0.]

[R,-0.2,0.3,0.,0.,0.]

[1.2,vR,0.2,0.,0.,0.]       vR = -1.+ 0.03*i
[1.2,vR,0.1,0.,0.,0.]       vR = -1.+ 0.03*i
[1.35,vR,0.1,0.,0.,np.pi/4] vR = -.2 + 0.02*i

[1.2,0.1,0.5,0.,vz,0.001]   vz = 0.03*i
[1.2,0.1,0.1,0.,vz,0.001]   vz = 0.03*i

[1.2,0.1,0.5,0.,0.,phi]     phi = np.pi/20*i
[1.2,0.1,0.1,0.,0.,phi]     phi = np.pi/200*i

in_rot:

[R,0.,0.009,0.,0.,0.001]    R = 0.1*i
[R,0.1,0.009,0.,0.,0.0]     R = 0.05*i
[R,.1,.009,0.,0.,np.pi/4]   R = 0.05*i
[R,0.01,0.01,0.,0.,0.]      R = 0.05*i
[R,0.01,0.01,0.,0.,np.pi/2] R = 0.05*i
[R,0.01,0.01,0.,0.,np.pi/4] R = 0.05*i

[0.7,vR,0.1,0.,0.,0.]       vR = -0.8+0.1*i
[0.5,vR,0.1,0.,0.,0.]       vR = -0.2+0.05*i
[0.5,vR,0.1,0.,0.,np.pi/4]  vR = -0.2+0.05*i

[0.2,0.2,vT,0.,0.,0.]       vT = -0.7+0.05*i
[0.3,0.2,vT,0.,0.,0.]       vT = -0.7+0.05*i
[0.4,0.2,vT,0.,0.,0.]       vT = -0.7+0.05*i
[0.5,0.2,vT,0.,0.,0.]       vT = -0.7+0.05*i
[0.6,0.2,vT,0.,0.,0.]       vT = -0.7+0.05*i
[0.7,0.2,vT,0.,0.,0.]       vT = -0.7+0.05*i
[0.8,0.2,vT,0.,0.,0.]       vT = -0.7+0.05*i

[0.7,0.2,0.6,0.,vz,0.001]   vz = 0.03*i
[0.5,0.1,0.6,0.,vz,0.]      vz = 0.03*i
[0.5,0.1,0.6,0.,vz,np.pi/4] vz = 0.03*i

[0.7,0.2,0.6,0.,0.,phi]     phi = np.pi/20*i
[0.5,0.2,0.4,0.,0.,phi]     phi = np.pi/20*i
[0.5,0.3,0.3,0.,0.,phi]     phi = np.pi/20*i, 
'''