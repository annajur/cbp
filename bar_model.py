######################################################################
# Numpy imported as np!!
# Even I have imported the division from __future__, I still use dots to
#   force the result in float data type
# Units still need to be converted to natural units
#
######################################################################
import numpy as np
from __future__ import division

from scipy.special import gamma
from math import factorial
from scipy.special import gegenbauer
from scipy.special import lpmv # http://docs.scipy.org/doc/scipy-0.17.0/reference/generated/scipy.special.lpmv.html
from scipy.special import eval_gegenbauer

from galpy.potential_src.Potential import Potential, kms_to_kpcGyrDecorator, \
    _APY_LOADED
if _APY_LOADED:
    from astropy import units

def diracdelta(m,n):
    if m == n:
        return 1
    else:
        return 0

class MWBarPotential(Potential):
    """Class that implements bar potential presented in http://adsabs.harvard.edu/abs/2012MNRAS.427.1429W
    .. math::
        
    """
    
    def __init__(self, amp=1.,.,ro = None, vo = None): #####rozdelane, dokumentace okopirovana
       
        """
        NAME:
           __init__
        PURPOSE:
           initialize Milky Way potential presented in http://adsabs.harvard.edu/abs/2012MNRAS.427.1429W
        INPUT: ######## opravit!!!!!!!!!!!!!!!!!!!!!!!
           amp - amplitude to be applied to the potential (default: 1); can be a Quantity with units of mass or Gxmass
           a - scale length (can be Quantity)
           b - scale height (can be Quantity)
           normalize - if True, normalize such that vc(1.,0.)=1., or, if given as a number, such that the force is this fraction of the force necessary to make vc(1.,0.)=1.
           ro=, vo= distance and velocity scales for translation into internal units (default from configuration file)
        OUTPUT:
           (none)
        HISTORY:
           2016-06-20 - Started - Bovy (NYU)
        """
        
        self.hasC= False
        self.hasC_dxdv= False    
        self._nemo_accname= 'MWBar'

    ##########
    def transform(self, r = False, theta = False, rho = False, z = False): # (r, theta, phi) <-> (rho, phi, z)
        if not theta and not r:
            return(np.sqrt(rho**2 + z**2), np.arccos(z/np.sqrt(rho**2 + z**2))) # (rho, z) transformed to (r, theta)
        elif not z and not rho:
            return(r*np.sin(theta), r*np.cos(theta)) # (r, theta) transformed to (rho, z)
        else:
            return('please specify coordinates you wish to transform') # should be unnecessary later on
    ##########

    Anlm = np.array([
            [ 0,  0,  0,  1.509],
            [ 1,  0,  0, -0.086],
            [ 2,  0,  0, -0.033],
            [ 3,  0,  0, -0.020],
            [ 0,  2,  0, -2.606],
            [ 1,  2,  0, -0.221],
            [ 2,  2,  0, -0.001],
            [ 0,  2,  2,  0.665],
            [ 1,  2,  2,  0.129],
            [ 2,  2,  2,  0.006],
            [ 0,  4,  0,  6.406],
            [ 1,  4,  0,  1.295],
            [ 0,  4,  2, -0.660],
            [ 1,  4,  2, -0.140],
            [ 0,  4,  4,  0.044],
            [ 1,  4,  4, -0.012],
            [ 0,  6,  0, -5.859],
            [ 0,  6,  2,  0.984],
            [ 0,  6,  4, -0.030],
            [ 0,  6,  6,  0.001],])

    ##########
    def Knl(n,l):
        return(0.5*n*(n + 4*l + 3) + (l + 1)*(2*l + 1))

    def Inl(n,l,m):
        if m<=l:
            return(Knl(n,l)*1/(2**(8*l + 6))*gamma(n + 4*l + 3)/(factorial(n)*(n + 2*l + 3/2) *\
                (gamma(2*l + 3/2))**2 )*(1 + diracdelta(m,0))*np.pi*2/(2*l + 1)*(factorial(l + m)/factorial(l - m)))
        else:
            return('m must be smaller or equal to l') # just a temporary solution for bad input

    def Phinlm(n,l,m,r,phi,z):
        r, theta = transform(rho = r, z = z)
        s = r/1. # 1. is in kpc!!!!!
        return(s**l/((1 + s)**(2*l + 1))*eval_gegenbauer(n, 2*l + 3/2, (s - 1)/(s + 1))*np.cos(m*phi)*lpmv(m, l, np.cos(theta)))
    print(Phinlm(1,3,1,3,1,1))
    ##########

    def _evaluate(self,R,phi,z,t=0.):
        """
        NAME:
            _evaluate
        PURPOSE:
            evaluate the potential at R, theta, phi
        INPUT:
            R - Galactcentric spherical radius
            phi - azimuth
            z - vertical height
            t - time
        OUTPUT:
            Phi(R,theta,phi)
        HISTORY:
            2016-06-20 - Started - Juranova (MU)
        """
        # the constant (GM_{bar}/r_s) is set to one... but is that correct...? Now I don't think so...
        Phi = 0
        for i in range(len(Anlm)):
            Phi += Anlm[i,3]*Phinlm(Anlm[i,0], Anlm[i,1], Anlm[i,2], R, phi, z)
        return(Phi)

    def _Rforce(self,R,theta,phi,t=0.):
        """
        NAME:
        _Rforce
        PURPOSE:
        evaluate the radial force for this potential
        INPUT:
        R - Galactocentric cylindrical radius
        phi - azimuth
        z - vertical height
        t - time
        OUTPUT:
        the radial force
        HISTORY:
        2016-06-20 - Started - Juranova (MU)
        """
        return !!!!!!!!!!!!!
        
    def _dens(self, R, phi, z, t=0., x_0 = 1.49, y_0 = 0.58, z_0 = 0.40, q = 0.6, rho_0 = 1): # units!!!!!
        """
        NAME:
            _dens
        PURPOSE:
            evaluate the density for this potential
        INPUT:
            R - Galactocentric cylindrical radius
            z - vertical height
            phi - azimuth
            t - time
            x_0 - ##########
            y_0 - ##########
            z_0 - ##########
            q - #########
            rho_0 - ##########
        OUTPUT:
            the density
        HISTORY:
        2016-06-20 - Written - Juranova (MU)
        """
        x = R*np.cos(phi)
        y = R*np.sin(phi)
        r_1 = (((x/x_0)**2 + (y/y_0)**2)**2 + (z/z_0)**4)**0.25
        r_2 = ((q**2*(x**2 + y**2) + z**2)/z_0**2)**0.5
        return(rho_0*(np.exp(-r_1[i]**2/2.) + r_2[i]**(-1.85)*np.exp(-r_2[i])))
        
    def _R2deriv(self,R,z,phi=0.,t=0.):
        """
        NAME:
            _R2deriv
        PURPOSE:
            evaluate the second radial derivative for this potential
        INPUT:
            R - Galactocentric cylindrical radius
            z - vertical height
            phi - azimuth
            t - time
        OUTPUT:
            the second radial derivative
        HISTORY:
            2011-10-09 - Written - Bovy (IAS)
        """