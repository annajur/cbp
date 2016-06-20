######################################################################
# Numpy imported as np!!
# Even I have imported the division from __future__, I still use dots to
#   force the result in float data type 
#
#
######################################################################
import numpy as np
from __future__ import division
from galpy.potential_src.Potential import Potential, kms_to_kpcGyrDecorator, \
    _APY_LOADED
if _APY_LOADED:
    from astropy import units

class MWBarPotential(Potential):
    '''Class that implements bar potential presented in http://adsabs.harvard.edu/abs/2012MNRAS.427.1429W
    .. math::
        
    '''
    
    
def __init__(self, amp=1.,.,):


def _evaluate(self,R,theta,phi,t=0.):
    """
    NAME:
        _evaluate
    PURPOSE:
        evaluate the potential at R, theta, phi
    INPUT:
        R - Galactcentric spherical radius
        theta - 
        phi - azimuth
        t - time
    OUTPUT:
        Phi(R,theta,phi)
    HISTORY:
        2016-06-20 - Started - Juranova (MU)
    """
    return !!!!!!!!!!!!!

def _Rforce(self,R,theta,phi,t=0.):
    """
    NAME:
       _Rforce
    PURPOSE:
       evaluate the radial force for this potential
    INPUT:
       R - Galactocentric cylindrical radius
       theta - 
       phi - azimuth
       t - time
    OUTPUT:
       the radial force
    HISTORY:
       2016-06-20 - Started - Juranova (MU)
    """
    return !!!!!!!!!!!!!
    
def _dens(self,R,z,phi=0.,t=0.):
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
    OUTPUT:
        the density
    HISTORY:
       2016-06-20 - Started - Juranova (MU)
    """
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