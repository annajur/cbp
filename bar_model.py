######################################################################
# vzorec pro phi 
# pozor na import numpy!!
# i pres import division pouzivam notaci s teckami pro normalni deleni
#
#
#
######################################################################
import numpy as nu
from __future__ import division
from galpy.potential_src.Potential import Potential, kms_to_kpcGyrDecorator, \
    _APY_LOADED
if _APY_LOADED:
    from astropy import units

class MWBarPotential(Potential):
    '''Class that implements bar potential presented in http://adsabs.harvard.edu/abs/2012MNRAS.427.1429W
    .. math::
        
    '''
    
    
def __init__(self, amp=1.,.,)