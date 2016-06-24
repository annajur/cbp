import matplotlib.pyplot as plt
from astropy import units
from galpy.potential import plotPotentials

from FerrersPotential import FerrersPotential as FP
from template import MiyamotoNagaiPotential as MNP

fp = FP()
mnp = MNP()

plotPotentials(fp)
