import matplotlib.pyplot as plt
from astropy import units
from galpy.potential import plotPotentials

from bar_model import MWBarPotential as MWBP
from template import MiyamotoNagaiPotential as MNP

mwbp = MWBP()
mnp = MNP()

plotPotentials(mwbp)
