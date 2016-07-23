import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import matplotlib

matplotlib.rcParams['figure.figsize'] = (10.0, 5.5)

''' RAW DATA PLOT '''

time = []
les = []

datafile = sys.argv[1]

readfile = open(datafile, 'r')
sepfile = readfile.read().split('\n')
readfile.close()
for i in range(5,len(sepfile)):
    if sepfile[i] != '':
        xandy = sepfile[i].split('\t')
        time.append(float(xandy[0]))
        les.append(float(xandy[1]))
        if time[i-5]<10:
                print(sepfile[i-5])

plt.loglog(time, les)
plt.xlabel(r'$t$', fontsize = 15)
plt.ylabel(r'$LE$', fontsize = 15)
plt.show()
