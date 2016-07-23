#!/u/juranova/anaconda2/bin/python
import sys
import os

import numpy as np

#input_filename = sys.argv[1]
#input_file = open(input_filename, 'r')
#dataline = input_file.read().split('\n')

R = []
#icon = []

for i in range(1,21):
    R.append(0.05*i)

counter = 0

for i in range(len(R)):
    for j in range(5): 
        vt, vz = 1, 1
        while (vt*vt+vz*vz)>1:
            vt = np.random.uniform(0.,1.)
            vz = np.random.uniform(0.,1.)
        #icon.append([round(R[i],2),0.,round(vt,3),0.,round(vz,3),0.])
        filename = 'icon_'+str(counter)+'.txt'
        outputfile = open(filename,'w')
        outputfile.write(str(round(R[i],2))+'\t0.\t'+str(vt)+'\t0.\t'+str(vz)+'\t0.' )
        outputfile.close()
        counter += 1