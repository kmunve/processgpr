'''
Created on 20.10.2010

@author: Karsten Mueller
'''
import os

import matplotlib.pyplot as plt

from processgpr.ext.fimbul import FIMBUL


wdir = r'C:\hubradata\Fimbul091230'
file = 'data1'
pro = FIMBUL()

pro.convertHDF(os.path.join(wdir, file), wdir,
               6000, combine=True, stack=3, filerange="all-all",
               envelope=True, gpsfile=os.path.join(wdir, file + '.cor'),
               projection='PS71')

print "Total Distance: ", pro.coords.total_distance, pro.traces
plt.scatter(pro.coords.x, pro.coords.y, c=pro.coords.elev, edgecolors='none')
plt.show()