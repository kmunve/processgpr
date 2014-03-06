'''
Created on 16. nov. 2010

@author: Karsten Mueller
'''

import os

import matplotlib.pyplot as plt

from processgpr.ext.fimbul import FIMBUL


wdir = r'C:\hubradata\Fimbul091230'
os.chdir(wdir)
file = 'data1.h5'
pro = FIMBUL()
pro.readH5(file)
print "Total Distance: ", pro.coords.total_distance, pro.traces
pro.distanceStack(10.0)
print "Total Distance: ", pro.coords.total_distance, pro.traces
pro.saveH5(os.path.splitext(file)[0] + "_ds10")
plt.scatter(pro.coords.x, pro.coords.y, c=pro.coords.elev, edgecolors='none')
plt.show()