'''
Created on 26. nov. 2010

@author: karstm
'''
import os

from processgpr.ext.hubra import HUBRA


wdir = r'/mn/terra/hubradata/Data2008_09/Wisdom/0301'
file = 'data2'
pro = HUBRA()

pro.convertHDF_prestack(os.path.join(wdir, file), wdir,
                        1000, combine=True, stack=20, filerange="all-all",
                        envelope=False)

print "Finished conversion for file: %s" % os.path.abspath(file)
try:
    print "Total Distance: %f (%i Traces)" % (pro.coords.total_distance, pro.traces)
except AttributeError:
    pass

os.rename(os.path.join(wdir, file + '.h5'),
          os.path.join(wdir, 'AT2009-01-03_1000ns_stk20.h5'))