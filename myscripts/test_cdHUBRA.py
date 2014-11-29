'''
Created on 29. nov. 2014

@author: kmu
'''
import os

from processgpr.ext.hubra import cdHUBRA


wdir = r'C:\Users\kmu\OneDrive\Dev\PYworkspace\chda\radar_test_data'
file = 'data5610'
pro = cdHUBRA()
pro.filename = os.path.join(wdir, file)
print "numfreq=", pro.numfreq
print "filename=", pro.filename
pro.readRaw()

#
# pro.convertHDF_prestack(os.path.join(wdir, file), wdir,
#                         1000, combine=True, stack=20, filerange="all-all",
#                         envelope=False)
#
# print "Finished conversion for file: %s" % os.path.abspath(file)
# try:
#     print "Total Distance: %f (%i Traces)" % (pro.coords.total_distance, pro.traces)
# except AttributeError:
#     pass
#
# os.rename(os.path.join(wdir, file + '.h5'),
#           os.path.join(wdir, 'AT2009-01-03_1000ns_stk20.h5'))
