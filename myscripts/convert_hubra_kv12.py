'''
Created on 26. nov. 2010

@author: karstm
'''
import os

from processgpr.ext.hubra import HUBRA


#wdir = r'/mn/terra/hubradata/Hubra_Svalbard/Kongsvegen_2008'
wdir = r'Y:\Hubra_Svalbard\Kongsvegen_2008'
file = 'data Thu Apr 24 13-15-42 2008 gated  sweep 1 MHz0'
pro = HUBRA()

pro.convertHDF(os.path.join(wdir, file), wdir,
               1200, combine=True, stack=1, filerange="all-all",
               envelope=False, gpsfile=os.path.join(wdir, file + '.cor'),
               projection='UTM', utmzone=33)

print "Finished conversion for file: %s" % os.path.abspath(file)
print "Total Distance: %f (%i Traces)" % (pro.coords.total_distance, pro.traces)
