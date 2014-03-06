'''
Created on 26. nov. 2010

@author: karstm
'''
import os

from processgpr.ext.hubra import HUBRA


wdir = r'/mn/terra/hubradata/Hubra_Svalbard/Kongsvegen_2008'
file = 'data Thu Apr 24 14-48-55 2008 gated  sweep 1 MHz0'
pro = HUBRA()

pro.convertHDF_prestack(os.path.join(wdir, file), wdir,
                        1200, combine=True, stack=10, filerange="all-all",
                        envelope=False)

os.rename(os.path.join(wdir, file + '.h5'),
          os.path.join(wdir, "KV08_stake67_pha_stk10.h5"))
#print "Finished conversion for file: %s" % os.path.abspath(file)
#print "Total Distance: %f (%i Traces)" % (pro.coords.total_distance, pro.traces)
