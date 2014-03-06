'''
Converts raw C-band data to time domain and stores it as '.h5' files.

Created on 20.10.2010

@author: Karsten Mueller
'''
import os

from processgpr.ext.cband import CBAND

procdir = r'Y:\Hubra_Svalbard\Kongsvegen_2006'
file = 'C_KV06_stake67_env'
os.chdir(procdir)
fout = os.path.join(procdir, file)
pro = CBAND()
pro.numfreq = 301
pro.convertHDF(fout, combine=False, filerange="31-31", envelope=True, stack=1,
               timens=1000)
