'''
Modification of the Hubra radar used in the Fimbulisen project.

Created on 20.10.2010

@author: Karsten Mueller
'''
# Built-in modules

# Third-party modules
import numpy as np

# Local modules
from processgpr.ext.hubra import HUBRA


class FIMBUL(HUBRA):
    '''Hubra radar used in Fimbulisen project'''

    def parameters(self):
        self.polarimetry = False  #: if True XX,XY,YX, and YY are loaded else only self.PP
        self.numfreq = 13000  #: number of samples in frequency domain
        self.gatefreq = 400e3  #: gating frequency
        self.Nfft = 2 ** 15  #:
        self.samples = None  #: number of samples in time domain, assigned after conversion to time domain
        self.traces = 0  #: number of traces in the profile
        self.timens = 22000.0  #: time window in ns
        self.lowF = 200000000  #: lowest system frequency
        self.highF = 500000000  #: highest system frequency
        self.fc = (self.highF + self.lowF) / 2.0  #: center frequency
        self.BW = self.highF - self.lowF  #: system bandwidth
        self.deltaF = self.BW / (self.numfreq - 1)  #: sampling step in the freq. domain
        self.nstart = int(self.lowF / self.deltaF)
        self.deltaT = 1.0 / (self.deltaF * self.Nfft)  #: sampling in the time domain in [ns]
        self.lam = self.cmed * 1.0e9 / self.fc  #: wave length in medium
        if self.samples:
            self.t = np.arange(0, self.deltaT * self.samples, self.deltaT) * 1e9  #: time axis in ns
            self.r = self.t * self.cmed / 2.0  #: range axis in m