'''
Defines functionality for Mala GPR systems.

Created on 20.10.2010

@author: Karsten Mueller
'''
import os

import numpy as np
import numpy.fft as fft

from processgpr.core.gpr import GPR


class MALA(GPR):
    ''' Impulse GPR data '''

    def parameters(self):
        pass


    def loadrd3(self, filename):
        headerfile = filename + '.rad'
        datafile = filename + '.rd3'
        self.filename = os.path.basename(filename)
        # Read the header file
        fH = open(headerfile, 'r')
        sH = fH.readlines()
        fH.close()

        self.samples = int(sH[0].split(':')[1].split('\\')[0])
        self.traces = int(sH[22].split(':')[1].split('\\')[0])
        self.timewindow = float(sH[18].split(':')[1].split('\\')[0])  # in ns
        self.deltaT = (self.timewindow * 1e-9) / (self.samples - 1)
        self.t = np.arange(self.samples) * self.deltaT * 1e9  # time in ns
        self.r = self.t * self.cmed / 2.0  #: range axis in m
        self.freq = fft.fftshift(fft.fftfreq(self.samples, self.deltaT))
        self.freqp = self.freq[self.samples / 2:]  # take only the positive
        print 'samples: ', self.samples
        print 'traces: ', self.traces
        print 'time window: ', self.timewindow
        print 'delta T: ', self.deltaT
        print 'freq: ', self.freq

        # Read the data file
        fD = open(datafile, 'rb')
        sD = fD.read()
        fD.close()
        data = np.fromstring(sD, 'int16')
        B = np.reshape(data, (self.traces, self.samples))
        #        Xaxis = np.arange(self.samples)
        #        Yaxis = np.arange(self.traces)
        self.data = np.transpose(B)


    def loadReflex(self, filename):
        ''' Reads the parameter and data files in Reflex format '''
        # Read parameter file
        pass


    def subtractDC(self):
        for i in range(self.traces):
            DCshift = np.mean(self.data[int(self.samples * 0.2):self.samples, i], axis=0)
            self.data[:, i] = self.data[:, i] - DCshift


    def rangeGain(self, a=1.0, b=2.0):
        import pylab

        g = a * self.r ** b
        pylab.plot(self.r, g)
        for i in range(self.traces):
            self.data[:, i] = self.data[:, i] * g


    def freqFilter(self, cut=1):
        Nfft = self.samples
        for i in range(self.traces):
            f = fft.fftshift(fft.fft(self.data[:, i], Nfft))
            for j in range(cut):
                f[j] = 0.0 + 0.0j
            self.data[:, i] = fft.ifft(f, Nfft)


    def bandpass(self, lowf, upf):
        Nfft = self.samples
        lows = np.where(np.abs(self.freq) < lowf * 1e6)
        ups = np.where(np.abs(self.freq) > upf * 1e6)
        for i in range(self.traces):
            f = fft.fftshift(fft.fft(self.data[:, i], Nfft))
            f[lows] = 0.0 + 0.0j
            f[ups] = 0.0 + 0.0j
            self.data[:, i] = fft.ifft(fft.ifftshift(f), Nfft)
