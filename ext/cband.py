__docformat__ = "reStructuredText"
'''
Extends the GPR class for FMCW C-band functionality.

:Created on 20.10.2010

:Author: Karsten Mueller
'''
import os
# Third-party modules
import numpy as np

from processgpr.core.gpr import GPR
from processgpr.core.gps import GPS


class CBAND(GPR):
    def parameters(self):
        self.polarimetry = False  #: if True XX,XY,YX, and YY are loaded else only XX
        self.subdirect = True
        self.PP = 'YY'
        self.top = 218  #: top defines the sample where the snow surface starts
        self.srange = 457
        self.numfreq = 451  #: number of frequencies
        self.Nfft = 2 ** 15  #: number of FFT samples
        self.samples = self.Nfft / 2  #: number of samples
        self.traces = 0  #: number of traces in the profile
        self.lowF = 4.8e9  #: lowest system frequency
        self.highF = 5.8e9  #: highest system frequency
        self.fc = 5.3e9  #: center frequency
        self.BW = self.highF - self.lowF  #: system bandwidth
        self.deltaF = self.BW / (self.numfreq - 1)  #: sampling step in the freq. domain
        self.nstart = int(self.lowF / self.deltaF)
        self.deltaF = self.BW / (self.numfreq - 1)  #: sampling step in the freq. domain
        self.deltaT = 1 / (self.deltaF * self.Nfft)  #: sampling in the time domain
        self.filename = None
        self.t = np.arange(0, self.deltaT * self.samples, self.deltaT) * 1e9  #: time axis in ns
        self.lam = self.cmed * 1.0e9 / self.fc  #: wave length in medium
        self.r = self.t * self.cmed / 2.0  #: range axis in m


    def readRaw(self, removeDW=True):
        """
        Reads the binary data from file 'filename'!
        Note! Enter filename without extension.
        """
        self.info('Reading data from %s ...' % self.filename)
        self.PP = self.filename[-2:]
        fid = open(self.filename, 'rb')
        self.data = np.fromstring(fid.read(), np.float32)
        fid.close()

        self.traces = len(self.data) / self.numfreq
        self.data.shape = (self.traces, self.numfreq)

        if removeDW:
            from scipy.io import loadmat

            data = loadmat('Y:\Hubra_Svalbard\Kongsvegen_2006\directwaveyy.mat')
            dwyy = data['directwaveyy']
            dwyy.shape = (1, dwyy.shape[0])
            self.data = self.data - dwyy

        self.done()


    def addCoords(self, gpsfile):
        self.info("Adding GPS coordinates from file: %s ..." % gpsfile)
        gps = GPS(projection='UTM', utmzone='33')

        if os.path.splitext(gpsfile)[1] == '.mat':
            gps.readMAT(gpsfile)
        else:
            gps.readGPS(gpsfile)
        gps.interpolateProj(self.traces * self.nstacks)
        if self.nstacks > 1:
            gps.stackCoords(self.traces)
        self.coords = gps
        self.hascoords = True
        self.done()


    def convertHDF(self, fout, timens=200, stack=1, combine=False, filerange="all-all",
                   gpsfile=None, envelope=True, ww=0, cf=0):
        from processgpr.core.h5file import H5FILE

        h5file = H5FILE()
        startcount = filerange.split("-")[0]
        stopcount = filerange.split("-")[1]
        if stopcount == "all":
            stopcount = "99999999999"

        if startcount == "all":
            filecounter = 0
            self.filename = "profile" + str(filecounter) + ".gYY"
            gpsfile = "profile" + str(filecounter) + ".gps"
        else:
            filecounter = int(startcount)
            self.filename = "profile" + startcount + ".gYY"
            if gpsfile == None:
                gpsfile = "profile" + startcount + ".gps"

        self.readRaw()
        self.applyWindow(ww=ww, cf=cf)
        self.zeropad()
        self.applyFFT(envelope=envelope)
        self.addCoords(gpsfile)
        #        self.staticCorrection(1070)
        if stack > 1:
            self.stackTraces(stack)

        self.savesamples = int((self.samples / self.t[-1]) * timens)
        if self.savesamples > self.samples:
            self.savesamples = self.samples
        h5file.new(fout, (self.savesamples, 0))
        h5file.set_metadata(self)
        h5file.append_earray_data(self.data[0:self.savesamples, :])
        h5file.append_earray_gps(self.coords)
        if not combine:
            h5file.close()

        filecounter += 1

        while (os.path.exists("profile" + str(filecounter) + ".gYY") and filecounter <= int(stopcount)):
            print "...continuing"
            self.filename = "profile" + str(filecounter) + ".gYY"
            gpsfile = "profile" + str(filecounter) + ".gps"
            print "Converting file: ", self.filename
            self.readRaw()
            self.applyWindow(ww=ww, cf=cf)
            self.zeropad()
            self.applyFFT(envelope=envelope)
            self.addCoords(gpsfile)
            if stack > 1:
                self.stackTraces(stack)
            #            if not combine:
            #                fout = os.path.join(dirout, fileinname+"continued"+str(filecounter))
            #                self.savesamples = int((self.samples/self.t[-1])*timens)
            #                h5file.new(fout, (self.savesamples, 0))
            #                h5file.set_metadata(self)
            h5file.append_earray_data(self.data[0:self.savesamples, :])
            h5file.append_earray_gps(self.coords)
            if not combine:
                h5file.close()
            filecounter += 1
        if combine:
            h5file.close()

        from fnmatch import fnmatch

        for file in os.listdir('.'):
            if fnmatch(file, 'tmp*.dat'):
                os.remove(file)


    def wscs(self):
        '''
        Calculating the weighted scattering cross section WSCS.
        Input is amplitude data from which the power is calculated P=A**2
        '''
        self.info('Calculating weighted RCS ...')
        Pr = self.data
        G = 10 ** (15 / 10)  # antenna gain in dB
        Pt = 1.0  # transmitted power in W
        beamangle = 61 * np.pi / 180  # in rad
        dr = self.deltaT * self.cmed * (1.0e9 / 2.0)  # depth increment
        Vs = np.pi * (2 * self.r * np.tan(beamangle * 0.5)) ** 2 * dr  # scattering volume
        C = (4 * np.pi) ** 3 / (Pt * G ** 2 * self.lam ** 2)  # constant in radar eq.
        sigw = np.zeros((self.samples, self.traces))
        for i in range(self.traces):
            sigw[:, i] = ((Pr[:, i] * self.r ** 2) / np.sqrt(Vs)) * C

        self.data = sigw
        self.procext = '_scs'
        self.done()
    
    