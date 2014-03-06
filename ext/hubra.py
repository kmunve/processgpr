'''
Specific methods for the I{HUBRA} radar.

The Hubra FMCW radar was developed at the Norwegian Defence Research Establishment.
It was applied on the 2nd stage of the Norwegian-U.S. Scientific Traverse across East Antarctica.

Created on 20.10.2010

@author: Karsten Mueller
'''
# Built-in modules
import os
import logging

# Third-party modules
import numpy as np

# Local modules
from processgpr.core.gpr import GPR
from processgpr.core.gps import GPS


class HUBRA(GPR):
    def parameters(self):
        self.polarimetry = False  #: if True XX,XY,YX, and YY are loaded else only self.PP
        self.numfreq = 10000  #: number of samples in frequency domain
        self.gatefreq = 1e6  #: gating frequency
        self.Nfft = 2 ** 15  #: number of samples for the FFT
        self.samples = None  #: number of samples in time domain, assigned after conversion to time domain
        self.traces = 0  #: number of traces in the profile
        self.timens = 2000.0  #: time window in ns
        self.lowF = 500000000  #: lowest system frequency
        self.highF = 3000000000  #: highest system frequency
        self.fc = (self.highF + self.lowF) / 2.0  #: center frequency
        self.BW = self.highF - self.lowF  #: system bandwidth
        self.deltaF = self.BW / (self.numfreq - 1)  #: sampling step in the freq. domain
        self.nstart = int(self.lowF / self.deltaF)
        self.deltaT = 1.0 / (self.deltaF * self.Nfft)  #: sampling in the time domain in [ns]
        self.lam = self.cmed * 1.0e9 / self.fc  #: wave length in medium
        if self.samples:
            self.t = np.arange(0, self.deltaT * self.samples, self.deltaT) * 1e9  #: time axis in ns
            self.r = self.t * self.cmed / 2.0  #: range axis in m

    def set_gatefreq(self, val):
        self.gatefreq = val


    def readRaw(self):
        '''
        Read HUBRA raw data.
        '''
        currentFileTraces = os.stat(self.filename)[6] / (
        self.numfreq * 2)  # the factor "2" is due to the 2-bytes of an int16
        print "%d traces" % (currentFileTraces)
        wholefile = np.fromfile(self.filename, 'int16', count=currentFileTraces * self.numfreq)
        self.info("Reading file %s ..." % self.filename)
        tmpdata = []
        for i in range(0, currentFileTraces):
            tmpdata.append(np.array(wholefile[i * self.numfreq:(i + 1) * self.numfreq]))
        self.data = np.fliplr(np.array(tmpdata, dtype=np.int16))
        self.traces = self.data.shape[0]
        #        self.OrigNumberOfTraces =
        self.done()


    #    def removeGating(self, view=False):
    #        self.info("Removing gating gain ...")
    #        N = int(self.Nfft/2)
    #        gategain = np.ones((N,1))
    #        gategain[:,0] = np.bartlett(N) # creates a triangle
    #        self.gategain = 1.0/gategain[self.samplendx]
    #
    #        if view==True:
    #            from pylab import show, plot, xlabel, ylabel
    #            x=np.linspace(0,2.0/self.gatefreq,N)*1e9
    #            print self.gatefreq, self.cmed
    #            plot(x[self.samplendx],10*np.log10(1.0/gategain[self.samplendx]))
    ##            plot(x[self.samplendx], self.gategain)
    #            xlabel("Time [ns]")
    #            ylabel("Gating Gain [dB]")
    #            show()
    #
    #        self.data = self.data*self.gategain
    #        self.done()

    def removeGating(self, view=False):
        self.info("Removing gating gain (gating freq.: %.2f Hz)..." % self.gatefreq)
        N = int(1.0 / (self.gatefreq * self.deltaT))  # samples for one gating gain cycle
        maxN = int(self.Nfft / 2)  # max number of samples
        gategain = np.ones((maxN, 1))  # vector with length=Nfft/2
        n = int(np.floor(maxN / N))  # how many gating cycles are included in time window
        rest = maxN - (n * N)  # how many samples are left after n full cycles?
        for i in range(n):
            gategain[i * N:(i + 1) * N, 0] = np.bartlett(N)  # creates one triangle for full cycle
        gategain[n * N:, 0] = np.bartlett(N)[:rest]  # creates part of triangle for incomplete cycle
        self.gategain = 1.0 / gategain[self.samplendx]

        if view == True:
            from pylab import show, plot, xlabel, ylabel

            plot(self.t[:len(self.samplendx)], 10 * np.log10(1.0 / gategain[self.samplendx]))
            xlabel("Time [ns]")
            ylabel("Gating Gain [dB]")
            show()

        self.data = self.data * self.gategain
        self.done()


    def compensateVivaldi(self, power=True):
        ''' System correction for the Vivaldi antennas (Plettemeier, 2009) in combination with Hubra 
        Returned is the (radar_volume_cross_section x total_power x extinction_function)!
        '''
        self.info('Compensating for radar system parameters ...')
        pi = np.pi
        wl = self.cmed / 1.75
        G = 10.0 ** (7.0 / 10.0)  # antenna gain of 7 dB, from Plettemeier (2009) Fig.13 @ 1.75 GHz
        beta = 25.0 * pi / 180.0  # antenna half-beamwidth in radians, from Plettemeier (2009) 3dB frop in Fig. 7 and 8
        dr = self.r[1] - self.r[0]  # depth increment
        Vs = pi * self.r ** 2 * np.sin(beta) ** 2 * dr

        r = self.r
        r.shape = (r.size, 1)
        Vs.shape = (Vs.size, 1)
        if power:
            self.data = self.data ** 2 * (4 * pi) ** 3 * r ** 4 / (wl ** 2 * G ** 2 * Vs)
        else:
            #            self.data = self.power_phase(self.data) * (4*pi)**3 * r**4 / (wl**2 * G**2 * Vs)
            self.data = self.data * np.sqrt((4 * pi) ** 3 * r ** 4 / (wl ** 2 * G ** 2 * Vs))
        self.done()


    def addCoords(self, projection="PS71", utmzone=None, gpsfile=None):
        '''Convenient function:
        Reads GPS coordinates from file and
        interpolates x,y,z coordinates in the specified projection to each trace.
        '''
        if gpsfile == None:
            gpsfile = self.filename + '.cor'
        cor = GPS(projection, utmzone=utmzone)
        cor.readCOR(gpsfile, correction='HUBRA')
        cor.interpolateProj(self.traces * self.nstacks)
        if self.nstacks > 1:
            cor.stackCoords(self.traces)
        self.coords = cor

        #        if projection == "PS71":
        #            proj_describe = 'Polar Stereographic Coordinates (lat_ts=-71)'
        #        elif projection == "UTM":
        #            proj_describe = "Universal Trans-Mercator, Zone %s" % utmzone
        #        else:
        #            proj_describe = "Lat/Long Degrees"

        self.hascoords = True


    def convertHDF(self, filein, dirout, timens=2000, stack=1, combine=False,
                   envelope=True, filerange="all-all", ww=0, cf=0,
                   projection=None, utmzone=None, gpsfile=None):
        '''
        Convenience function:
        Transfers the HUBRA raw data to time domain and
        stores the output in HDF5 file format.
        '''
        self.log()

        from processgpr.core.h5file import H5FILE

        h5file = H5FILE()
        startfile = filerange.split("-")[0]
        stopfile = filerange.split("-")[1]
        if stopfile == "all":
            stopfile = 99999999999
        else:
            stopfile = int(stopfile)

        if startfile == "all":
            filecounter = 0
            self.filename = filein
        else:
            filecounter = int(startfile)
            self.filename = filein + "continued" + str(filecounter)

        fileinname = os.path.basename(filein)
        fout = os.path.join(dirout, fileinname)

        self.readRaw()
        self.applyWindow(ww=ww, cf=cf)
        self.zeropad()
        self.applyFFT(dirout, envelope)
        if stack > 1:
            self.stackTraces(stack)
        self.savesamples = int((self.samples / self.t[-1]) * timens)
        print self.savesamples, self.samples, self.data[0:self.savesamples, :].shape
        if self.savesamples > self.samples:
            self.savesamples = self.samples
        h5file.new(fout, (self.savesamples, self.traces))
        h5file.set_metadata(self)
        h5file.append_earray_data(self.data[0:self.savesamples, :])
        if not combine:
            self.addCoords(projection=projection, utmzone=utmzone)
            h5file.append_earray_gps(self.coords)
            h5file.close()

        while (os.path.exists(filein + "continued" + str(filecounter)) and filecounter <= stopfile):
            self.filename = filein + "continued" + str(filecounter)
            print "Converting file: ", self.filename
            self.readRaw()
            self.applyWindow(ww=ww, cf=cf)
            self.zeropad()
            self.applyFFT(dirout, envelope)
            if stack > 1:
                self.stackTraces(stack)
            if not combine:
                fout = os.path.join(dirout, fileinname + "continued" + str(filecounter))
                self.savesamples = int((self.samples / self.t[-1]) * timens)
                if self.savesamples > self.samples:
                    self.savesamples = self.samples
                h5file.new(fout, (self.savesamples, 0))
                h5file.set_metadata(self)
            h5file.append_earray_data(self.data[0:self.savesamples, :])
            if not combine:
                self.addCoords(projection=projection, utmzone=utmzone)
                h5file.append_earray_gps(self.coords)
                h5file.close()
            filecounter += 1
        if combine:
            print self.traces, self.nstacks
            self.traces = h5file.earray_data.shape[1]
            print self.traces, self.nstacks
            self.addCoords(projection=projection, utmzone=utmzone,
                           gpsfile=gpsfile)
            h5file.append_earray_gps(self.coords)
            h5file.close()

        logging.info('''
        File %s processed over range %s
        Written to %s
        Time window stored: %f
        Number of traces stacked: %i
        Envelope data: %s
        Clipped frequency band - Center frequency: %f GHz, Bandwidth: %f GHz''' %
                     (filein, filerange, dirout, timens, stack, envelope, cf, ww))

        from fnmatch import fnmatch

        for file in os.listdir(os.path.dirname(os.path.abspath(filein))):
            if fnmatch(file, 'tmp%s*.dat' % os.path.basename(filein)):
                os.remove(os.path.join(os.path.dirname(os.path.abspath(filein)),
                                       file))

    def convertHDF_prestack(self, filein, dirout, timens=2000, stack=1,
                            combine=False, envelope=True, filerange="all-all",
                            ww=0, cf=0,
                            projection=None, utmzone=None, gpsfile=None):
        '''
        Convenient function:
        Transfers the HUBRA raw data to time domain and
        stores the output in HDF5 file format.
        '''
        self.log()

        from processgpr.core.h5file import H5FILE

        h5file = H5FILE()
        startfile = filerange.split("-")[0]
        stopfile = filerange.split("-")[1]
        if stopfile == "all":
            stopfile = 99999999999
        else:
            stopfile = int(stopfile)

        if startfile == "all":
            filecounter = 0
            self.filename = filein
        else:
            filecounter = int(startfile)
            self.filename = filein + "continued" + str(filecounter)

        fileinname = os.path.basename(filein)
        fout = os.path.join(dirout, fileinname)

        self.readRaw()
        if stack > 1:
            self.prestackTraces(stack)
        self.applyWindow(ww=ww, cf=cf)
        self.zeropad()
        self.applyFFT(dirout, envelope)
        self.savesamples = int((self.samples / self.t[-1]) * timens)
        print self.savesamples, self.samples, self.data[0:self.savesamples, :].shape
        if self.savesamples > self.samples:
            self.savesamples = self.samples
        h5file.new(fout, (self.savesamples, self.traces))
        h5file.set_metadata(self)
        h5file.append_earray_data(self.data[0:self.savesamples, :])
        if not combine:
            self.addCoords(projection=projection, utmzone=utmzone)
            h5file.append_earray_gps(self.coords)
            h5file.close()

        while (os.path.exists(filein + "continued" + str(filecounter)) and filecounter <= stopfile):
            self.filename = filein + "continued" + str(filecounter)
            print "Converting file: ", self.filename
            self.readRaw()
            if stack > 1:
                self.prestackTraces(stack)
            self.applyWindow(ww=ww, cf=cf)
            self.zeropad()
            self.applyFFT(dirout, envelope)
            if not combine:
                fout = os.path.join(dirout, fileinname + "continued" + str(filecounter))
                self.savesamples = int((self.samples / self.t[-1]) * timens)
                if self.savesamples > self.samples:
                    self.savesamples = self.samples
                h5file.new(fout, (self.savesamples, 0))
                h5file.set_metadata(self)
            h5file.append_earray_data(self.data[0:self.savesamples, :])
            if not combine:
                self.addCoords(projection=projection, utmzone=utmzone)
                h5file.append_earray_gps(self.coords)
                h5file.close()
            filecounter += 1
        if combine:
            print self.traces, self.nstacks
            self.traces = h5file.earray_data.shape[1]
            print self.traces, self.nstacks
            self.addCoords(projection=projection, utmzone=utmzone,
                           gpsfile=gpsfile)
            h5file.append_earray_gps(self.coords)
            h5file.close()

        #            self.addCoords(projection="PS71")#, utmzone="33")
        #            h5file.append_earray_gps(self.coords)
        #            h5file.close()

        logging.info('''
        File %s processed over range %s
        Written to %s
        Time window stored: %f
        Number of traces stacked: %i
        Envelope data: %s
        Clipped frequency band - Center frequency: %f GHz, Bandwidth: %f GHz''' %
                     (filein, filerange, dirout, timens, stack, envelope, cf, ww))

        from fnmatch import fnmatch

        for file in os.listdir(os.path.dirname(os.path.abspath(filein))):
            if fnmatch(file, 'tmp%s*.dat' % os.path.basename(filein)):
                os.remove(os.path.join(os.path.dirname(os.path.abspath(filein)),
                                       file))
          
