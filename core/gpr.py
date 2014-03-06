'''
Created on 20.10.2010

@author: Karsten Mueller
'''

# Built-in modules
import sys
import os
import logging

# Third-party modules
import tables
import numpy as np
import numpy.fft as fft

# Local modules
from processgpr.core.gps import GPS
from processgpr.core.trace import TRACE



#------------------------------------------------------------------------------ 
class GPR:
    def __init__(self):
        '''
        Initialize GPR class and set global variables.
        '''
        self.eps0 = 8.85418782e-12  #: permittivity of free space [m-3 kg-1 s4 A2]
        self.mu0 = 1.2566371e-6  #: permeability  of free space [H/m]
        self.c0 = 0.299792458  #: speed of light in m/ns
        self.cmed = 0.23  #: speed in bulk medium in m/ns

        self.data = None
        self.samples = None
        self.deltaT = 0
        self.ylabel = "Time"

        self.nstacks = 1
        self.hascoords = False
        self.filename = ""
        self.cmin = 0.0
        self.cmax = 1.0
        self.procext = '_org'  #: added to the filename when saved to avoid file deletion
        self.procstring = ""
        self.proclog = []  #: lists the processing steps in detail
        self.parameters()

    def parameters(self):
        pass

    def info(self, text):
        '''
        Writes info text to standard out.
        
        @todo: Set up for logging and H5.history.
        '''
        sys.stdout.write(text)
        self.procstring = text[:-4]

    def done(self):
        sys.stdout.write(' Done!')
        sys.stdout.flush()
        print ""
        self.proclog.append(self.procstring)

    def log(self):
        HomeDir = os.path.expanduser('~')
        log_file = os.path.join(HomeDir, "processing.log")
        logging.basicConfig(filename=log_file, level=logging.INFO)

    def writeBin(self, filename):
        self.info('Writing data to file %s ...' % filename)
        self.data.tofile(filename)  # data will be written in 'double' format
        self.done()

    def writeAscii(self, filename):
        from scipy.io import write_array

        self.info('Writing data to file: ' + filename + '.asc')
        write_array(filename + ".asc", np.transpose(self.data[1:-1, :]), separator=' ', linesep='\n')
        self.done()

    def writeH5(self, filename):
        self.info('Writing data to file: %s.h5 ...' % filename)
        if np.dtype(self.data[0, 0]) != 'float32':
            self.data = self.data.astype('float32')
        h5file = tables.openFile(filename + '.h5', mode='w', title='GPR data file')
        profile = h5file.createGroup(h5file.root, 'profile', 'GPR profile ' + filename)
        h5file.createArray(profile, 'traces', self.data, 'GPR traces')

        if self.hascoords:
            coordinates = h5file.createGroup(h5file.root, 'coordinates', 'Polar Stereographic coordinates (lat_ts=-71)')
            h5file.createArray(coordinates, 'Y', self.coords.y, 'Trace y location')
            h5file.createArray(coordinates, 'X', self.coords.x, 'Trace x location')
            h5file.createArray(coordinates, 'Elevation', self.coords.elev, 'Trace elevation')
            h5file.createArray(coordinates, 'Distance', self.coords.dist,
                               'Distance from previous trace')  # NEEDS TO BE RECALCULATED AFTER SATACK

        h5file.root.profile._v_attrs.polarity = self.PP
        h5file.flush()
        print h5file
        h5file.close()
        self.done()

    def saveH5(self, filename):
        from processgpr.core.h5file import H5FILE

        h5file = H5FILE()
        h5file.new(filename, self.data.shape)
        h5file.set_metadata(self)
        h5file.append_earray_data(self.data)
        if self.hascoords:
            h5file.append_earray_gps(self.coords)
        h5file.close()

    def readH5(self, filename, tracemin=0, tracemax=-1,
               samplemin=0, samplemax=-1, roffset=0.0):
        '''
        Reads GPR data files in HDF5 format that have been created by L{writeH5}.
        
        @param filename: The name of the input file.
        @param tracemin: The first trace to used.
        @param tracemax: The last trace to be used.
        @param samplemin: The first sample to be used.
        @param samplemax: The last sample to be used.
        @param roffset: The offset in meters between I{samplemin} and the\
        antenna. 
        '''

        self.info('Reading data from file: %s ...' % filename)
        ''' Get data '''
        self.filename = os.path.splitext(filename)[0]
        h5file = tables.openFile(filename, mode='r')
        try:
            self.totalsamples, self.totaltraces = h5file.root.data.traces.shape
        except AttributeError:
            self.totalsamples, self.totaltraces = h5file.root.profile.traces.shape
        print "\n Total", self.totalsamples, self.totaltraces
        if tracemax == -1:
            tracemax = self.totaltraces
        if samplemax == -1:
            samplemax = self.totalsamples

        self.samplendx = np.arange(samplemin, samplemax)
        #        print self.samplendx.shape
        self.tracendx = np.arange(tracemin, tracemax + 1)
        #        print self.tracendx
        try:
            self.data = h5file.root.data.traces[samplemin:samplemax,
                        tracemin:tracemax]
            self.samples, self.traces = np.shape(self.data)
        except AttributeError:
            self.data = h5file.root.profile.traces[samplemin:samplemax,
                        tracemin:tracemax]  # old format
            self.samples, self.traces = np.shape(self.data)
        if h5file.root.__contains__('coordinates'):
            self.hascoords = True
            self.coords = GPS()
            self.coords.y = h5file.root.coordinates.Y[tracemin:tracemax]
            self.coords.x = h5file.root.coordinates.X[tracemin:tracemax]
            self.coords.elev = h5file.root.coordinates.Elevation[tracemin:tracemax]
            self.coords._getDistance()
            self.d = np.linspace(0, np.add.reduce(self.coords.dist), self.traces)
        else:
            print "WARNING: No coordinates found"
        try:
            self.gatefreq = h5file.root.data.traces.attrs.gatefreq
        except AttributeError:
            print "WARNING: Missing attribute - Gating frequency not defined"
        try:
            self.nstacks = h5file.root.data.traces.attrs.nstacks
        except AttributeError:
            print "WARNING: Missing attribute - Number of stacks not defined"
        try:
            self.deltaT = h5file.root.data.traces.attrs.deltaT
        except AttributeError:
            print "WARNING: Missing attribute - Time increment not defined"
        try:
            self.cmed = h5file.root.data.traces.attrs.cmed
        except AttributeError:
            print "WARNING: Missing attribute - Medium velocity not defined"
        h5file.close()
        ''' Set parameters '''
        if samplemax == -1:
            samplemax = self.samples
        #        self.t = self.deltaT*self.samplendx*1e9
        self.t = np.arange(0.0, self.deltaT * self.samples, self.deltaT) * 1e9
        self.r = (self.t * self.cmed / 2.0) + roffset
        self.done()


    def applyWindow(self, window="hanning", ww=0, cf=0):
        '''
        Apply window function to frequency domain data
        cf: the frequency the window is centered over [Hz]
        ww: the window width [Hz], if ww equals 0 the window covers the full range
        '''
        self.info("Applying %s window ..." % window)
        if window == "hanning":
            if ww == 0:
                w = np.hanning(self.numfreq)
            else:
                pos = int((cf - self.lowF) / self.deltaF)
                halfwidth = int(ww / (2.0 * self.deltaF))
                w = np.zeros(self.numfreq)
                w[pos - halfwidth:pos + halfwidth] = np.hanning(2 * halfwidth)
        elif window == "hamming":
            if ww == 0:
                w = np.hamming(self.numfreq)
            else:
                pos = int((cf - self.lowF) / self.deltaF)
                halfwidth = int(ww / (2.0 * self.deltaF))
                w = np.zeros(self.numfreq)
                w[pos - halfwidth:pos + halfwidth] = np.hamming(2 * halfwidth)
        elif window == "blackman":
            if ww == 0:
                w = np.blackman(self.numfreq)
            else:
                pos = int((cf - self.lowF) / self.deltaF)
                halfwidth = int(ww / (2.0 * self.deltaF))
                w = np.zeros(self.numfreq)
                w[pos - halfwidth:pos + halfwidth] = np.blackman(2 * halfwidth)
        self.data = self.data * w
        self.done()


    def zeropad(self):
        '''
        Zero-pad to move samples to the correct position in the FD
        '''
        self.info("Zero-padding ...")
        tmparray = np.zeros((self.traces, self.Nfft), dtype=np.float32)
        for trace in xrange(self.traces):
            tmparray[trace, self.nstart:self.nstart + self.numfreq] = self.data[trace, :]
        self.data = tmparray
        self.numfreqpadded = self.data.shape[1]
        del (tmparray)
        self.done()


    def applyFFT(self, wdir=os.getcwd(), envelope=True):
        '''
        FFT transfers the FD data to TD data
        '''
        self.info("Applying FFT ...")
        workfile = os.path.join(wdir, "tmp" + os.path.basename(self.filename) + ".dat")
        fid = open(workfile, "wb")
        if envelope:
            for trace in xrange(self.traces):
                np.abs(fft.fft(self.data[trace])).tofile(fid)  # envelope
        else:
            for trace in xrange(self.traces):
                np.real(fft.fft(self.data[trace])).tofile(fid)  # incl. phase

        fid.close()
        fid = open(workfile, "rb")
        try:  # Zero-padding applied
            self.data = np.zeros((self.numfreqpadded, self.traces), dtype=np.float32)
            for i in range(self.traces):
                fid.seek(self.numfreqpadded * i * 8)
                self.data[:, i] = np.fromfile(fid, "<f8", self.numfreqpadded)
        except AttributeError:  # Zero-padding not applied
            self.data = np.zeros((self.numfreq, self.traces), dtype=np.float32)
            for i in range(self.traces):
                fid.seek(self.numfreq * i * 8)
                self.data[:, i] = np.fromfile(fid, "<f8", self.numfreq)

        fid.close()
        self.samples = self.data.shape[0]
        self.traces = self.data.shape[1]
        self.t = np.arange(0, self.deltaT * self.samples, self.deltaT) * 1e9  #: time axis in ns
        self.r = self.t * self.cmed / 2.0  #: range axis in m
        self.done()


    def stackTraces(self, nstacks=2):
        '''
        Stack "nstacks" traces in time domain.
        '''
        self.info('Stacking profile over %s traces ...' % nstacks)
        newtraces = int(self.traces / nstacks)
        S = np.zeros((self.samples, newtraces), np.dtype(self.data[0, 0]))
        n = 0
        for i in range(newtraces):
            S[:, i] = np.mean(self.data[:, n:n + nstacks], axis=1)
            if n + nstacks <= self.traces:
                n = n + nstacks
        self.data = S
        self.traces = newtraces
        try:
            tracendxmin = int(self.tracendx.min() / nstacks)
            tracendxmax = int(self.tracendx.min() / nstacks) + newtraces
            self.tracendx = np.arange(tracendxmin, tracendxmax + 1)
        except AttributeError:
            self.tracendx = np.arange(self.traces)
        if self.hascoords:
            self.coords.stackCoords(self.traces)
        self.procext += '_stk'
        self.nstacks = nstacks
        self.done()

    def prestackTraces(self, nstacks=2):
        '''
        Pre-stack "nstacks" traces in frequency domain.
        '''
        self.info('Pre-stacking profile over %s traces ...' % nstacks)
        newtraces = int(self.traces / nstacks)
        S = np.zeros((newtraces, self.numfreq), np.dtype(self.data[0, 0]))
        n = 0
        for i in range(newtraces):
            S[i, :] = np.mean(self.data[n:n + nstacks, :], axis=0)
            if n + nstacks <= self.traces:
                n = n + nstacks
        self.data = S
        self.traces = newtraces

        #        if self.hascoords:
        #            self.coords.stackCoords(self.traces)
        self.procext += '_prestk'
        self.nstacks = nstacks
        self.done()

    def clip(self):
        '''
        Clips unwanted neg. values.
        '''
        self.data = np.clip(self.data, 0.0, self.data.max())


    def reverse_profile(self):
        ''' Flips the profile left to right '''
        self.info('Flipping profile...')
        self.data = np.fliplr(self.data)
        self.done()


    def power_phase(self, data):
        mask = data < 0.0
        tmpdata = data ** 2
        tmpdata[mask] = -tmpdata[mask]
        return tmpdata

    def backgroundRemoval(self):
        self.info('Applying background removal ...')
        med = np.mean(self.data, axis=1)
        for i in range(self.traces):
            self.data[:, i] = self.data[:, i] - med
        self.procext = '_bkg'
        self.done()

    def subtractMean(self, Lhalf=50):
        self.info('Subtracting mean over %4.0f traces...' % (Lhalf * 2))
        reftrace = self.data[:, 70]
        tmp_data = np.zeros(self.data.shape)
        med_first = np.mean(self.data[:, 0:2 * Lhalf], axis=1)
        med_last = np.mean(self.data[:, self.traces - (2 * Lhalf):self.traces], axis=1)
        med = []
        for i in range(self.traces):
            if i < Lhalf:
                tmp_data[:, i] = self.data[:, i] - med_first
            elif i > self.traces - Lhalf:
                tmp_data[:, i] = self.data[:, i] - med_last
            elif Lhalf <= i <= (self.traces - Lhalf):
                print i - Lhalf, i, i + Lhalf
                med.append(np.mean(self.data[:, (i - Lhalf):(i + Lhalf + 1)], axis=1))
                tmp_data[:, i] = self.data[:, i] - med[-1]
        self.data = tmp_data
        self.procext = '_sbt'
        import pylab

        pylab.plot(reftrace)
        #        pylab.plot(med_first)
        #        pylab.plot(med_last)
        #        pylab.plot(med[50])
        pylab.plot(self.data[:, 70])
        pylab.plot(med[70])

        pylab.show()
        self.done()

    def linearGain(self, maxgain=1.2):
        self.info('Applying linear gain (Max. Gain = %2.2f) ...' % maxgain)
        gain = np.linspace(1.0, maxgain, self.samples)
        for i in range(self.traces):
            self.data[:, i] *= gain
        self.procext = '_lg'
        self.done()

    def compensateSpread(self, a=2.0, power=True):
        self.info('Compensating for spreading losses (simple) ...')
        r = self.r
        r.shape = (r.shape[0], 1)
        if power:
            self.data = (self.data ** 2) * r ** a
        else:
            self.data = self.data * np.sqrt(r ** a)
        self.done()

    def compensateVol(self, h0=0.65, ba=61.0, power=True, view=False):
        ''' From Davis and Poznyak, IEEE Trans., 1993 '''
        self.info('Compensating for volume spreading ...')
        beamangle = np.deg2rad(ba) / 2.0
        dr = self.deltaT * self.cmed * (1.0e9 / 2.0)  # depth increment
        Vs = 2 * np.pi * (1 - np.cos(beamangle)) / 3 * ((h0 + self.r + dr) ** 3 - (h0 + self.r) ** 3)
        r = self.r
        r.shape = (r.size, 1)
        Vs.shape = (Vs.size, 1)
        gain = 1.0 / Vs
        if view == True:
            from pylab import show, plot, xlabel, ylabel

            plot(self.t, gain)
            xlabel("Time [ns]")
            ylabel("Volume spreading correction")
            show()

        # One has to make sure that the traces start at the surface and not at the antenna
        # otherwise the initial height h0 needs to be corrected.
        if power:
            self.data = (self.data ** 2 * gain)
        else:
            self.data = (self.data * gain)
        self.done()


    def compensateVolSpread(self, h0=0.65, ba=61.0, power=True, view=False):
        ''' From Davis and Poznyak, IEEE Trans., 1993 '''
        self.info('Compensating for volume spreading losses ...')
        beamangle = np.deg2rad(ba) / 2.0
        dr = self.deltaT * self.cmed * (1.0e9 / 2.0)  # depth increment
        Vs = 2 * np.pi * (1 - np.cos(beamangle)) / 3 * ((h0 + self.r + dr) ** 3 - (h0 + self.r) ** 3)
        r = self.r
        r.shape = (r.size, 1)
        Vs.shape = (Vs.size, 1)

        if view == True:
            from pylab import show, plot, xlabel, ylabel

            if power:
                plot(self.t, ((h0 + r) ** 4 / Vs))
            else:
                plot(self.t, ((h0 + r) ** 2 / np.sqrt(Vs)))
            xlabel("Time [ns]")
            ylabel("Volume spreading correction")
            show()

        # One has to make sure that the traces start at the surface and not at the antenna
        # otherwise the initial height h0 needs to be corrected.
        if power:
            self.data = (self.data ** 2 * ((h0 + r) ** 4 / Vs))
        else:
            self.data = (self.data * ((h0 + r) ** 2 / np.sqrt(Vs)))
        self.done()


    def equalizeTraces(self):
        '''Calculates the mean value of the whole array and adjust each trace to have this mean value'''
        self.info('Equalizing traces ...')
        ndx = np.argwhere(self.data != np.nan)
        sum = np.nansum(self.data)
        s = np.shape(ndx)
        allmean = sum / (s[0] * s[1])
        print allmean
        for i in range(self.traces):
            tndx = np.argwhere(self.data[:, i] != np.nan)
            tsum = np.nansum(self.data[:, i])
            ts = np.shape(tndx)
            tmean = tsum / (ts[0] * ts[1])
            diff = allmean - tmean
            #            print allmean, tmean, diff, ts
            self.data[:, i] += diff
        self.procext = '_eqz'
        self.done()


    def distanceStack(self, distint=25.0):
        '''
        Bins traces with regard to there distance from each other.
        Apply on raw data file that were converted to the time domain!
        Inputs:
        dist is an array that contains the distances between the traces - distEach output from GPSprofile.correctTraceNo()
        distint is the distance interval over which traces should be stacked 
        '''
        self.info('Averaging traces over %f meter ...' % distint)
        averaged = []
        avglat = []
        avglong = []
        avgelev = []
        m = 0
        stackint = 0.0
        #        print self.traces
        for n in range(self.traces):
            if stackint + self.coords.dist[n] < distint:
                stackint = stackint + self.coords.dist[n]
            #                print stackint
            else:
                averaged.append(np.mean(self.data[:, m:n], axis=1))
                avglat.append(np.mean(self.coords.y[m:n], axis=0))
                avglong.append(np.mean(self.coords.x[m:n], axis=0))
                avgelev.append(np.mean(self.coords.elev[m:n], axis=0))
                m = n
                stackint = 0.0
            #        print m, np.add.reduce(self.coords.dist[m:-1])
        averaged.append(np.mean(self.data[:, m:-1], axis=1))
        avglat.append(np.mean(self.coords.y[m:-1], axis=0))
        avglong.append(np.mean(self.coords.x[m:-1], axis=0))
        avgelev.append(np.mean(self.coords.elev[m:-1], axis=0))
        self.data = np.transpose(np.asarray(averaged))
        self.coords.y = np.asarray(avglat)
        self.coords.x = np.asarray(avglong)
        self.coords.elev = np.asarray(avgelev)
        self.coords._getDistance()
        self.traces = np.shape(self.data)[1]
        self.procext = '_dstk'
        self.done()


    def removeTraces(self, remove, extend=True):
        '''
        Removes traces defined in array 'remove' in Tdata!
        NEEDS to be updated to remove the correct traces in the HDF file, too.
        '''
        if extend:
            remove = np.arange(remove[0], remove[1] + 1)
        self.info('Removing %i traces ...' % len(remove))
        old_tracendx = self.tracendx
        R = np.delete(self.data, remove, axis=1)
        if self.hascoords:
            Rx = np.delete(self.coords.x, remove, axis=0)
            Ry = np.delete(self.coords.y, remove, axis=0)
            Relev = np.delete(self.coords.elev, remove, axis=0)
        self.data = R
        self.traces = np.shape(self.data)[1]
        tracendxmin = int(self.tracendx.min())
        tracendxmax = tracendxmin + self.traces
        self.tracendx = np.arange(tracendxmin, tracendxmax + 1)
        if self.hascoords:
            self.coords.x = Rx
            self.coords.y = Ry
            self.coords.elev = Relev
            self.coords._getDistance()
        self.done()
        print "Removed traces: %i-%i" % (old_tracendx[remove[0]], old_tracendx[remove[-1]])


    def staticCorrection(self, s):
        '''
        Removes all samples above I{s}. Used to clip the data before the 
        direct wave.
        
        @param s: Number of samples to be removed from top.
        @type s: Integer
        '''
        self.info('Applying static correction ...')
        self.data = self.data[s:-1, :]
        self.samples = np.shape(self.data)[0]
        self.t = np.arange(0, self.deltaT * self.samples, self.deltaT) * self.BW
        self.r = self.t * self.cmed / 2.0
        self.procext = '_scor'
        self.done()


    def envelope(self):
        '''
        NEEDS TESTING
        
        Computes the envelope of the traces using the Hilbert transform. 
        '''
        from scipy.fftpack import hilbert

        self.info('Applying Hilbert transform ...')
        for trace in range(self.traces):
            self.data[:, trace] = np.abs(hilbert(self.data[:, trace]))
        self.done()


    def correctTopography(self, smoothit=False, sw=5):
        '''
        Places each trace at its correct elevation.
        
        @param smoothit: Bool(false), decides if the elevation values will be\
        smoothed.
        @param sw: int(5) is the windowlength of the Hanning window
        '''
        self.info('Applying topography correction ...')
        if smoothit:
            from processgpr.util.smooth import smooth

            self.coords.elev = smooth(self.coords.elev, window_len=sw, window='hanning')
        if self.hascoords:
            zmin = self.coords.elev.min() - self.r[-1]
            zmax = self.coords.elev.max()
            dz = self.r[1] - self.r[0]
            zsamples = int((zmax - zmin) / dz + 1)
            tmpdata = np.ones((zsamples, self.traces)) * -9999.0
            for n in range(self.traces):
                zspos1 = int((zmax - self.coords.elev[n]) / dz)
                zspos2 = zspos1 + self.samples
                tmpdata[zspos1:zspos2, n] = self.data[:, n]
            self.data = tmpdata
            self.samples = zsamples
            self.e = np.flipud(np.arange(zmin, zmax, dz))
            self.ylabel = "Elevation"
        else:
            print "No elevation data found!"
        self.done()

    def get_phasecenter(self, amp_threshold=1.0, freq_range=[500, 3000],
                        stacked=1, save=False, view=False):
        '''
        Calculates the phase center along the profile.
        The phase center is -delta_phi/(2 delta_k) 
        delta_k = delta_f * c0 / (2pi)
        The factor 2 in front of delta_k is due to the measuring geometry - the wave travels the distance twice 
        '''
        self.info("Calculating the phase center for each trace ...")
        self.phasecenter = np.zeros((self.traces,))
        for n in range(self.traces):
            trace = TRACE(self.data[:, n], self.deltaT, Nfft=2 ** 15)
            trace.get_phasecenter(amp_threshold, freq_range, self.cmed)
            self.phasecenter[n] = trace.zphi
            if self.hascoords:
                scoords = self.coords
        if stacked > 1:
            from processgpr.util.congrid import congrid

            newlength = int(self.phasecenter.shape[0] / stacked)
            self.phasecenter = congrid(self.phasecenter, (newlength,))
            scoords = self.coords
            scoords.stackCoords(newlength)

        if save:
            if self.hascoords:
                np.savez('zphi_coords_' + self.filename, zphi=self.phasecenter,
                         x=scoords.x, y=scoords.y, z=scoords.elev)
            else:
                np.savez('zphi_' + self.filename, zphi=self.phasecenter)

        if view:
            from pylab import show, plot

            plot(self.phasecenter)
            show()

        self.done()

    def get_penetrationdepth(self, save=False):
        """
        Calculates the penetration depth for each trace.
        """
        self.info("Calculating the penetration depth along profile ...")
        zp = []
        for trace in range(self.traces):
            T = TRACE(self.data[:, trace], self.deltaT)
            T.get_penetrationdepth(self.r)
            zp.append(T.zp)

        self.zp = np.asarray(zp, np.float32)
        if save:
            if self.hascoords:
                np.savez('zp_' + self.filename, zp=self.zp,
                         x=self.coords.x, y=self.coords.y, z=self.coords.elev)
            else:
                np.savez('zp_' + self.filename, zp=self.zp)
        self.done()

    def lin2dB(self, power=True):
        self.info('Converting to decibel ...')
        if power:
            self.data = (10 * np.log10(self.data))
        else:  # intensity needs to be squared: P = I^2
            self.data = (10 * np.log10(self.data ** 2))
        self.procext = '_dB'
        self.done()

    def dB2lin(self):
        self.info('Converting to linear power scale ...')
        self.data = 10. ** (self.data / 10.)
        self.procext = '_lin'
        self.done()

    def norm2max(self):
        self.info('Normalizing to maximum value (lin) ...')
        for n in range(self.traces):
            max_power = np.max(self.data[:, n])
            self.data[:, n] = self.data[:, n] / max_power
        self.procext = '_norm'
        self.done()

    def dBnorm2max(self):
        self.info('Normalizing to maximum value (dB) ...')
        for n in range(self.traces):
            max_power = np.max(self.data[:, n])
            self.data[:, n] = self.data[:, n] - max_power
        self.procext = '_dBnorm'
        self.done()

    def normalize(self):
        self.info('Normalizing linear scale ...')
        for n in range(self.traces):
            total_power = np.sum(self.data[:, n])
            self.data[:, n] = self.data[:, n] / total_power
        self.procext = '_norm'
        self.done()

    def dBnormalize(self):
        self.info('Normalizing dB scale ...')
        for n in range(self.traces):
            total_power = 10 * np.log10(np.add.reduce(10 ** (self.data[:, n] / 10.0)))
            self.data[:, n] = self.data[:, n] - total_power
        self.procext = '_dBnorm'
        self.done()

    def powerDecay(self):
        self.info('Estimating power decay ...')
        for n in range(self.traces):
            total_power = np.sum(self.data[:, n])
            z_power = np.cumsum(self.data[:, n])
            self.data[:, n] = (total_power - z_power) / total_power
        self.procext = '_Pdecay'
        self.done()


    def compensateExtinction(self, k0, A):
        self.info('Correcting extinction term ...')
        E = np.zeros_like(self.data)
        for n in range(self.traces):
            E[:, n] = np.exp((2 * k0[n] * self.r[:] + A[n] * self.r[:] ** 2))

        from pylab import show, imshow

        imshow(E, aspect='auto')
        show()

        self.data = 10.0 * np.log10(10.0 ** (self.data / 10.0) * E)
        self.done()

    def set_numfreq(self, val):
        self.numfreq = val
        self.deltaF = self.BW / (self.numfreq - 1)  #: sampling step in the frequency domain
        self.nstart = int(self.lowF / self.deltaF)
        self.deltaT = 1.0 / (self.deltaF * self.Nfft)  #: sampling in the time domain in [ns]
        if self.samples:
            self.t = np.arange(0, self.deltaT * self.samples, self.deltaT) * 1e9  #: time axis in ns
            self.r = self.t * self.cmed / 2.0  #: range axis in m


    def set_cmed(self, val):
        self.cmed = val
        self.r = self.t * self.cmed / 2.0
