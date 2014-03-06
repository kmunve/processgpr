'''
Defines functionality for a single GPR trace.

Created on 20.10.2010

@author: Karsten Mueller
'''
# Third-party modules
import numpy as np
import numpy.fft as fft


class TRACE():
    def __init__(self, trace, deltaT, Nfft=None):
        '''
        Initializes the trace class.
        '''
        self.values = trace
        self.samples = len(self.values)
        self.deltaT = deltaT
        if Nfft == None:
            self.Nfft = self.samples
        else:
            self.Nfft = Nfft


    def envelope(self):
        '''
        Computes the envelope of the trace. 
        '''
        from scipy.fftpack import hilbert

        self.values = 20 * np.log10(np.abs(hilbert(self.values)))


    def subtractDC(self):
        '''
        Subtracts the DC shift from the trace.
        '''
        DCshift = np.mean(self.values[int(self.samples * 0.2):self.samples], axis=0)
        self.values = self.values - DCshift


    def spectrum(self):
        '''
        Computes the frequency spectrum of the trace.
        '''
        import pylab

        freq = fft.fftfreq(self.samples, self.deltaT)[:self.samples / 2.0]
        print freq.shape, self.samples
        f = fft.fft(self.values, self.samples) / (self.samples / 2.0)
        magnitude = 10 * np.log10(np.abs(f)[:self.samples / 2.0])
        power = (np.abs(f)[:self.samples / 2.0]) ** 2
        phase = np.angle(f)[:self.samples / 2.0]
        fig = pylab.figure()
        amp = fig.add_subplot(311)
        pow = fig.add_subplot(312)
        pha = fig.add_subplot(313)
        # Amplitude spectrum
        amp.plot(freq, magnitude)
        # Power spectrum
        pow.plot(freq, power)
        # Phase spectrum
        pha.plot(freq, phase)

    def cspectrum(self):
        '''
        Computes the complex spectrum
        '''

        freq = fft.fftfreq(self.Nfft, self.deltaT)
        Y = fft.fft(self.values, self.Nfft)
        f = freq[0:self.Nfft / 2] * 1e-6  # [MHz]
        y1 = np.abs(Y)[0:self.Nfft / 2]
        mag = y1 / y1.max()  # normalized amplitude spectrum
        pha = np.unwrap(np.angle(Y)[0:self.Nfft / 2])  # unwrapped phase spectrum

        self.freq = f
        self.mag = mag
        self.pha = pha


    def get_phasecenter(self, amp_threshold=1.0,
                        freq_range=[500, 3000],
                        c=0.299792458, view=False):
        '''
        Calculates the phase center along the profile.
        The phase center is -delta_phi/(2 delta_k) 
        delta_k = delta_f * c0 / (2pi)
        The factor 2 in front of delta_k is due to the measuring geometry - the wave travels the distance twice 
        '''
        c = c * 1e9
        self.cspectrum()
        if amp_threshold < 1.0:
            # Select range depending on the normalized amplitude spectrum
            ndx = np.where(self.mag > amp_threshold)[0]
            ndx1 = ndx[0]
            ndx2 = ndx[-1]
        else:
            # Select range depending on the frequency
            ndx1 = np.where(self.freq > freq_range[0])[0][0]
            ndx2 = np.where(self.freq > freq_range[1])[0][0]
        p1 = np.polyfit(self.freq[ndx1:ndx2] * 1e6, self.pha[ndx1:ndx2], 1)
        self.zphi = -p1[0] * c / (4 * np.pi)


    def get_penetrationdepth(self, r=None, fit=False, view=False):
        """
        Calculates the penetration depth of the trace.
        """
        from scipy.optimize import leastsq

        if r == None:
            pass

        total_power = np.sum(self.values)
        z_power = np.cumsum(self.values)
        y = (total_power - z_power) / total_power

        try:
            ndx1 = np.where(y <= np.exp(-2))[0][0]
            self.zp = r[ndx1]  # penetration depth
        except IndexError:
            print "WARNING: The curve does not decay by 1/e in the given time window!"
            self.zp = [-1.0]  # penetration depth

        if fit:
            fitfunc = lambda p: np.exp(-2 * p[0] * r - p[1] * r ** 2)  # Target function
            errfunc = lambda p: fitfunc(p)[:, 0] - y  # Distance to the target function
            p = np.asarray([-0.005, 0.005], np.float32)  # Initial guess for the parameter

            param, success = leastsq(errfunc, p)  #@UnusedVariable
            bestfit = fitfunc(param)[:, 0]
            #        print param, "-", success
            #        print "%.3e + %.3e z" % (param[0], param[1])
            pn = np.poly1d([param[1] * 0.5,
                            param[0],
                            -1])
            [root1, root2] = pn.r
            if root1 > 0.0 and root2 > 0:
                if root1 > root2:
                    self.penetration_depth = root2
                else:
                    self.penetration_depth = root1
            else:
                if root1 > root2:
                    self.penetration_depth = root1
                else:
                    self.penetration_depth = root2

            try:
                ndx2 = np.where(bestfit <= np.exp(-2))[0][0]
                self.mzp = r[ndx2]  # modeled penetration depth
            except IndexError:
                print "WARNING: The curve does not decay by 1/e in the given time window!"

                self.mzp = [-1.0]  # modeled penetration depth
            cc = np.corrcoef(y, bestfit)[0, 1]
            print "R = %2.2f" % cc

        if view:
            from pylab import plot, scatter, hold, show

            plot(r, y, c="0.5")
            hold(True)
            plot(r, bestfit, c='r')
            scatter(r[ndx1], y[ndx1], marker=">")
            scatter(r[ndx2], bestfit[ndx2], marker="x")
            show()
        
