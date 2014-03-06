'''
Handles GPS coordinate calculations.

Created on 20.10.2010

@author: Karsten Mueller
'''
import numpy as np

from pyproj import Proj


class GPS():
    def __init__(self, projection='PS71', **kwargs):
        self.projection = projection
        self.elev = None
        self.date = None
        self.time = None
        self.traceNo = None
        self.dist = None

        if self.projection == 'PS71':
            '''
            Initialize the profile coordinates.
            X and Y are in polar stereographic coordinates with true meridian at S71.0 and center at S90, E0
            '''
            self.y = None
            self.x = None

        elif self.projection == 'UTM':
            '''
            Initialize the profile coordinates.
            X and Y are in UTM coordinates
            REQUIRES **kwargs - utmzone
            '''
            self.y = None
            self.x = None
            try:
                self.utmzone = kwargs['utmzone']
            except:
                print ''' Missing argument 'utmzone' '''

        else:
            '''
            Initialize the profile coordinates.
            Lat and Long are in geographic coordinates
            '''
            self.Lat = None
            self.Long = None


    def readCOR(self, filename, correction=None):
        dt = np.dtype([('traceNo', np.int), ('date', 'S10'), ('time', 'S8'),
                       ('Lat', np.float), ('NS', 'S1'),
                       ('Long', np.float), ('EW', 'S1'),
                       ('Elev', np.float), ('M', 'S1'), ('factor', np.float)])
        cordata = np.loadtxt(filename, dt)
        if correction == 'HUBRA':
            i = 0
            for (lat, long) in zip(cordata['Lat'], cordata['Long']):
                cordata['Lat'][i] = float(str(lat)[:2]) + float(str(lat)[2:]) / 60.0
                if len(str(long).split('.')[0]) == 1:
                    cordata['Long'][i] = float(str(long)) / 60.0
                elif len(str(long).split('.')[0]) == 2:
                    cordata['Long'][i] = float(str(long)) / 60.0
                elif len(str(long).split('.')[0]) == 3:
                    cordata['Long'][i] = float(str(long)[:1]) + float(str(long)[1:]) / 60.0
                elif len(str(long).split('.')[0]) == 4:
                    cordata['Long'][i] = float(str(long)[:2]) + float(str(long)[2:]) / 60.0
                else:
                    cordata['Long'][i] = float(str(long)[:3]) + float(str(long)[3:]) / 60.0
                i += 1
        if cordata['NS'][0] == 'S':
            cordata['Lat'] = -cordata['Lat']
        if cordata['EW'][0] == 'W':
            cordata['Long'] = -cordata['Long']

        self.date = cordata['date']
        self.time = cordata['time']
        self.elev = cordata['Elev']
        self.traceNo = cordata['traceNo']

        if self.projection == 'UTM':
            project = Proj(proj='utm', zone=self.utmzone, ellps='WGS84')
            self.x, self.y = project(cordata['Long'], cordata['Lat'])
            self._getDistance()

        elif self.projection == 'PS71':
            project = Proj(proj='stere', ellps='WGS84', lat_0=-90.0, lon_0=0.0,
                           lat_ts=-71.0, x_0=0.0, y_0=0.0)
            self.x, self.y = project(cordata['Long'], cordata['Lat'])
            self._getDistance()

        else:
            self.Lat = cordata['Lat']
            self.Long = cordata['Long']


    def readGPS(self, filename):
        '''Read a C-band .gps file'''
        dt = np.dtype([('traceNo', np.int),
                       ('Lat', np.float), ('Long', np.float),
                       ('date', 'S10'), ('gpstime', 'S8'),
                       ('pctime', np.float)])
        gpsdata = np.loadtxt(filename, dtype=dt, skiprows=1)

        self.date = gpsdata['date']
        self.time = gpsdata['gpstime']
        self.traceNo = gpsdata['traceNo']
        self.elev = np.zeros_like(self.traceNo)

        if self.projection == 'UTM':
            project = Proj(proj='utm', zone=self.utmzone, ellps='WGS84')
            self.x, self.y = project(gpsdata['Long'], gpsdata['Lat'])
            self._getDistance()

        elif self.projection == 'PS71':
            project = Proj(proj='stere', ellps='WGS84', lat_0=-90.0, lon_0=0.0, lat_ts=-71.0, x_0=0.0, y_0=0.0)
            self.x, self.y = project(gpsdata['Long'], gpsdata['Lat'])
            self._getDistance()

        else:
            self.Lat = gpsdata['Lat']
            self.Long = gpsdata['Long']

    def readMAT(self, filename):
        from scipy.io import loadmat

        data = loadmat(filename)
        self.traceNo = data['cgps'][:, 0]
        self.x = data['cgps'][:, 2]
        self.y = data['cgps'][:, 3]
        self.elev = data['cgps'][:, 4]
        self._getDistance()

    def _getDistance(self):
        '''
        Calculate distance between traces.
        Coordinates must be in a projected system, e.g. UTM or PS
        '''
        self.dist = np.empty(self.x.shape)
        for n in range(self.x.size - 1):
            self.dist[n + 1] = np.sqrt((self.x[n + 1] - self.x[n]) ** 2 + \
                                       (self.y[n + 1] - self.y[n]) ** 2)
            if np.isnan(self.dist[n + 1]):
                # if it is NaN
                try:
                    self.dist[n + 1] = self.dist[n] + (self.dist[n] - self.dist[n - 1])
                except:
                    print "Warning: Distance array contains NaN!"
        self.total_distance = np.sum(self.dist)


    def getDistanceArray(self):
        '''
        Calculate distance between traces.
        Coordinates must be in a projected system, e.g. UTM or PS
        '''
        self.distarray = np.zeros(self.x.shape)
        for n in range(self.x.size - 1):
            self.distarray[n + 1] = self.distarray[n] + \
                                    np.sqrt((self.x[n + 1] - self.x[n]) ** 2 + (self.y[n + 1] - self.y[n]) ** 2)


    def interpolateProj(self, totaltraces):
        '''
        Adds points between two GPS points in a projected coord. sys.
        depending on self.traceNo
        '''

        tmpx = []
        tmpy = []
        tmpelev = []
        tmptraceNo = []

        # insert coords for traces 0 to first traces with GPS measurement 
        dy = self.y[1] - self.y[0]
        dx = self.x[1] - self.x[0]
        dh = self.elev[1] - self.elev[0]
        tdiff = int(self.traceNo[0])
        for i in range(tdiff):
            tmpx.append(self.x[0] - dx / tdiff * i)
            tmpy.append(self.y[0] - dy / tdiff * i)
            tmpelev.append(self.elev[0] - dh / tdiff * i)
            tmptraceNo.append(i)

        # insert coords for traces from first to last GPS measurement
        for i in range(self.x.size - 1):
            tdiff = int(self.traceNo[i + 1] - self.traceNo[i])
            dx = self.x[i + 1] - self.x[i]
            dy = self.y[i + 1] - self.y[i]
            dh = self.elev[i + 1] - self.elev[i]
            for k in range(tdiff):
                tmpx.append(self.x[i] + dx / tdiff * k)
                tmpy.append(self.y[i] + dy / tdiff * k)
                tmpelev.append(self.elev[i] - dh / tdiff * k)
                tmptraceNo.append(self.traceNo[i] + k)

        # insert coords for last GPS measurement until last trace
        dy = self.y[-2] - self.y[-1]
        dx = self.x[-2] - self.x[-1]
        dh = self.elev[-2] - self.elev[-1]
        tdiff = int(totaltraces - self.traceNo[-1])
        for i in range(tdiff):
            tmpx.append(self.x[-1] + dx / tdiff * i)
            tmpy.append(self.y[-1] + dy / tdiff * i)
            tmpelev.append(self.elev[-1] - dh / tdiff * i)
            tmptraceNo.append(self.traceNo[-1] + i)

        self.x = np.array(tmpx)
        self.y = np.array(tmpy)
        self.elev = np.array(tmpelev)
        self.traceNo = np.array(tmptraceNo)
        self._getDistance()


    def stackCoords(self, newlength):
        '''
        Recalculate the trace coordinates when stacking 'n' times.
        '''
        from processgpr.util.congrid import congrid

        self.x = congrid(self.x, (newlength,))
        self.y = congrid(self.y, (newlength,))
        self.elev = congrid(self.elev, (newlength,), method='spline')
        self.traceNo = np.arange(newlength)
        self._getDistance()


def test1():
    import pylab

    cor = GPS(projection='PS71')
    cor.readCOR(r'C:\hubradata\Fimbul091230\data1.cor', 'HUBRA')
    print 'total distance', cor.total_distance, len(cor.traceNo)
    #    pylab.scatter(cor.x, cor.y, c='b', marker='>', s=50)
    pylab.hold(True)
    #    cor.interpolateProj(200)
    #    print 'total distance', cor.total_distance, len(cor.traceNo)
    #    cor.stackCoords(46)
    #    print 'total distance', cor.total_distance, len(cor.traceNo)
    pylab.scatter(cor.x, cor.y, c=cor.elev, edgecolors='none')
    pylab.show()


if __name__ == '__main__':
    test1()        