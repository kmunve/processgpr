'''
The I{PICKS} class handles picked data.

Created on 20.10.2010

@author: Karsten Mueller
'''


class PICKS():
    def __init__(self):
        self.trace = []
        self.twt = []
        self.x = []
        self.y = []
        self.elev = []
        self.tell = 0  # counts the picks
        self.pltid = 0  # controls if picks are plotted

    def append(self, trace, twt, x, y, elev):
        self.trace.append(trace)
        self.twt.append(twt)
        self.x.append(x)
        self.y.append(y)
        self.elev.append(elev)
        self.tell += 1

    def save(self, file):
        fid = open(file, 'w')
        fid.write("#trace, twt, x, y, elev\n")
        for n in range(self.tell):
            fid.write('%i, %f, %f, %f, %f\n' % (self.trace[n], self.twt[n],
                                                self.x[n], self.y[n],
                                                self.elev[n]))
        fid.close()

    def load(self, file):
        from numpy import loadtxt, dtype, int, float

        dt = dtype([('trace', int), ('twt', float),
                    ('x', float), ('y', float), ('elev', float)])
        data = loadtxt(file, dtype=dt, delimiter=',')
        self.trace = data['trace']
        self.twt = data['twt']
        self.x = data['x']
        self.y = data['y']
        self.elev = data['elev']