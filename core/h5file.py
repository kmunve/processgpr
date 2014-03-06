'''
Defines the standard file format in I{processgpr}.

Created on 20.10.2010

@todo: Make name of coordinate system variable - default="unknown"

@author: Karsten Mueller
'''
# Third-party modules
from tables import openFile, Float32Atom  #@UnresolvedImport


class H5FILE():
    def open(self):
        pass

    def close(self):
        self.h5file.close()

    def display(self):
        print self.h5file

    def new(self, filename=None, dims=(1024, 5000), gpsdims=(0,)):
        expected = dims[1]
        dims = (dims[0], 0)
        self.h5file = openFile(filename + '.h5', mode='w',
                               title='GPR data file')
        self.group_data = self.h5file.createGroup(self.h5file.root, 'data',
                                                  'GPR data')
        self.earray_data = self.h5file.createEArray(self.group_data, 'traces',
                                                    Float32Atom(), dims,
                                                    expectedrows=expected)
        self.earray_data.attrs.processing_log = ""

        self.coordinates = self.h5file.createGroup(self.h5file.root,
                                                   'coordinates',
                                                   'x,y,z coordinates of the traces')

        self.earray_x = self.h5file.createEArray(self.coordinates, 'X',
                                                 Float32Atom(), gpsdims,
                                                 expectedrows=5000)
        self.earray_y = self.h5file.createEArray(self.coordinates, 'Y',
                                                 Float32Atom(), gpsdims,
                                                 expectedrows=5000)
        self.earray_elev = self.h5file.createEArray(self.coordinates,
                                                    'Elevation', Float32Atom(),
                                                    gpsdims, expectedrows=5000)
        self.earray_dist = self.h5file.createEArray(self.coordinates,
                                                    'Distance', Float32Atom(),
                                                    gpsdims, expectedrows=5000)  # NEEDS TO BE RECALCULATED AFTER STACK


    def append_earray_data(self, data):
        self.earray_data.append(data)
        self.h5file.flush()

    def append_earray_gps(self, coords):
        self.earray_x.append(coords.x)
        self.earray_y.append(coords.y)
        self.earray_elev.append(coords.elev)
        self.earray_dist.append(coords.dist)
        self.h5file.flush()

    def set_metadata(self, GPR):
        self.earray_data.attrs.deltaT = GPR.deltaT
        try:
            self.earray_data.attrs.gatefreq = GPR.gatefreq
        except AttributeError:
            pass
        self.earray_data.attrs.cmed = GPR.cmed
        self.earray_data.attrs.nstacks = GPR.nstacks
        self.earray_data.attrs.hascoords = GPR.hascoords
        self.earray_data.attrs.filename = GPR.filename
        for log in GPR.proclog:
            self.earray_data.attrs.processing_log += "%s\n" % log

    def save(self, GPR):
        pass
        