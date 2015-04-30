__author__ = 'kmu'

"""
Storing GPR traces in an SQLite database.

Check out:
http://stackoverflow.com/questions/18621513/python-insert-numpy-array-into-sqlite3-database
"""

import sqlite3
import numpy as np
import io


# Defining the trace object
class TraceDb():

    def __init__(self, values, deltaT=0, Nfft=None):
        """
        Init the TraceDb object.

        :param values: a data array
        :param deltaT: time step between samples in ns (future use)
        :param Nfft: number of samples used for FFT (future use)
        :return:
        """
        self.values = values
        self.deltaT = deltaT
        self.direction = None # either up, down, cross_up, cross_down
        self.quality = 0

    def __str__(self):
        """

        :return: string output of TraceDb object
        """
        return "Direction: {0} [Q: {2}]\n{1}".format(self.direction, self.values, self.quality)

    def set_direction(self, d):
        """
        Sets the direction parameter
        :param d: string equal to "up", "down", "cross_up", or "cross_down"
        :return:
        """
        directions = ["up", "down", "cross_up", "cross_down"]
        if d in directions:
            self.direction = d
        else:
            print "Choose one o the following directions:"
            for direc in directions:
                print direc


# Preparing sqlite to store and retrieve arrays
def adapt_array(arr):
    out = io.BytesIO()
    np.save(out, arr)
    out.seek(0)
    return buffer(out.read())


def convert_array(text):
    out = io.BytesIO(text)
    out.seek(0)
    return np.load(out)

# Converts np.array to TEXT when inserting
sqlite3.register_adapter(np.ndarray, adapt_array)

# Converts TEXT to np.array when selecting
sqlite3.register_converter("ARRAY", convert_array)

################
# Doing a test #
################

# Creating a test TraceDb object
val = np.exp(np.linspace(-1., -4.0, 50)) * np.sin(np.linspace(0.0, 10*np.pi, 50))
trace = TraceDb(val)
trace.set_direction("up")

# Creating and connecting to a database
con = sqlite3.connect("test.db", detect_types=sqlite3.PARSE_DECLTYPES)
cur = con.cursor()

# Removing an existing table with name "test"
cur.execute("DROP TABLE IF EXISTS test")

# Adding a table called "test" to the DB
cur.execute("create table test(Id INTEGER PRIMARY KEY AUTOINCREMENT, direction TEXT, val ARRAY, quality INT);")

# Inserting three traces to test ID-autoincrement
cur.execute("INSERT INTO test(direction, val, quality) VALUES(?, ?, ?)", (trace.direction, trace.values, trace.quality))
cur.execute("INSERT INTO test(direction, val, quality) VALUES(?, ?, ?)", (trace.direction, trace.values, trace.quality))
cur.execute("INSERT INTO test(direction, val, quality) VALUES(?, ?, ?)", (trace.direction, trace.values, trace.quality))

# Commiting the changes to the DB
con.commit()

# Retrieving data from the DB
cur.execute("select * from test")
data = cur.fetchone()[2]

# Check output
print(data)
print(type(data))
if data.all() == trace.values.all():
    print("YES")





