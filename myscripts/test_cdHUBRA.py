'''
Created on 29. nov. 2014

@author: kmu
'''
import os
import numpy as np
import matplotlib.pylab as plt
from processgpr.ext.hubra import cdHUBRA
from processgpr.core.h5file import H5FILE

start = 5610
end = 5611
traces = (end-start)+1

stacked_data = np.zeros((2**15, traces)) # !!!hard-coded!!!

wdir = r'C:\Users\kmu\OneDrive\Dev\PYworkspace\chda\radar_test_data'

for i in range(start, end+1):


    file = 'data{0}'.format(i)
    pro = cdHUBRA()
    pro.filename = os.path.join(wdir, file)

    pro.readRaw()

    #plt.plot(pro.data[3,:])

    pro.applyWindow()
    pro.zeropad()
    pro.applyFFT(envelope=False)

    #unstacked_trace = pro.data[:, 30]

    pro.maxStack()
    #pro.saveH5('test_st.h5')
    #stacked_data[:, i-start] = pro.data[:, 0]
    # plt.figure()
    # plt.plot(pro.t, pro.data[:, 0])
    # plt.hold(True)
    # plt.plot(pro.t, unstacked_trace)

    if i == start:
        h5file = H5FILE()
        h5file.new(r'test_st', pro.data.shape)
        h5file.set_metadata(pro)
    h5file.append_earray_data(pro.data)

h5file.close()
plt.imshow(stacked_data, aspect='auto')


plt.show()
