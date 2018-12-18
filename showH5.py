
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as P
import math as m
import h5py


h5file = sys.argv[1]

#h5file = 'hdf5_cxif5315/r0107-eqarc/r0107-detector0-class1-sum.h5'
dataGroup = 'data/data'

h5import = h5py.File(h5file, 'r')

h5data = np.array(h5import[dataGroup])

plt.imshow(h5data, cmap='gnuplot', vmin = 0)
plt.colorbar()
plt.show()

