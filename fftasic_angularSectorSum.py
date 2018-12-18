import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as P
import math as m
import h5py
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp

# import uncorrected (geom) data
h5file = sys.argv[1] #singular/cxim2716-r0129_1120.h5 
dataGroup = 'data/data'
h5import = h5py.File(h5file, 'r')
h5data = np.array(h5import[dataGroup])

# select out a single asic
asd = 185 #approxAsicDim
cr = 5 #crop edge region width
m = 2 #asic index
n = 1 #asic index
panelRegion = h5data[((m-1)*asd+1+cr):(m*asd-cr),((n-1)*asd+1+cr):(n*asd-cr)]   

# frequency 
fpr = np.fft.fft2(panelRegion)
fpr = np.fft.fftshift(fpr)
amp = np.absolute(fpr)

# calculate power spectrum, (intensity radially averaged)
y, x = np.indices(amp.shape)
center = np.array([np.ceil((x.max()-x.min())/2.0), np.ceil((x.max()-x.min())/2.0)])
r = np.hypot(x - center[0], y - center[1])
#r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
r_int = r.astype(int)
tbin = np.bincount(r_int.ravel(), amp.ravel())
nr = np.bincount(r_int.ravel())
radial_prof = tbin / nr

#miny = min(radial_prof[60:(len(radial_prof)-5)])
#print miny
#radial_prof = radial_prof - miny # zero-base

n = len(radial_prof)
x = np.arange(n) 
y = radial_prof

a = 200000
b= -0.3
c = 12000	
		
def powerlaw(x,a,b,c):
	return a*x**b + c

def exponentialDecay(x,a,b,c):
	return a*np.exp(b*x) + c
		
popt,pcov = curve_fit(exponentialDecay,x[5:len(x)],y[5:len(x)],p0=[a,b,c])

#print popt[0], popt[1], popt[2]
#radial_prof = radial_prof - a*np.exp(b*x) + c 
#powerlawFit = popt[0] * x**popt[1] + popt[2]
subtractFit = 0
exponentialFit = popt[0]*np.exp(popt[1]*x) + popt[2]
radial_prof_mfit = radial_prof - subtractFit * exponentialFit
neg_values_indices = radial_prof_mfit < 0
radial_prof_mfit[neg_values_indices] = 0.01

# calculate polar sums of intensities within sector bins
y, x = np.indices(amp.shape)
center = np.array([np.ceil((x.max()-x.min())/2.0), np.ceil((x.max()-x.min())/2.0)])
x = x - center[0]
y = y - center[1]
r = np.hypot(x, y)

pixelAngles = np.degrees(np.arctan2(y,x)) # relative to the asic center
#pixelAngles[(y >= 0) & (x < 0)] = pixelAngles[(y >= 0) & (x < 0)] + 90
#pixelAngles[(y < 0) & (x < 0)] = pixelAngles[(y < 0) & (x < 0)] + 180
#pixelAngles[(y < 0) & (x >= 0)] = pixelAngles[(y < 0) & (x >= 0)] + 270 

print np.min(pixelAngles)

sectorSize = 2 # degrees
nSectors = 180/sectorSize
sectorSums = [0] * nSectors
sectorAngles = np.arange(nSectors) * sectorSize
print sectorAngles

for i in range(nSectors): 
	sectorIndices = [(pixelAngles > i * sectorSize) & (pixelAngles < (i+1) * sectorSize)]
	sectorSums[i] = np.sum(amp[sectorIndices])

print sectorSums
#plt.plot(sectorAngles, sectorSums)
plt.imshow(pixelAngles, cmap='Greys')#, vmin = 0)
plt.show()
sys.exit()

#radial_prof = radial_prof**2

# plot data
#ax1 = plt.subplot2grid((3,4), (0,0), colspan=3, rowspan=3 )
#ax1.imshow(h5data, cmap='gnuplot', vmin = 0)
#ax1.set_title('File=%s'%(h5file))
#ax2 = plt.subplot2grid((3,4), (0, 3))
#ax2.imshow(panelRegion, cmap='gnuplot', vmin = 0)
#ax2.set_title('Asic, D(%s,%s)'%(m,n))
#ax3 = plt.subplot2grid((3,4), (1, 3))
#ax3.imshow(np.log(amp), cmap='gnuplot')#, vmin = 0)
#ax3.set_title('log(FFT(D))')
#ax4 = plt.subplot2grid((3,4), (2, 3))
#ax4.semilogy(radial_prof)
#ax4.set_title('Power Spectrum')
#ax4.set_xlim([0,center[0]])
#plt.show()

ax1 = plt.subplot2grid((4,6), (0,1), colspan=4, rowspan=2 )
ax1.imshow(h5data, cmap='gnuplot', vmin = 0)
ax1.set_title('File=%s'%(h5file))
ax2 = plt.subplot2grid((4,6), (2, 0), colspan=2, rowspan=2 )
ax2.imshow(panelRegion, cmap='gnuplot', vmin = 0)
ax2.set_title('Asic, D(%s,%s)'%(m,n))
ax3 = plt.subplot2grid((4,6), (2, 2), colspan=2, rowspan=2 )
ax3.imshow(np.log(amp), cmap='gnuplot')#, vmin = 0)
ax3.set_title('log(Amp.)')
ax4 = plt.subplot2grid((4,6), (2, 4), colspan=2, rowspan=2 )
ax4.semilogy(radial_prof_mfit**2)
ax4.set_title('Power Spectra')
ax4.set_xlim([0,center[0]])
plt.show()





