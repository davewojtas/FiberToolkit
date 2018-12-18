import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as P
import math as m
import h5py
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
import get_events as ge
import pyqtgraph as pg

# import data from psana
e = ge.Events()
e.events(1, start=12, roi=True)
fd = e.do_fft(e.frames[0])

amp = fd

# calculate polar sums of intensities within sector bins
x,y = np.indices(amp.shape)
center = np.array([(x.max()-x.min())/2.0, (x.max()-x.min())/2.0])
x = x - center[0]
y = y - center[1]

#print center
r = np.hypot(x, y)

pixelAngles = np.degrees(np.arctan2(y,x))

rayHW = 0 # pixels
angleIncr = 1 # degrees
nRays = 360/angleIncr + 1
angleSeries = np.zeros(nRays) #[0] * nRays
sumSeries = np.zeros(nRays) #[0] * nRays

#sectorAngles = np.arange(nRays) * sectorSize
#print sectorAngles
fS = 15
fV = 70
diff  = fV - fS
r_s = center[0] - fS # ray starting x index
r_f = center[0] - fV # ray finishing x index
# construct ray mask
rayXValues = x[(center[0]-rayHW):(center[0]+rayHW+1), r_f:r_s]
#print rayXValues
rayYValues = y[(center[0]-rayHW):(center[0]+rayHW+1), r_f:r_s]
#print rayYValues
rayRValues = np.sqrt(rayXValues**2 + rayYValues**2)  
#print rayRValues
rayAngles = np.arctan2(rayYValues,rayXValues)


amp2 = np.copy(amp)

polarPlot = np.zeros((diff, nRays))

for i in range(nRays): 
	angle = np.pi * i * angleIncr / 180
	rayXValues_rot = np.round((rayRValues) * np.cos(rayAngles + angle) + center[0]) 
	rayYValues_rot = np.round((rayRValues) * np.sin(rayAngles + angle) + center[1]) 
	amp2[rayXValues_rot.astype(int),rayYValues_rot.astype(int)] = amp2.max()
	angleSeries[i] = 180 * angle/np.pi
	rayVals = amp[rayXValues_rot.astype(int),rayYValues_rot.astype(int)]
	sumSeries[i] = np.sum(rayVals)
	sumSeries[i] = sumSeries[i] / len(rayXValues)
	#polarPlot[:, i] = amp[rayXValues_rot.astype(int),rayYValues_rot.astype(int)]
	
ssMin = sumSeries.min()
ssMax = sumSeries.max()
sumSeries = (sumSeries - ssMin) / (ssMax - ssMin)


#plt.imshow(polarPlot, cmap='gray_r')
#plt.show()
#sys.exit()

#pixelAngles = pixelAngles.flatten()
#pixelAngles = pixelAngles.reshape(174,174)
#amp = ampFlatten.reshape(174,174)
#print raySums
#plt.plot(angleSeries, sumSeries)
#plt.imshow(amp2, cmap='gnuplot')#, vmin = 0)
#plt.show()
#sys.exit()

# plot data
ax1 = plt.subplot2grid((2,3), (0,0))
ax1.imshow(np.log(amp), cmap='gnuplot')
ax1.set_title('log(FFT(D)')
ax2 = plt.subplot2grid((2,3), (1,0))
ax2.imshow(np.log(amp2), cmap='gnuplot')
ax2.set_title('ray coverage')
ax3 = plt.subplot2grid((2,3), (0,1), colspan=2, rowspan=2)
ax3.plot(angleSeries, sumSeries)
ax3.set_title('norm. px average vs polar angle')
ax3.set_xlim([0,360])
plt.show()

sys.exit()




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



# polar plot of det space intensity
sectorSums = np.zeros(nRays)

for i in range(nRays): 
	sectorCount = pixelAngles[(pixelAngles > i * angleIncr) & (pixelAngles < (i+1) * angleIncr)].size
	sectorSums[i] = np.sum(amp[(pixelAngles > i * angleIncr) & (pixelAngles < (i+1) * angleIncr)])
	sectorSums[i] = sectorSums[i] / sectorCount
plt.plot(sectorSums)
plt.show()
sys.exit()
