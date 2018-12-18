import get_events as ge
import matplotlib.pyplot as plt
from scipy.signal import argrelmax
import frame_orientate as fo
import numpy as np
import sys

# run131 Working for ...
# single 141, 941, 54, 130
# multiple 12, 2843
# Not working for ...
# single
# multiple
e=ge.Events(run=131)
e.events(1, start=1087, roi=True) 

fd = e.do_fft(e.frames[0])
rd = e.sum_roi(e.frames[0])
#hits = e.get_hitframes('True')

f=fo.frame_orientate(np.log(fd))

f.determine_orientation()

#p=f.angular_profile_rayave(1, 1, 20, 70)
# normalize profile
#p[1] = (p[1] - p[1].min() ) / (p[1].max() - p[1].min())

#use moving_stddev to find peak regions
#lag = 30;       # lag of the algorithm
#diff = 3.5;     # number of standard deviations from the mean
#infl = 0.0;  # influence parameter (convergence to signals)
#prof = np.zeros(len(p[1]) + lag -1)
#prof[0:len(p[1])] = p[1]
#prof[len(p[1]):len(prof)] = p[1, 0:(lag-1)] 
#signal_mask=f.define_peak_regions(prof, lag, diff, infl)
#signal_mask=np.ones(len(p[0])) 

#find max in series
#b = p[1]
#globMax = p[0,b == b.max()]

#lower threshold
#lt = 0.4
#b = (b - lt)
#neg_indices = b < 0
#b[neg_indices] = 0
#b = b/(1 - lt)

#smooth
#window_len=5
#s=np.r_[b[window_len-1:0:-1],b,b[-1:-window_len:-1]]
#w=np.ones(window_len,'d')
#b=np.convolve(w/w.sum(),s,mode='valid')
#buff = np.floor(window_len/2)
#s=b[buff:(len(b)-buff)]
#b=s

#apply signal mask
#b = b*signal_mask

#lower threshold
#lt = 0.1
#b = (b - lt)
#neg_indices = b < 0
#b[neg_indices] = 0
#b = b/(1 - lt)

#find peaks
#pk = np.zeros(len(b))
#k = np.r_[True, b[1:] > b[:-1]] & np.r_[b[:-1] > b[1:], True]
#pk[k]= b[k]

#apply angle exclusion rules relative to the series max
#print pk


ax1 = plt.subplot2grid((2,3), (0,0))
ax1.imshow(rd, cmap='gnuplot2')
ax1.set_title('ROI sum')
ax2 = plt.subplot2grid((2,3), (1,0))
ax2.imshow(f.array, cmap='gnuplot2')
ax2.set_title('log(fft(ROI))')
ax3 = plt.subplot2grid((2,3), (0,1), colspan=2, rowspan=2)
plt.plot(f.ang_prof[0], f.ang_prof[1], 'b-', f.ang_prof[0], f.peak_prof, 'g-')
ax3.set_title('norm. px average vs polar angle')
ax3.set_xlim([0,180])
plt.show()




