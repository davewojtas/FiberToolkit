import numpy as np
import sys
from scipy.signal import lfilter

class frame_orientate():
    def __init__(self, array):
        self.array = array
	self.r_map = np.zeros(self.array.shape)
        self.orientation = 0
        self.is_single = 'False' 

    def gen_radave(self):
        '''
	dsdsd
	'''
	x,y = np.indices(self.array.shape)
	center = np.array([(x.max()-x.min())/2.0, (x.max()-x.min())/2.0])
	x = x - center[0]
	y = y - center[1]
	r = np.hypot(x, y).astype(int)
	#pixelAngles = np.degrees(np.arctan2(y,x))
	tbin = np.bincount(r.ravel(), self.array.ravel())
	nr = np.bincount(r.ravel())
	rad_ave = tbin / nr
	for i in range(0, r.max()+1):
		self.r_map[r==i] = rad_ave[i]
 
    def angular_profile_rayave(self, ray_width=1, ang_incr=1, ray_ir=0, ray_or=70):
        '''
        Returns the angular profile of pixel averages within rays.
        Arguments:
            ray_width = Width in pixels of the ray to be summed over (default=1)
            ang_incr = Angular increment of ray iterative rotation (default=1)
            ray_ir = Ray start position as pixels from center (default=0)
            ray_or = Ray end position as pixels from center (default=70)
        Returns:
            2 x 360 / ang_incr array of angles and pixel averages
        '''
	self.gen_radave()
	self.array = self.array / self.r_map
	x,y = np.indices(self.array.shape)
	center = np.array([(x.max()-x.min())/2.0, (x.max()-x.min())/2.0])
	x = x - center[0]
	y = y - center[1]
	r = np.hypot(x, y)
	px_angle = np.degrees(np.arctan2(y,x))
	ray_hw = np.floor(ray_width / 2)
	n_rays = 180/ang_incr + 1
	r_s = center[0] - ray_ir # ray starting x index
	r_f = center[0] - ray_or # ray finishing x index
	ray_x = x[(center[0]-ray_hw):(center[0]+ray_hw+1), r_f:r_s]
	ray_y = y[(center[0]-ray_hw):(center[0]+ray_hw+1), r_f:r_s]
	ray_r = np.sqrt(ray_x**2 + ray_y**2)  
	ray_ang = np.arctan2(ray_y,ray_x)
	output = np.zeros((2,n_rays))
	for i in range(n_rays): 
		angle = np.pi * i * ang_incr / 180
		ray_x_rt = np.round((ray_r) * np.cos(ray_ang + angle) + center[0]) 
		ray_y_rt = np.round((ray_r) * np.sin(ray_ang + angle) + center[1]) 
		ray_sum = np.sum(self.array[ray_x_rt.astype(int),ray_y_rt.astype(int)])
		output[0, i] = 180 * angle/np.pi
		output[1, i] = ray_sum / len(ray_x)
	return output

    def determine_orientation(self, kern_win=5, low_t=0.4 ):        
        #average angular sums
        ang_prof = self.angular_profile_rayave(1, 1, 20, 70)
        #normalize
        norm_prof = (ang_prof[1] - ang_prof[1].min()) / (ang_prof[1].max() - ang_prof[1].min())
        self.orientation = ang_prof[0, norm_prof == norm_prof.max()]
        #threshold
        norm_prof = norm_prof - low_t
        neg_indices = norm_prof < 0
        norm_prof[neg_indices] = 0
        norm_prof = norm_prof / ( 1 - low_t)
        #window smooth
        s=np.r_[norm_prof[kern_win-1:0:-1], norm_prof, norm_prof[-1:-kern_win:-1]]
        w=np.ones(kern_win,'d')
        b=np.convolve(w/w.sum(), s, mode='valid')
        buff=np.floor(kern_win/2)
        norm_prof = b[buff:(len(b)-buff)]
        #find maxima
        peaks=np.zeros(len(norm_prof))
        k=np.r_[True, norm_prof[1:]>norm_prof[:-1]] & np.r_[norm_prof[:-1]>norm_prof[1:], True]
        if np.sum(k) == 1:
            self.is_single = 'True'
        print np.sum(k), self.orientation, self.is_single
        

    def moving_average(self, p, lag=3):
        '''
        Calculate the moving average of the data series p using a moving window. 
        Arguments:
            p = Data series
            lag = Window size to average over (default=3)
        Returns:
            ret = Array of averages
        '''
        ret = np.cumsum(p, dtype=float)
        ret[lag:] = ret[lag:] - ret[:-lag]
        return ret[lag - 1:] / lag

    def moving_stddev(self, x, k):
	'''
	Calculate the standard deviation of a data series using a moving window.
	Arguments:
	    p = Data series
	Returns:
	    ret = Array of standard deviations
	'''
 	#Movingstd uses filter to compute the standard deviation, using
 	#the trick of std = sqrt((sum(x.^2) - n*xbar.^2)/(n-1)).

	n = len(x)
	x = x - np.mean(x)
	x2 = x**2
	A = 1;
	#backward filtering
	#B = np.ones(k)
	#ret = np.sqrt((lfilter(B / k, A, x2) - (lfilter(B / k, A, x)**2)*(1/k))/(k-1))
        #for i in range(0, (k-1)):
      	#    ret[i] = np.std(x[0:i])
	B = np.ones(2*k+1);
	ret = np.sqrt((lfilter(B/k,A,x2) - (lfilter(B/k,A,x)**2)*(1/(2*k+1)))/(2*k))
	ret[k:(n-k)] = ret[(2*k):len(ret)]
	for i in range(0,k):
      	    ret[i] = np.std(x[0:(k+i)])
	    ret[n-k+i] = np.std(x[(n-2*k+i):n])
	return ret

    def define_peak_regions(self, p, lag=30, diff=3, influence=0):
	output_mean  = self.moving_average(p,lag)
	output_stdev = self.moving_stddev(p,lag)
	nw = len(output_mean)
	print nw
	newMean  = np.zeros(nw)
	newStdev = np.zeros(nw)
	signals  = np.zeros(nw)
	newMean[lag-1]  = output_mean[lag]
	newStdev[lag-1] = output_stdev[lag]
	for i in range(lag,nw):
	   if p[i] > newMean[i-1]+diff*newStdev[i-1]:
	      newMean[i]  = (newMean[i-1]  + influence*p[i])/(1+influence)
	      newStdev[i] = (newStdev[i-1] + influence*np.sqrt((p[i]-newMean[i-1])**2))/(1+influence)
	      signals[i]  = 1
	   else:
	      newMean[i]  = (newMean[i-1]+p[i])/2
	      newStdev[i] = (newStdev[i-1] + np.sqrt((p[i]-newMean[i-1])**2))/2
	      signals[i]  = 0
	return signals
	#figure;
	#hold all;
	#plot(p, ':r', 'LineWidth', 1, 'Color', 'black');
	#plot(signals, 'LineWidth', 2, 'Color', 'blue');
	#plot(newMean, 'LineWidth', 2, 'Color', 'red');
	#plot(newMean+newStdev, 'LineWidth', 2, 'Color', 'green');


