import numpy as np
from scipy.signal import argrelextrema
import pyqtgraph as pq
import matplotlib.pyplot as plt
import sys

class frame_orientate():
    def __init__(self, array, det_cen):
        self.array = array
	self.r_map = np.ones(self.array.shape)
        self.orientation = 0
        self.orientation_snr = 0
        self.is_hit_s1 = False
        self.is_hit_s2 = False 
        self.n_orientations = 0
        self.peak_prof = np.zeros((1,180))
        self.ang_prof = np.zeros((2,180))
        self.center = det_cen
        self.plot_graphs = False
        pq.setConfigOption('background', 'w')
        pq.setConfigOption('foreground', 'k')

    def calc_r_map(self):
        '''
	Calculates the average radial intensity distribution as a class 2d array r_map
        Arguments:
        Returns:
	'''
	x,y = np.indices(self.array.shape)
	center = np.array([(x.max()-x.min())/2.0, (x.max()-x.min())/2.0]) 
	x = x - center[0]
	y = y - center[1]
	r = np.hypot(x, y).astype(int)
	tbin = np.bincount(r.ravel(), self.array.ravel())
	nr = np.bincount(r.ravel())
	rad_ave = tbin / (nr + 1)
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
	self.calc_r_map()
	self.array = self.array / self.r_map
        self.log_array = np.log(self.array)
        if self.plot_graphs == True:
            pq.show(self.log_array)
            #wait = input("THIS PRESS ENTER TO FINISH..")
	x,y = np.indices(self.log_array.shape)
	center = np.array([(x.max()-x.min())/2.0, (x.max()-x.min())/2.0])
	x = x - center[0]
	y = y - center[1]
	r = np.hypot(x, y)
	px_angle = np.degrees(np.arctan2(y,x))
	ray_hw = np.floor(ray_width / 2)
	n_rays = 180/ang_incr 
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
		ray_sum = np.sum(self.log_array[ray_x_rt.astype(int),ray_y_rt.astype(int)])
		output[0, i] = 180 * angle/np.pi
		output[1, i] = ray_sum / len(ray_x)
	self.ang_prof = output
        #normalize
        self.ang_prof[1] = (self.ang_prof[1] - self.ang_prof[1].min()) / (self.ang_prof[1].max() - self.ang_prof[1].min())  

    def calc_autocorrelation(self):
        autocor = np.abs(np.fft.fftshift(np.fft.ifft2(self.array)))/(len(self.array)*len(self.array))

        if self.plot_graphs == True:
            pq.show(autocor)
            wait = input("AUTO_COR: PRESS ENTER TO FINISH..")    

    def calc_autocorrelation_1d(self, angle, ray_width=1, peak_threshold=2):
        '''
        Returns the radial profile of fft intensity for a given angle.
        Arguments:
            angle = angle at which to calculate the radial intensity profile
            ray_width = Width in pixels of the ray to be summed over (default=1)
        Returns:
            2 x 360 / ang_incr array of angles and pixel averages
        '''
        ray_ir = 0
        ray_or = 70
	y,x = np.indices(self.array.shape)
	center = np.array([(x.max()-x.min())/2.0, (x.max()-x.min())/2.0])
	x = x - center[0]
	y = y - center[1]
	r = np.hypot(x, y)
	px_angle = np.degrees(np.arctan2(y,x))

	ray_hw = np.floor(ray_width / 2)
	r_s = center[0] - ray_ir # ray starting x index
	r_f = center[0] - ray_or # ray finishing x index
	ray_x = x[(center[0]-ray_hw):(center[0]+ray_hw+1), r_f:r_s]
	ray_y = y[(center[0]-ray_hw):(center[0]+ray_hw+1), r_f:r_s]
	ray_r = np.sqrt(ray_x**2 + ray_y**2)          
	ray_ang = np.arctan2(ray_y,ray_x)     
        
        ray_x_rt = np.round((ray_r) * np.cos(ray_ang + np.pi*angle/180) + center[0]) 
	ray_y_rt = np.round((ray_r) * np.sin(ray_ang + np.pi*angle/180) + center[1])
	ray_r_rt = np.sqrt((ray_x_rt- center[0])**2 + (ray_y_rt-center[1])**2).astype(int)
        ray_x_rt_pi = np.round((ray_r) * np.cos(ray_ang + np.pi*(angle+180)/180) + center[0]) 
	ray_y_rt_pi = np.round((ray_r) * np.sin(ray_ang + np.pi*(angle+180)/180) + center[1])	
        ray_r_rt_pi = np.sqrt((ray_x_rt_pi- center[0])**2 + (ray_y_rt_pi-center[1])**2).astype(int)
        ray_x_rt = np.concatenate((ray_x_rt, ray_x_rt_pi),axis=0)
        ray_y_rt = np.concatenate((ray_y_rt, ray_y_rt_pi),axis=0)	
        ray_r_rt = np.concatenate((ray_r_rt, ray_r_rt_pi),axis=0)

        ray_x_rt = ray_x_rt.astype(int)
        ray_y_rt = ray_y_rt.astype(int)
        ray_r_rt = ray_r_rt.astype(int)


	ray_I = self.array[ray_y_rt.astype(int),ray_x_rt.astype(int)]
	tbin = np.bincount(ray_r_rt.ravel(), ray_I.ravel())
	nr = np.bincount(ray_r_rt.ravel())
        n = len(nr)
        x = np.arange(n) 
        non_zero_counts = nr > 0
        nr = nr[non_zero_counts]
        tbin = tbin[non_zero_counts]
        x = x[non_zero_counts]
        radial_prof = tbin / nr              

        #smooth radial profile and remove background trend
        kern_win = 15
        s=np.r_[radial_prof[kern_win-1:0:-1], radial_prof, radial_prof[-1:-kern_win:-1]]
        w=np.ones(kern_win,'d')
        b=np.convolve(w/w.sum(), s, mode='valfinid')
        buff=np.floor(kern_win/2)
        radial_prof = radial_prof - b[buff:(len(b)-buff)]
        kern_win = 5
        s=np.r_[radial_prof[kern_win-1:0:-1], radial_prof, radial_prof[-1:-kern_win:-1]]
        w=np.ones(kern_win,'d')
        b=np.convolve(w/w.sum(), s, mode='valfinid')
        buff=np.floor(kern_win/2)
        radial_prof = b[buff:(len(b)-buff)] 

        autocor = 2*np.abs(np.fft.ifft(radial_prof))/(len(radial_prof))

        output = np.zeros((len(autocor),2))
	output[:,0] = x 
	output[:,1] = autocor 

        maxima=argrelextrema(autocor[0:len(autocor)/2], np.greater, mode='wrap')
        n_maxima = 0
        print maxima
        maxima = np.asarray(maxima[0])
        maxima = maxima[output[maxima, 0] >= 0]
        for peak_angle in maxima:
            n_maxima += 1


        glob_max = np.max(maxima)
        print maxima
        print glob_max
        if autocor[7] >= peak_threshold:
            #print maxima
            #print autocor[7]
            self.is_hit_s2 = True

        #print autocor

        if self.plot_graphs == True:
            plt.plot(radial_prof)
            plt.show()
            plt.plot(autocor, 'k')
            #plt.plot(output[:,0], fft_profile, 'k')
            plt.xlabel('r')
            plt.ylabel('Autocorrelation')
            plt.xlim([0,35])
            plt.show()

            #wait = input("AUTO_COR_1D: PRESS ENTER TO FINISH..")    

    def radial_profile_peak_detect(self, angle, ray_width=1, min_second_peak_amp=0.1, peak_to_mean_ratio=3):
        '''
        Returns the radial profile of fft intensity for a given angle.
        Arguments:
            angle = angle at which to calculate the radial intensity profile
            ray_width = Width in pixels of the ray to be summed over (default=1)
        Returns:
            2 x 360 / ang_incr array of angles and pixel averages
        '''
        ray_ir = 0
        ray_or = 70
        self.is_hit_s2 = False
	y,x = np.indices(self.log_array.shape)
	center = np.array([(x.max()-x.min())/2.0, (x.max()-x.min())/2.0])
	x = x - center[0]
	y = y - center[1]
	r = np.hypot(x, y)
	px_angle = np.degrees(np.arctan2(y,x))

	ray_hw = np.floor(ray_width / 2)
	r_s = center[0] - ray_ir # ray starting x index
	r_f = center[0] - ray_or # ray finishing x index
	ray_x = x[(center[0]-ray_hw):(center[0]+ray_hw+1), r_f:r_s]
	ray_y = y[(center[0]-ray_hw):(center[0]+ray_hw+1), r_f:r_s]
	ray_r = np.sqrt(ray_x**2 + ray_y**2)          
	ray_ang = np.arctan2(ray_y,ray_x)     
        
        ray_x_rt = np.round((ray_r) * np.cos(ray_ang + np.pi*angle/180) + center[0]) 
	ray_y_rt = np.round((ray_r) * np.sin(ray_ang + np.pi*angle/180) + center[1])
	ray_r_rt = np.sqrt((ray_x_rt- center[0])**2 + (ray_y_rt-center[1])**2).astype(int)
        ray_x_rt_pi = np.round((ray_r) * np.cos(ray_ang + np.pi*(angle+180)/180) + center[0]) 
	ray_y_rt_pi = np.round((ray_r) * np.sin(ray_ang + np.pi*(angle+180)/180) + center[1])	
        ray_r_rt_pi = np.sqrt((ray_x_rt_pi- center[0])**2 + (ray_y_rt_pi-center[1])**2).astype(int)
        ray_x_rt = np.concatenate((ray_x_rt, ray_x_rt_pi),axis=0)
        ray_y_rt = np.concatenate((ray_y_rt, ray_y_rt_pi),axis=0)	
        ray_r_rt = np.concatenate((ray_r_rt, ray_r_rt_pi),axis=0)

        ray_x_rt = ray_x_rt.astype(int)
        ray_y_rt = ray_y_rt.astype(int)
        ray_r_rt = ray_r_rt.astype(int)

	ray_I = self.log_array[ray_y_rt.astype(int),ray_x_rt.astype(int)]
	tbin = np.bincount(ray_r_rt.ravel(), ray_I.ravel())
	nr = np.bincount(ray_r_rt.ravel())
        n = len(nr)
        x = np.arange(n) 
        non_zero_counts = nr > 0
        nr = nr[non_zero_counts]
        tbin = tbin[non_zero_counts]
        x = x[non_zero_counts]
        radial_prof = tbin / nr              

        #smooth radial profile and remove background trend
        kern_win = 9
        s=np.r_[radial_prof[kern_win-1:0:-1], radial_prof, radial_prof[-1:-kern_win:-1]]
        w=np.ones(kern_win,'d')
        b=np.convolve(w/w.sum(), s, mode='valfinid')
        buff=np.floor(kern_win/2)
        radial_prof = radial_prof - b[buff:(len(b)-buff)]
        kern_win = 3
        s=np.r_[radial_prof[kern_win-1:0:-1], radial_prof, radial_prof[-1:-kern_win:-1]]
        w=np.ones(kern_win,'d')
        b=np.convolve(w/w.sum(), s, mode='valfinid')
        buff=np.floor(kern_win/2)
        radial_prof = b[buff:(len(b)-buff)] 
        #radial_prof = (radial_prof - radial_prof.min()) / (radial_prof.max() - radial_prof.min())     
        

        output = np.zeros((len(radial_prof),2))
	output[:,0] = x 
	output[:,1] = radial_prof 

        fft_profile = np.abs(np.fft.fft(output[:,1])) 
        #fft_profile = (fft_profile - fft_profile.min()) / (fft_profile.max() - fft_profile.min())     

        #below not used - just for diagnostics
        offset=3
        maxima=argrelextrema(radial_prof, np.greater, mode='wrap')
        n_maxima = 0
        maxima = np.asarray(maxima[0])
        maxima = maxima[output[maxima, 0] >= offset]
        for peak_angle in maxima:
            n_maxima += 1
        #print 'Rad prof = ', output[maxima,0]
        n_mode_hits = 0
        node_spacing = 9.3 
        for mode in maxima:
            for i in range(1,4):
                if output[mode, 0] / i >= (node_spacing-1) and output[mode,0] / i <= (node_spacing+1):
                    n_mode_hits += 1
        
        #if n_mode_hits >= 3:
        #    self.is_hit_s2 = True
        #else:
        #    self.is_hit_s2 = False

        peak_delta = np.zeros(len(maxima)-1)
        peak_cnt = 0
        for i in range(len(maxima)):
            if i > 0:
                peak_delta[i-1] = output[maxima[i],0] - output[maxima[i-1],0]

        peak_delta_median = np.median(peak_delta)
        
        #if peak_delta_median >= 6 and n_mode_hits >= 3:
        #    self.is_hit_s2 = True
        #else:
        #    self.is_hit_s2 = False

        #print self.is_hit_s2
        #fft_profile = np.fft.fftshift(fft_profile)
        #smooth fft profile
        #kern_win = 3
        #s=np.r_[fft_profile[kern_win-1:0:-1], fft_profile, fft_profile[-1:-kern_win:-1]]
        #w=np.ones(kern_win,'d')
        #b=np.convolve(w/w.sum(), s, mode='valfinid')
        #buff=np.floor(kern_win/2)
        #fft_profile = b[buff:(len(b)-buff)]

        #find maxima - we insert nan for zeros in the data series so that the comparitor greater_equal can be used without detecting false maxmima from the zeros
        #zeroVals = self.peak_prof == 0
        #self.peak_prof[zeroVals] = np.nan
        offset=0
        maxima=argrelextrema(fft_profile[0:len(fft_profile)/2], np.greater, mode='wrap')
        n_maxima = 0
        maxima = np.asarray(maxima[0])
        maxima = maxima[output[maxima, 0] >= offset]
        for peak_angle in maxima:
            n_maxima += 1
        
        fftprof_mean = (np.mean(fft_profile[0:5]) + np.mean(fft_profile[15:len(fft_profile)/2]))/2
        fftprof_mean = np.mean(fft_profile)

        #print 'fft mean = ', fftprof_mean 
        #fft_profile[fft_profile > 1] = 1
        #print np.mean(peak_delta)
        #maxima = maxima[fft_profile[maxima] < 1]
        #maxima = maxima[0:len(maxima)/2]
        peak_within_window = maxima[output[maxima, 0] >= 7]
        peak_within_window = peak_within_window[output[peak_within_window,0] <= 11]
        
        #window = fft_profile[7:12]
        #print window
        
        #max_peak_val = np.max(fft_profile[maxima[0:len(maxima)]])
        #max_peak_pos = output[fft_profile == max_peak_val, 0]
        #if max_peak_pos[0].astype(int) >= 7 and max_peak_pos[0].astype(int) <= 11 : #and max_peak_val / fftprof_mean >= peak_to_mean_ratio:
        #    self.is_hit_s2 = True
        #else:
        #    self.is_hit_s2 = False
        if len(peak_within_window) > 0 :
            window_peak_val = fft_profile[peak_within_window[0]]
            window_peak_pos = output[peak_within_window[0],0]
            window_peak_ratio = window_peak_val / 1 #fftprof_mean
            #sm = window_peak_val / fftprof_mean
            #print 'r2 = ', window_peak_ratio
            if window_peak_pos.astype(int) >= 7 and window_peak_pos.astype(int) <= 11 and window_peak_ratio >= peak_to_mean_ratio:
                self.is_hit_s2 = True
                #print window_peak_ratio
            else:
                self.is_hit_s2 = False

        if self.plot_graphs == True:
            print 'mode hits=', n_mode_hits
            print 'delta=', peak_delta, peak_delta_median
            #print max_peak_pos[0], max_peak_val
            print window_peak_pos, window_peak_val
            print output[maxima,0]
            #win = pq.GraphicsWindow()
            #p1 = win.addPlot(title="Radial profile (ave. subtracted)")
            #p2 = win.addPlot(title="FFT of radial profile")
            #p1.plot(output[:,0], output[:,1])
            #p2.plot(output[:,0], fft_profile)

            #pq.show(self.array)

            #plotWidget = pq.plot()
            #plotWidget.plot(output[:,0], output[:,1], pen='k')
            #plotWidget.setLabel('left', "Average intensity")
            #plotWidget.setLabel('bottom', "R")
            #plotWidget = pq.plot()
            #plotWidget.plot(output[:,0], fft_profile, pen='k')
            #plotWidget.setLabel('left', "FFT")
            #plotWidget.setLabel('bottom', "r")
            #plotWidget.setXRange(0, 35, padding=0)

            plt.plot(output[:,0], output[:,1], 'k')
            plt.xlabel('r')
            plt.ylabel('F(r)')
            plt.show()

            #plt.plot(output[:,0], fft_profile, 'k')
            #plt.xlabel('r')
            #plt.ylabel('F(F(r)')(r'$\theta$')
            #plt.xlim([0,35])
            #plt.show()

            #plt.setXRange(0, 35, padding=0)
            #pq.plot(output[:,0], output[:,1])
            #wait = input("PRESS ENTER TO FINISH..")    
            #plotWidget = pq.plot(title="Angular profiles")
            #plotWidget.plot(self.ang_prof[0], self.ang_prof[1],pen='b')
            #plotWidget.plot(self.ang_prof[0,:], self.peak_prof,pen='c')
            #plotWidget.plot(self.ang_prof[0], 0*self.ang_prof[0]+ prof_mean,pen='r')
            #plotWidget.plot(self.ang_prof[0,:], peaks_t, pen='g') 

	self.rad_prof= output #radial_prof.astype(float)
        #if max_peak_pos[0].astype(int) >= 7 and max_peak_pos[0].astype(int) <= 10 and peak_delta_median >= 6 and n_mode_hits >= 2:
        #    self.is_hit_s2 = True
        #else:
        #    self.is_hit_s2 = False

        #print self.rad_prof.shape
        #output
        #normalize
        #self.rad_prof[1] = (self.rad_prof[1] - self.rad_prof[1].min()) / (self.rad_prof[1].max() - self.rad_prof[1].min())
        # introduce r_binning
        

    def determine_orientation(self, kern_win=5, peak_to_mean_ratio=3, ann_ir=15, ann_or=50):        
        '''
        For a given frame, the prominant orientation angle and whether there are other orientation angles is determined.
        Arguments:
            kern_win = size of kernal window convolved with data for smoothing
            low_t = lower threshold value 0<low_t<1
        Returns:
        '''
        #average angular sums
        self.angular_profile_rayave(1, 1, ann_ir, ann_or)
        ang_maxsignal = self.ang_prof[0, self.ang_prof[1] == self.ang_prof[1].max()]
        if len(ang_maxsignal) > 1:
            return
        
        mean_signal = np.mean(self.ang_prof)
        #print ang_maxsignal
        #threshold
        #low_t = np.mean(self.ang_prof[1])
        self.peak_prof = self.ang_prof[1]
        #neg_indices = self.peak_prof < 0
        #self.peak_prof[neg_indices] = 0
        #self.peak_prof = self.peak_prof / ( 1 - low_t)  
        #filter special cases
        if ang_maxsignal == 0 or ang_maxsignal == 90 :
            return  
        self.orientation = 180 - ang_maxsignal
        
        #if ang_maxsignal >= 90:
        #    self.orientation = ang_maxsignal - 90
        #else:
        #    self.orientation = ang_maxsignal + 90
        #window smooth
        s=np.r_[self.peak_prof[kern_win-1:0:-1], self.peak_prof, self.peak_prof[-1:-kern_win:-1]]
        w=np.ones(kern_win,'d')
        b=np.convolve(w/w.sum(), s, mode='valid')
        buff=np.floor(kern_win/2)
        self.peak_prof = b[buff:(len(b)-buff)]      
        #prof_mean = np.mean(self.ang_prof[1])
        prof_mean = np.mean(self.peak_prof)
        #plotWidget.plot(output[:,0], output[:,1])
        #win = pq.GraphicsWindow()
        #p1 = win.addPlot()
        #p2 = win.addPlot()
        #p1.plot(output[:,0], output[:,1])
        #p2.plot(output[:,0], fft_profile)
        #pq.show(self.array)
        #pq.plot(output[:,0], output[:,1])
        #find maxima - we insert nan for zeros in the data series so that the comparitor greater_equal can be used without detecting false maxmima from the zeros
        #zeroVals = self.peak_prof == 0
        #self.peak_prof[zeroVals] = np.nan
        xdata = self.ang_prof[0]
        ydata = self.peak_prof
        ydata = (ydata - np.min(ydata)) /(np.max(ydata) - np.min(ydata))
        prof_mean = np.mean(ydata)
        self.orientation_snr = 1 / prof_mean;
        
        peaks = np.copy(ydata) - prof_mean
        neg_indices = peaks < 0 
        peaks[neg_indices] = 0 #low threshold to remove small maxima
        #peaks = peaks / ( 1 - prof_mean) #normalization - not necessary
        #filter special cases
        maxima=argrelextrema(peaks, np.greater, mode='wrap')
        self.n_orientations = 0
        peaks_t = 0* peaks
        #peaks_t = peaks_t.astype(np.float64)
        for peak_angle in maxima[0]:
            if peak_angle != 0 & peak_angle != 90 :
                peak_snr = ydata[peak_angle] / prof_mean 
                #print ydata[peak_angle], prof_mean, peak_snr
                if peak_snr >= peak_to_mean_ratio:
                    self.is_hit_s1 = True
                if peak_snr >= peak_to_mean_ratio - 1:
                    self.n_orientations = self.n_orientations + 1     
                    peaks_t[peak_angle.astype(int)] = 1
        
        #if self.is_hit_s1 == True:
            #self.is_hit_s1 = True
            #calculate average angular seperation
        #    peak_indices = self.ang_prof[0,peaks_t == 1]
        #    orientation_angles = self.ang_prof[0,peak_indices.astype(int)]
            #print orientation_angles
        #    self.mean_ang_separation = 0
        #    for i in range (1, len(orientation_angles)):
        #        self.mean_ang_separation += np.abs(orientation_angles[i-1] - orientation_angles[i])
        #    if len(orientation_angles) > 1:
        #        self.mean_ang_separation /= (len(orientation_angles)-1)
        #else:
        #    self.is_hit_s1 = False
        #xdata = self.ang_prof[0]
        #ydata = self.peak_prof
        #ydata = (ydata - np.min(ydata)) /(np.max(ydata) - np.min(ydata))
        #self.orientation_snr = 1 / np.mean(ydata);
        if self.plot_graphs == True:
            print self.is_hit_s1

            #plotWidget = pq.plot()
            #plotWidget.plot(self.ang_prof[0], self.ang_prof[1],pen='b')
 
            #plotWidget.plot(xdata, ydata, pen='k')
            #plotWidget.plot(xdata, 0*xdata + np.mean(ydata), pen=(150,150,150))
            
            #plotWidget.setLabel('left', "Average intensity")
            #plotWidget.setLabel('bottom', "Polar Angle")
            mean_line = 0*xdata + np.mean(ydata)
            plt.figure()
            plt.plot(xdata,ydata, 'k', xdata, mean_line, 'k--')
            plt.ylabel(r'I($\theta$)')
            plt.xlabel(r'$\theta$')
            plt.show()
            #plotWidget.setXRange(0, 180, padding=0)
            #plotWidget.setYRange(0, 1, padding=0)
            #plotWidget.plot(self.ang_prof[0,:], peaks_t, pen='g')
        #wait = input("PRESS ENTER TO FINISH..")  
        #print self.n_orientations, self.orientation, self.is_single

   
