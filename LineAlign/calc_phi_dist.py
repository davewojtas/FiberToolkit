import numpy as np
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp


def fit_gaussian(x_data, fit_data):
    '''
    Calculates best-fit gaussian distribution for the input data 
    Arguments: x, fit_data
    Returns: array of statistical firt parameters
    '''
    n = len(fit_data)
    a = np.max(fit_data)
    mean = x_data[fit_data == a]
    sigma = (mean - x_data[0])/2
    m = (fit_data[n-1] - fit_data[0]) / (x_data[n-1] - x_data[0])
    c = fit_data[0] - m * x_data[0]	
		
    #def gaussOnGradient(x,a,x0,sigma,m,c):
	#return a*exp(-(x-x0)*(x-x0)/(2*sigma*sigma)) + m*(x-x0) + c

    def gauss(x_data,a,x0,sigma):
        return a*exp(-(x_data-x0)*(x_data-x0)/(2*sigma*sigma))
		
    popt,pcov = curve_fit(gauss,x_data,fit_data,p0=[a,mean,sigma])

    return popt


#phi_gauss = fit_gaussian(phi_bin_pos, phi_line)


class Hit_Stat:
    def __init__(self, run_num, on_silicon, frame_num, orientation, peak_snr, n_orientations, px_mean):
        self.run = int(run_num)
        self.on_silicon = on_silicon
        self.frame = int(frame_num)
        self.phi = orientation
        self.phi_peak_to_mean = peak_snr
        self.n_phi_vals = int(n_orientations)
        self.roi_px_mean = px_mean
        self.sym_res_ll = 0
        self.sym_res_ave = 0
    def Print(self):        
        print self.frame, self.roi_px_mean, self.phi, '(',self.phi_peak_to_mean, self.n_phi_vals,')  ', self.sym_res_ave
    def PrintVerbose(self):        
        print self.run, self.on_silicon, self.frame, self.roi_px_mean, self.phi, '(',self.phi_peak_to_mean, self.n_phi_vals,')  '

restrict_phi = 0
phi_peak_to_mean_min = 0
plot_roi_vs_pos = 1
plot_phi_vs_pos = 0
plot_phi_hist = 0
run_numbers = [242]
phi_values = []
roi_px_values = []
frame_numbers = []
max_sym_res_values = []

for i in range(0, len(run_numbers)):
    file_name = 'hits_run'
    if run_numbers[i] < 100:
        file_name += '0'
    file_name += str(run_numbers[i]) + '.npy'
    #file_name += str(run_numbers[i]) + '_relaxed.npy'
    print 'loading', file_name, '...'
    run_hits = np.load(file_name)
    print "n_hits=", len(run_hits)
    for h in range(0, len(run_hits)):
        #run_hits[h].Print()
        if run_hits[h].phi_peak_to_mean >= phi_peak_to_mean_min:
            phi_values = np.append(phi_values, run_hits[h].phi)
            roi_px_values = np.append(roi_px_values, run_hits[h].roi_px_mean)
            frame_numbers = np.append(frame_numbers, run_hits[h].frame)
            #max_sym_res_values  = np.append(max_sym_res_values, run_hits[h].sym_res_ave)

if plot_phi_vs_pos:
    frame_numbers = frame_numbers[phi_values > 1]
    roi_px_values = roi_px_values[phi_values > 1]
    phi_values = phi_values[phi_values > 1]
    if restrict_phi:
        frame_numbers = frame_numbers[phi_values > 60]
        roi_px_values = roi_px_values[phi_values > 60]
        phi_values = phi_values[phi_values > 60]
        frame_numbers = frame_numbers[phi_values < 140]
        roi_px_values = roi_px_values[phi_values < 140]
        phi_values = phi_values[phi_values < 140]

#print frame_numbers
print "n_hits in range=", len(phi_values)
max_roi_px_value = np.max(roi_px_values)
roi_px_values[roi_px_values > max_roi_px_value/20] = max_roi_px_value / 20

if plot_roi_vs_pos:

    upscale_factor = 1
    n_windows = 21
    n_holes = 14 # per window
    n_holes_per_row = 297
    window_gap = 11 # in holes
    accel_pulses = 7
    window_numbers = frame_numbers / (n_holes_per_row)
    row_number = np.floor(window_numbers)
    print window_numbers
    window_numbers -= np.floor(window_numbers)
    #window_numbers *= (n_windows)
    window_numbers = np.floor(window_numbers)
    print np.min(window_numbers), np.max(window_numbers)
    x_pos = frame_numbers - row_number*n_holes_per_row + (window_numbers) * (window_gap) #- 3*row_number#+ np.ceil(row_number /2) * 0.5 * 
    y_pos = row_number #* np.sqrt(3) / 2

    x_pos *= upscale_factor
    y_pos *= upscale_factor

    x_max = np.max(x_pos) + 1 #n_holes_per_row #n_windows * n_holes + (n_windows - 1) * window_gap
    y_max = np.max(y_pos) + 1
    print x_max
    y_mesh, x_mesh = np.mgrid[slice(0, y_max, 1), slice(0, x_max, 1)]    

    window_mesh = 0*y_mesh
    for i in range(0,n_windows):
        x_start =  i*(n_holes + window_gap) + accel_pulses + 1
        x_fin = x_start + n_holes - 1
    
        x_start -= 1 #zero base
        x_fin -=  1 #zero base
    #    print x_start, x_fin
        window_mesh[(x_mesh >= x_start) & (x_mesh <= x_fin)] = 1

    px_min = np.min(roi_px_values)
    px_max = np.max(roi_px_values)
    
    z_mesh = 0*y_mesh + px_min + -1 # default white
    #z_mesh[window_mesh == 1] = px_max + 1 # in window, set grey
    in_window = window_mesh == 1
    print "no. in window=", len(in_window)
    
    z_mesh[y_pos.astype(int), x_pos.astype(int)] = roi_px_values 
   
    plt.figure()
    plt.axis('on')
    ax = plt.gca()
    cmap = plt.cm.gnuplot2
    #cmap.set_over("silver")
    cmap.set_under("white")
    im = ax.pcolor(x_mesh, y_mesh, z_mesh, cmap='gnuplot2', vmin = px_min, vmax=px_max)
    #print 'im=', im
    plt.xlim([0,np.max(x_mesh)])
    plt.ylim([0,np.max(y_mesh)])
    plt.xlabel('y')
    plt.ylabel('x')
    plt.axes().set_aspect('equal')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", size="15%", pad=0.2)
    #plt.colorbar(im, cax=cax)
    #plt.colorbar(im)
    cbar = plt.colorbar(im, cax=cax, orientation='horizontal')
    cbar.ax.tick_params(axis='x',direction='in', labeltop='on', labelbottom='off')  # horizontal colorbar
    plt.show()

if plot_phi_vs_pos:

    upscale_factor = 1
    n_windows = 21
    n_holes = 14 # per window
    n_holes_per_row = 297
    window_gap = 11 # in holes
    accel_pulses = 7
    window_numbers = frame_numbers / (n_holes_per_row)
    row_number = np.floor(window_numbers)
    #print window_numbers
    window_numbers -= np.floor(window_numbers)
    #window_numbers *= (n_windows)
    window_numbers = np.floor(window_numbers)
    #print np.min(window_numbers), np.max(window_numbers)
    x_pos = frame_numbers - row_number*n_holes_per_row + (window_numbers) * (window_gap) #- 3*row_number#+ np.ceil(row_number /2) * 0.5 * 
    y_pos = row_number #* np.sqrt(3) / 2

    x_pos *= upscale_factor
    y_pos *= upscale_factor

    x_max = np.max(x_pos) + 1 #n_holes_per_row #n_windows * n_holes + (n_windows - 1) * window_gap
    y_max = np.max(y_pos) + 1
    print x_max
    y_mesh, x_mesh = np.mgrid[slice(0, y_max, 1), slice(0, x_max, 1)]    

    window_mesh = 0*y_mesh
    for i in range(0,n_windows):
        x_start =  i*(n_holes + window_gap) + accel_pulses + 1
        x_fin = x_start + n_holes - 1
    
        x_start -= 1 #zero base
        x_fin -=  1 #zero base
    #    print x_start, x_fin
        window_mesh[(x_mesh >= x_start) & (x_mesh <= x_fin)] = 1
    
    z_mesh = 0*y_mesh + -1 # default white
    z_mesh[window_mesh == 1] = 182 # in window, set grey
    z_mesh[y_pos.astype(int), x_pos.astype(int)] = phi_values

    blob_region = (z_mesh > 60) & (z_mesh <180) & (x_mesh > 200)
    print "n_hits in blob=", np.sum(blob_region)
    
    plt.figure()
    plt.axis('off')
    ax = plt.gca()
    cmap = plt.cm.gnuplot2
    cmap.set_over("silver")
    cmap.set_under("white")
    im = ax.pcolor(x_mesh, y_mesh, z_mesh, cmap='gnuplot2', vmin = 0, vmax=180)
    #print 'im=', im
    plt.xlim([0,np.max(x_mesh)])
    plt.ylim([0,np.max(y_mesh)])
    plt.xlabel('y')
    plt.ylabel('x')
    plt.axes().set_aspect('equal')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", size="15%", pad=0.2)
    #plt.colorbar(im, cax=cax)
    #plt.colorbar(im)
    cbar = plt.colorbar(im, cax=cax, ticks=[0, 45, 90, 135, 180], orientation='horizontal')
    cbar.ax.tick_params(axis='x',direction='in', labeltop='on', labelbottom='off')  # horizontal colorbar
    plt.show()



#print max_sym_res_values
if plot_phi_hist:
    plt.figure(figsize=(5.2,3.9))
    fig = plt.gca()

    n_bins = 40
    
    phi_step_size = (np.max(phi_values) - np.min(phi_values) + 1)/ n_bins
    phi_bins = np.arange(np.min(phi_values), np.max(phi_values), phi_step_size)
    n, bins, patches = plt.hist(phi_values, phi_bins, facecolor='k', alpha=0.3)
    print n
    print bins
    print phi_bins
    phi_line,x = np.histogram(phi_values, phi_bins)
    phi_bin_pos = np.linspace(np.min(phi_bins), np.max(phi_bins), num=len(phi_line))

    plt.xlabel(r'$\phi\,(^\circ)$', fontsize=16)
    plt.ylabel('Count', fontsize=14)

    phi_gauss = fit_gaussian(phi_bin_pos, phi_line)

    #fig.axes.get_yaxis().set_visible(False)

    #plt.plot(phi_bin_pos, phi_line,'k-')
    #plt.ylabel('Count')
    #plt.title('Histogram of Orientations')
    #plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
    plt.axis([60, 140, 0, 40])
    plt.axis([0, 180, 0, 20]) #r049 ratio 3
    #plt.grid(True)
    plt.tight_layout()
    plt.show()
