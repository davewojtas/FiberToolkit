from psana import *
import numpy as np
import h5py as h5
import os.path
import sys
import pyqtgraph as pq

class Events():
    def __init__(self, ds_str='exp=cxils2616:run=67:idx', det_str='DsaCsPad', run=None, det_cen=[867,869]):
        if run is not None:
            ds_str = 'exp=cxils2616:run=%d:idx' % run
        ds = DataSource(ds_str)
        print ds_str
        self.epics = ds.env().epicsStore()
        self.calib = ds.env().calibStore()
        self.run = ds.runs().next()
        self.times = self.run.times()
        self.evt = self.run.event(self.times[0])
        self.det = Detector(det_str)
        self.ix, self.iy = self.det.indexes_xy(self.evt)
        self.det_cen = det_cen
        self.ROI = False
 
    def read_background(self):
        '''read in '''
	base_path = '/reg/d/psdm/cxi/cxim2716/res/carolin/'
	file_name = 'hits_r' + str(self.run.run()).zfill(4) + '.h5'
	file_path = base_path + file_name
	if os.path.isfile(file_path) is False:
            sys.stderr.write('Hit file does not exist for run %s.\n' % str(self.run.run()))
	    return

    def events(self, num=-1, start=0, roi=False):
        '''
        Get events from the given run for this class.
        Arguments:
            num = Number of frames to extract, set <=0 for all events (default=-1)
            start = Index of first frame (default=0)
            roi = Extract inner asics only (default=False)
        Returns:
            None. Sets self.frames as the given set of frames.
        '''
        if num <= 0:
            end = len(self.times)
        else:
            end = min(len(self.times), start+num)
        
        self.ROI = roi
        if roi:
            self.frames = np.zeros((end-start,) + (4,183,183))
        else:
            self.frames = np.zeros((end-start,) + self.ix.shape)
        
        for i in range(start, end):
            if roi:
                self.frames[i-start] = self.get_roi(self.det.calib(self.run.event(self.times[i])))
            else:
                self.frames[i-start] = self.det.calib(self.run.event(self.times[i]))
            sys.stderr.write('\r%d/%d' % (i-start+1,end-start))
        sys.stderr.write('\n')

    def get_frameslist(self, hit_list, roi=False):
        '''
        Get all frames in hitlist for the given run. Hits are selected by user .
        Arguments:
            hitlist = array of frames indices to be extracted
	    roi = Extract inner asics only (default=false)
        Returns:
            None. Sets self.frames as the given set of hit frames.
        '''
	self.hit_indices = hit_list
	n_hits = len(self.hit_indices)
	if n_hits == 0:
            sys.stderr.write('Hitlist contains no frames.\n')
	    return	
        self.ROI = roi
        if roi:
            self.frames = np.zeros((n_hits,) + (4,183,183))
        else:
            self.frames = np.zeros((n_hits,) + self.ix.shape)

        frame_iter=0
        for i in self.hit_indices:
            if roi:
                self.frames[frame_iter] = self.get_roi(self.det.calib(self.run.event(self.times[i])))
            else:
                self.frames[frame_iter] = self.det.calib(self.run.event(self.times[i]))
	    frame_iter = frame_iter + 1
            sys.stderr.write('\r%d/%d' % (frame_iter,n_hits))
        sys.stderr.write('\n')

    def get_hitframes(self, roi=False):
        '''
        Get all hits for the given run for this class. Hits are selected by partitioning the histogram of frame intensity sums.
        Arguments:
	    roi = Extract inner asics only (default=false)
        Returns:
            None. Sets self.frames as the given set of hit frames.
        '''
	base_path = '/reg/d/psdm/cxi/cxim2716/res/carolin/'
	file_name = 'hits_r' + str(self.run.run()).zfill(4) + '.h5'
	file_path = base_path + file_name
	if os.path.isfile(file_path) is False:
            sys.stderr.write('Hit file does not exist for run %s.\n' % str(self.run.run()))
	    return
	hits = h5.File(file_path, 'r')
	self.hit_indices = np.array(hits['hits/indices'])
	n_hits = len(self.hit_indices)
	if n_hits == 0:
            sys.stderr.write('Hit file for run %s contains no hits.\n' % str(self.run.run()))
	    return	
        self.ROI = roi
        if roi:
            self.frames = np.zeros((n_hits,) + (4,183,183))
        else:
            self.frames = np.zeros((n_hits,) + self.ix.shape)

        frame_iter=0
        for i in self.hit_indices:
            if roi:
                self.frames[frame_iter] = self.get_roi(self.det.calib(self.run.event(self.times[i])))
            else:
                self.frames[frame_iter] = self.det.calib(self.run.event(self.times[i]))
	    frame_iter = frame_iter + 1
            sys.stderr.write('\r%d/%d' % (frame_iter,n_hits))
        sys.stderr.write('\n')

    def get_roi(self, array):
        '''
        Get inner four ASICs of the given PSANA-shaped array.
        The second and fourth ASICS are rotated by 90 degrees to make the 
        orientations consistent.
        Arguments:
            array = PSANA shaped CSPAD array with shape (32,185,388)
        Returns:
            roi = Square regions of the inner 4 ASICs. Shape (4,183,183)
        '''
        if array.shape != self.ix.shape:
            sys.stderr.write('Array has the wrong shape to get roi\n')
	    return
        roi = array[1::8,1:184,1:184]
        roi[1] = np.rot90(roi[1])
        roi[3] = np.rot90(roi[3])
        return roi

    def do_fft(self, array):
        '''
        Average Fourier transform intensities of the inner 4 ASICs of a PSANA
        shaped array.
        Calls get_roi() internally to get the ROI.
        Arguments:
            array = PSANA shaped CSPAD array with shape (32,185,388)
                    or (4,183,183) shaped ROI array
        Returns:
            fft = Average FT intensities. Shape (183,183)
        '''
        if array.shape == self.ix.shape:
            return np.array([np.fft.fftshift(np.absolute(np.fft.fftn(i)*np.conjugate(np.fft.fftn(i)))) for i in self.get_roi(array)]).sum(axis=0)
        elif array.shape == (4,183,183):
            return np.array([np.fft.fftshift(np.absolute(np.fft.fftn(i)*np.conjugate(np.fft.fftn(i)))) for i in array]).sum(axis=0)

    def sum_roi(self, array):
        '''
        Sums ASIC intensities of the ROI of CSPAD array.
        Arguments:
            array = PSANA shaped CSPAD array with shape (32,185,388)
                    or (4,183,183) shaped ROI array
        Returns:
            fft = Average FT intensities. Shape (183,183)
        '''
        if array.shape == self.ix.shape:
            return np.array([i for i in self.get_roi(array)]).sum(axis=0)
        elif array.shape == (4,183,183):
            return np.array([i for i in array]).sum(axis=0)


    def assemble_frame(self, array):
        '''
        Assemble PSANA shaped frame into a 2D array.
        Arguments:
            array = PSANA shaped array
        Returns:
            assem = Assembled 2D array
        '''
        if self.ROI:
            sys.stderr.write('Event frame only contains ROI and cannot be assembled (call events() with roi=False) \n')
        if array.shape != self.ix.shape:
            sys.stderr.write('Array has the wrong shape to be assembled\n')
        assem = np.zeros((self.ix.max()+1, self.iy.max()+1))        
        assem[self.ix, self.iy] = array
        return assem

    def assemble_frame_w_geom(self, array, trim_asics):
        '''
        Assemble the psana array using 
        '''
        #print 'applying detector geometry mapping'
	base_path = '/reg/neh/home/dwojtas/cxils2616/scratch/dwojtas/'
	file_name = 'cxi_nano_metrol_v2_cxix29016_v0.h5'
	file_path = base_path + file_name
	if os.path.isfile(file_path) is False:
            sys.stderr.write('Geometry file %s does not exist.\n' % file_path)
	    return
	geom = h5.File(file_path, 'r')
	x_location = np.array(geom['/x'])
	y_location = np.array(geom['/y'])
        r_location = np.sqrt(x_location*x_location + y_location*y_location)
        px_size = 0.00011
        x_pix = x_location/px_size
        x_min = np.floor(x_pix.min())
        x_pix = x_pix + np.abs(x_min)
        x_max = np.ceil(x_pix.max())
        y_pix = y_location/px_size
        y_min = np.floor(y_pix.min())
        y_pix = y_pix + abs(y_min)
        y_max = np.ceil(y_pix.max())

        #print 'Centre? :(', x_max-np.abs(x_min), ',', y_max-np.abs(y_min), ')' 
        
        assem = np.zeros((x_max.astype(int),y_max.astype(int)));
        #roi = array[1::8,1:184,1:184]
        #array = event['det'].calib(event['evt'])
        cspad_np_og = array.reshape((4, 8, 185, 388))
        cspad_ij = np.zeros((1480,1552), dtype=cspad_np_og.dtype)
        
        for i in range(cspad_np_og.shape[0]):
            if trim_asics == True:
                for j in range(cspad_np_og.shape[1]):
                    cspad_np_og[i, j, 0, :] = 0
                    cspad_np_og[i, j, 184, :] = 0
                    cspad_np_og[i, j, :, 0] = 0
                    cspad_np_og[i, j, :, 387] = 0
                    cspad_np_og[i, j, :, 193] = 0
                    cspad_np_og[i, j, :, 194] = 0
                    
            #print cspad_np_og.shape[0]
            cspad_ij[:, i * cspad_np_og.shape[3]: (i+1) * cspad_np_og.shape[3]] = cspad_np_og[i].reshape((cspad_np_og.shape[1] * cspad_np_og.shape[2], cspad_np_og.shape[3]))            

        #self.det_array = cspad_ij
        #cspad_ij = np.transpose( cspad_ij)
        #x_pix = np.transpose( x_pix)
        #y_pix = np.transpose( y_pix)
        assem[ x_pix.astype(int), y_pix.astype(int)] = cspad_ij
        
        r_map = np.zeros((x_max.astype(int),y_max.astype(int)))
        r_map[x_pix.astype(int), y_pix.astype(int)] = r_location

        #print np.min(r_location)
        #print np.min(r_map)
        #origin_i =  np.where(r_map == np.min(r_location))
        #print origin_i
        #print r_map[origin_i[0], origin_i[1]]
        #print r_map[869, 924]

        #x_center = x_pix[r_map == np.min(r_map)]
        #y_center = y_pix[r_map == np.min(r_map)]
        
        #print np.where(r_map == np.min(r_map))
        #x_center_i = x_location == x_center
        #y_center_i = y_location == y_center
        ##print "Det cen: [", x_center, y_center, "]"
        #print "Det cen i: [", x_center_i, y_center_i, "]"
        r_map[self.det_cen[0], self.det_cen[1]] = 2 * np.max(r_map)
        #pq.show(r_map)
        
        #wait = input("PRESS ENTER TO FINISH..")
        #print "DONE"

        #x_pix = x_pix.flatten()
        #y_pix = y_pix.flatten()
        #array = array.flatten()

        #for ii in range(0, len(x_pix)):
        #    det_array[np.round(y_pix[ii]),np.round(x_pix[ii])]=array[ii]

        #print det_array.shape
        return assem
        
