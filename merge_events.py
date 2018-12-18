import numpy as np
import scipy as sp
from scipy import ndimage
import sys
import pyqtgraph as pq

class merge_events():
    def __init__(self, wavelength, det_dist, data_size,det_cen):
        self.wl = wavelength
        self.detd = det_dist #85mm
        self.q_size = 1701
        n_pix = 2296960
        self.sum_aligned_data = np.zeros(data_size)
        self.det_cen = det_cen
        #self.det_cen = [np.floor(data_size[0]/2), np.floor(data_size[1]/2)]
        self.merged_data = np.zeros((self.q_size,self.q_size))
        self.cumm_sum = np.zeros(data_size)
        self.frame_count = 0

    def merge_RZ(self, array, phi=0, beta=0, mirror=True):
        '''
        Map the frame 'array' into reciprocal space and merge into class 'merged_data'
        Arguments:
            phi = in-plane angle of fiber relative to detector y-axis
            array = hit frame data to be merged 
        Returns:
        '''
        subtract_rad_ave = False
        # calculate scale factor using bg mask
        sp = np.sin(phi * np.pi / 180)
        cp = np.cos(phi * np.pi / 180)
        sb = np.sin(beta * np.pi / 180)
        cb = np.cos(beta * np.pi / 180)
        px_size = 0.00011
        DETD = 913.27
        MSCALE = DETD
        size = 1701
        px_x, px_y = np.indices(array.shape)
        px_x = px_x - self.det_cen[1]
        px_y = px_y - self.det_cen[0]

        #subtract the radial average intensity distribution
        if subtract_rad_ave is True:
            r_dist = np.zeros(array.shape)
            px_r = np.hypot(px_x, px_y).astype(int)
            tbin = np.bincount(px_r.ravel(), array.ravel())
            nr = np.bincount(px_r.ravel())
            rad_ave = tbin / nr
            for i in range(0, px_r.max()+1):
                r_dist[px_r==i] = rad_ave[i]
            array = array - r_dist
            neg_indices = array < 0
            array[neg_indices] = 0

        #calculate polarization factor
        px_pol = (1. - px_y*px_y / (px_x*px_x + px_y*px_y + DETD*DETD)) 
        #print self.det_cen[1]
        #self.merged_data = px_y
        #return
        # for each pixel, calculate (R, Z) and interpolate
        #for i in range(0, n_pix):
        px_xr = px_x*cp - px_y*sp
	px_yr = px_x*sp + px_y*cp 
        # calculate 3D position
        factor = np.sqrt(px_xr*px_xr + px_yr*px_yr + DETD*DETD)         
	qx = MSCALE * (px_xr / factor)
	qy = MSCALE * (px_yr / factor)
	qz = MSCALE * (DETD / factor - 1) 
        # Rotate Ewald sphere by beta
	rqy = cb*qy - sb*qz
	rqz = sb*qy + cb*qz
	# Calculate (R, Z) coordinates
	r = np.sqrt(qx*qx + rqz*rqz) + size/2
	z = rqy + size/2
        #self.merged_data = z
        #return
        #print "rshape", r.shape
        #print "zshape", z.shape
        #r = r.flatten()
        #z = z.flatten()
        #array = array.flatten()
        #self.merged_data = np.zeros(self.q_size * self.q_size)
        
        ir = np.floor(r)
        ir = ir.astype(int)
        iz = np.floor(z)
        iz = iz.astype(int)
        #self.merged_data = z
        #return
        #print ir.max()
        #print iz.max()
        #print ir.min()
        #print iz.min()

        #ob_mask = np.ones(array.shape)
        #ob_mask[ir > size-2] = 0
        #ob_mask[iz > size-2] = 0
        #ob_mask[iz < 0] = 0
        #ob_mask = ob_mask.astype(int)

        #print ir.shape
        #inBounds = ir < size-1
        #print ir[outOfBounds]
        #r_sub = r[inBounds]
        #ir_sub = ir[inBounds]# =0
        #print ir_sub.shape
        #z_sub = z[inBounds]
        #iz_sub = iz[inBounds]# =0
        #print iz_sub.shape
        #array_sub = array[inBounds]# =0
        #print array_sub.shape
        '''inBounds = iz_sub >= 0
        r_sub = r_sub[inBounds]
        ir_sub = ir_sub[inBounds]# =0
        z_sub = z_sub[inBounds]
        iz_sub = iz_sub[inBounds]# =0
        array_sub = array_sub[inBounds]# =0
        inBounds = iz_sub < size-1
        r_sub = r_sub[inBounds]
        ir_sub = ir_sub[inBounds]# =0
        z_sub = z_sub[inBounds]
        iz_sub = iz_sub[inBounds]# =0
        array_sub = array_sub[inBounds]# =0'''
        #ir[outOfBounds] =0
        #iz[outOfBounds] =0
        #array[outOfBounds] =0
        
	# fixup below
        #val = data[t] * pix_pol[t] ;
			
	#array[r > size-2] = 0 
        #array[z < 0 || z > size-2)
	#	continue ;
			
	#ir = np.floor(r)
        #ir = r#.astype(int)
	#iz = np.floor(z)
        #iz = z#.astype(int)
        fr = r - ir 
	fz = z - iz
	cr = 1 - fr
	cz = 1 - fz			
        #self.merged_data[ir*ob_mask, iz*ob_mask] += cz*cr*array*ob_mask
        #self.merged_data[(ir+1)*ob_mask, iz*ob_mask] += cz*fr*array*ob_mask
	#self.merged_data[ir*ob_mask, (iz+1)*ob_mask] += fz*cr*array*ob_mask
	#self.merged_data[(ir+1)*ob_mask, (iz+1)*ob_mask] += fz*fr*array*ob_mask

        array = array * px_pol
        self.merged_data[ir, iz] += cz*cr*array
        self.merged_data[(ir+1), iz] += cz*fr*array
	self.merged_data[ir, (iz+1)] += fz*cr*array
	self.merged_data[(ir+1), (iz+1)] += fz*fr*array
	# Merge `negative' r
	r = size - r
        ir = np.floor(r)
        ir = ir.astype(int)
        #outOfBounds = ir > size-2
        #print r[outOfBounds]
        #ir[outOfBounds] =0
        #z[outOfBounds] =0
        #array[outOfBounds] =0
	#ir = np.floor(r)
        #ir = r.astype(int)
	fr = r - ir 
	cr = 1 - fr			
	self.merged_data[ir, iz] += cz*cr*array
	self.merged_data[(ir+1), iz] += cz*fr*array
	self.merged_data[ir, (iz+1)] += fz*cr*array
	self.merged_data[(ir+1), (iz+1)] += fz*fr*array

    def rotateImage(self, img, angle):
        pivot = self.det_cen
        #print 'Det. Centre = (', pivot[0], pivot[1], ")"
        #img_copy = np.copy(img)
        #img[pivot[0],pivot[1]] = np.max(img_copy)
        #pq.show(img_copy)
        padX = [img.shape[0] - pivot[0], pivot[0]]
        padY = [img.shape[1] - pivot[1], pivot[1]]
        #padX = [img.shape[1] - pivot[0], pivot[0]]    
        #padY = [img.shape[0] - pivot[1], pivot[1]]
        #imgP = np.pad(img, [padY, padX], 'constant')
        imgP = np.pad(img, [padX, padY], 'constant')
        imgR = ndimage.rotate(imgP, angle, reshape=False)
        croppedImg = imgR[padX[0] : -padX[1], padY[0] : -padY[1]]
        #croppedImg = imgR[padY[0] : -padY[1], padX[0] : -padX[1]]
        croppedImg[pivot[0],pivot[1]] = np.max(croppedImg)
        #pq.show(croppedImg)
        
        #wait = input("PRESS ENTER TO FINISH..")
        #print "DONE"
        return croppedImg

    def cummul_sum_det(self, phi, array):
        '''
        Sum rotated frames in detector space. Orientation angle is assumed to be purely phi i.e. beta = 0
        Arguments:
            phi = in-plane angle of fiber relative to detector y-axis
            array = hit frame data to be merged d
        Returns:
        '''
        #rotate input array
        #array_rot = np.rot(array, phi)
        #add to cummulative sum
        # calculate scale factor using bg mask
        #if self.frame_count == 0:
            #phi = phi - 180
        #sp.ndimage.interpolation.rotate(array, phi, reshape=False, output=array, order=3, mode='constant', cval=0.0, prefilter=True)   

        #array[self.det_cen[0], self.det_cen[1]] = np.max(array)
        #array_copy = np.copy(array)
        #pq.show(array_copy)
        #array = self.rotateImage(array, phi)
        #array[self.det_cen[0], self.det_cen[1]] = np.max(array)
        #pq.show(array)   
        
        #print array_rot.shape
        #if self.frame_count == 0:
        #    array = (array - array.min())/(array.max() - array.min()) 
        #    array = array + 1        
        #if self.frame_count == 1:
        #    array = (array - array.min())/(array.max() - array.min()) 
        #    array = 1 - array            
        self.frame_count = self.frame_count + 1
        self.cumm_sum = self.cumm_sum + array
