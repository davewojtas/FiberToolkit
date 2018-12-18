import get_events as ge
import merge_events as me
import frame_orientate as fo
import pyqtgraph as pq
import numpy as np
import sys
import matplotlib.pyplot as p

substrate_tilt = 0
wavelength = 0
det_dist = 0
print "Getting events"
e=ge.Events(run=57)

#e.get_hitframes(False) # get frames identified as potential hits
#hitlist = np.array([54, 622, 3868, 3873 ]) #angles 22-27
#hitlist = np.array([54, 159, 1120, 1227, 1228, 3868]) #159
#hitlist = np.array([1120,1227,3785,3873,3874]) #angles 120-135
hitlist = np.array([59]) #angles 120-135
print "Getting frames from event.. "
e.get_frameslist(hitlist, False)
#e.events(4351, 0, roi=False) # get all frames in run (4351 frames for run 131)

r1 = [2,2,2,2,3,3,3,3,4,4,4,4]
r2 = [0.005, 0.01, 0.015, 0.02, 0.005, 0.01, 0.015, 0.02, 0.005, 0.01, 0.015, 0.02]

r1 = [3]
r2 = [0.02]

#sys.stdout = open("data_output.txt", "w")

max_intensity = 80000
for i in range(0, 1):
    print "------------------"
    print "----", r1[i], "  ", r2[i], "----"
    m=me.merge_events(wavelength, det_dist, e.assemble_frame_w_geom(e.frames[0]).shape)
    max_frames = len(e.frames)
    for h in range(0,max_frames):
        current_frame = np.copy(e.frames[h]) # copy for performing rotation and fft
        #pq.show(e.assemble_frame(current_frame))
        #wait = input("PRESS ENTER TO FINISH..")       
        mean_px_intensity = np.mean(current_frame)
        if mean_px_intensity < max_intensity:
            hit_ft = e.do_fft(current_frame)
            f = fo.frame_orientate(hit_ft)
            f.plot_graphs = 0
            f.determine_orientation(kern_win=5, peak_to_mean_ratio=r1[i], ann_ir=30, ann_or=70)
            #f.calc_autocorrelation() 
            f.calc_autocorrelation_1d(f.orientation, ray_width=5, peak_threshold=r2[i])
            #pq.show(f.r_map)
            #wait = input("PRESS ENTER TO FINISH..")
            #f.radial_profile_peak_detect(f.orientation, 1, 0.1, r2[i])
            #if f.is_hit_s2 is True and f.n_orientations >= 1: 
            if f.is_hit_s1 is True :# and f.is_hit_s2 is True :
                print h, f.orientation[0], f.n_orientations, f.mean_ang_separation, mean_px_intensity, 1, 1
                #print h
                #f.orientation = 0
                #m.merge_RZ( e.assemble_frame_w_geom(e.frames[h]), -f.orientation, 0, False)
                m.cummul_sum_det(-f.orientation, e.assemble_frame_w_geom(e.frames[h]) )
  

#m.cumm_sum[ 875,:] = np.max(m.cumm_sum) / 2
#pq.show(m.cumm_sum)

p.figure()
maxval = np.max(m.cumm_sum)
thresh = 0.05 * mean_px_intensity / 200
if thresh < 0.05:
    thresh = 0.05
m.cumm_sum[m.cumm_sum > thresh * maxval] = thresh * maxval 
m.cumm_sum[875,887] = thresh *maxval
p.imshow(m.cumm_sum, cmap='gnuplot', vmin=0, origin="lower")
p.show()

#p = pq.show(m.merged_data)
#pq.mkQApp()
#w = pq.GraphicsLayoutWidget()
#w.show()
#vb = w.addViewBox()
#img = pq.show(m.merged_data)# pq.ImageItem(np.random.normal(size=(100,100))) #m.merged_data)
#vb.addItem(img)
#def mouseMoved(pos):
#    print "Image position:", img.mapFromScene(pos)

#w.scene().sigMouseMoved.connect(mouseMoved)

#pq.SignalProxy(p.scene().sigMouseMoved, rateLimit=60, slot=callback)
wait = input("PRESS ENTER TO FINISH..")
#print "DONE"
