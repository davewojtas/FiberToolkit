import get_events as ge
import merge_events as me
import frame_orientate as fo
import pyqtgraph as pq
import numpy as np
import sys
import matplotlib.pyplot as p

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

class Layerline_Profile:
    def __init__(self, layerline, data):
        self.layerline = int(layerline)
        self.data = data
        self.data_incl = 0*data + 1 # no exclusions 
        self.data_symres = np.zeros(len(data)/2 + 1)
    def Show(self):
        x = np.arange(-len(Data)/2, len(data))
        plt.figure()
        plt.plot(x, self.data)
        

verbose_mode = True
substrate_tilt = 0
run_numbers = [242] #227, 139
n_frames = 10000
start_frame = 0
shooting_silicon = 0
silicon_mean_offset = 15
det_cen = [870,884] #determined from oleksandr's geom
mode = 'tmv' # adeno, tmv, nanotube
if mode=='nanotube' or mode=='abeta':
    min_mean = 30 + shooting_silicon * silicon_mean_offset
    max_mean = 1025 + shooting_silicon * silicon_mean_offset
    determine_orientation = False
    max_n_orientations = 200
    min_orientation_snr = 0.0
    detect_asymmetry = False
    detect_layerline = True
    merge_data = False
    #rotate_frames = False
elif mode=='tmv':
    min_mean = -1000000*25 + shooting_silicon * silicon_mean_offset
    max_mean = 4000000 + shooting_silicon * silicon_mean_offset #225
    determine_orientation = False
    max_n_orientations = 500#25 # 5 # 500
    min_orientation_snr = 0#3 # 4 # 0
    detect_asymmetry = False
    detect_layerline = False
    merge_data = False
    #rotate_frames = True


#get geom shape and beam characteristics
e=ge.Events(run=run_numbers[0],det_cen=det_cen)
beam_energy = e.epics.value('BEND:DMP1:400:BDES')
wavelength = e.epics.value('SIOC:SYS0:ML00:AO192')
det_dist = e.epics.value('CXI:DS2:MMS:06.RBV')
e.events(1, 0, roi=False)
m=me.merge_events(wavelength, det_dist, e.assemble_frame_w_geom(e.frames[0],False).shape, det_cen)

for r in range(0, len(run_numbers)):
    print "Getting events for run " + str(run_numbers[r])
    e=ge.Events(run=run_numbers[r],det_cen=det_cen)
    e.events(n_frames, start_frame, roi=True)

    mean_vals = np.zeros([len(e.frames), 2])
 
    for i in range(0, len(e.frames)):
        mean_vals[i,0] = i
        mean_vals[i,1] = np.mean(e.frames[i])

    #log_mean_vals = np.sqrt(mean_vals[:,1] - np.min(mean_vals[:,1]))

    #p.figure()
    #n, bins, patches = p.hist(log_mean_vals, 50, normed=0)
    #p.show()

    hits0 = mean_vals[mean_vals[:,1] > min_mean]
    hits1 = []
    for i in range(0, len(hits0)):
        current_frame = np.copy(e.frames[hits0[i,0]])
        #if np.max(current_frame) > 2000:
            #print np.max(current_frame)
            #continue
        if determine_orientation == True:
            hit_ft = e.do_fft(current_frame)
            f = fo.frame_orientate(hit_ft, det_cen)
            f.plot_graphs = 0
            f.determine_orientation(kern_win=7, peak_to_mean_ratio=min_orientation_snr, ann_ir=30, ann_or=70)    
            if f.n_orientations <= max_n_orientations and f.is_hit_s1:
                hit = Hit_Stat(run_numbers[r], shooting_silicon, hits0[i,0]+start_frame, f.orientation, f.orientation_snr, f.n_orientations, hits0[i,1])
                #hit.Print()
                hits1.append(hit) 
        else:
            hit = Hit_Stat(run_numbers[r], shooting_silicon, hits0[i,0]+start_frame, 0, 0, 0, hits0[i,1])
            #hit.Print()
            hits1.append(hit)
     
    
    
    #save stage 1 hits to .npy
    file_name = 'hits_run'
    if run_numbers[r] < 100:
        file_name = file_name + '0' + str(run_numbers[r])
    else:
        file_name = file_name + str(run_numbers[r])
    if shooting_silicon == 1:
        file_name = file_name + '.npy'
    else:
        file_name = file_name + '.npy'

    np.save(file_name, hits1)
    wait = input("PRESS ENTER TO FINISH..")
    print "DONE"

    print 'N. hits=', len(hits1)
    hit_frames = np.zeros((len(hits1)),dtype=int)

    for i in range(0,len(hits1)):
        hit_frames[i] = int(hits1[i].frame)

    e.get_frameslist(hit_frames, False) 

    if merge_data == True or detect_asymmetry == True or verbose_mode == True:

        for i in range(0,len(hits1)):
            
            if merge_data == True or detect_asymmetry == True:
                assembled_frame = e.assemble_frame_w_geom(e.frames[i], True)
                if determine_orientation == True:
                    assembled_frame = m.rotateImage(assembled_frame, -hits1[i].phi)    
                dims = np.shape(assembled_frame)
        
            if detect_asymmetry == True:
        
                iw = 50 # integration_width / 2 - 1
                ih = 4 # integration_height / 2 - 1
                ll_px = 22.5 #layerline_spacing_px
           
                off_centre_px = 2*(det_cen[0] - np.floor(dims[0]/2.0)) - 1
                lhs_start = np.max([0, off_centre_px])
                lhs_finish = det_cen[0]
                rhs_start = det_cen[0]
                rhs_finish = np.min([dims[0], dims[0] + off_centre_px])
                lhs_image = assembled_frame[lhs_start:lhs_finish, :]
                rhs_image = assembled_frame[rhs_start:rhs_finish, :]
                #rhs_image = np.flipud(rhs_image)
                residual_asym = np.copy(assembled_frame)
                residual_asym[lhs_start:lhs_finish, :] -= np.flipud(rhs_image)
                residual_asym[rhs_start:rhs_finish, :] -= np.flipud(lhs_image)

                ll_int = [-7, -6, -3, 3, 6, 7]
                ll_symres_aves = np.zeros(len(ll_int))
                #ll_int = [-3, 3]
                #ll_profs = [Layerline_Profile()

                fig = p.figure()
                max_val = np.max(assembled_frame)
                #asym_values = ''
                for l in range(0, len(ll_int)):
                    h_multiplier = np.abs(ll_int[l] / 3.)
                    ih_m = ih * h_multiplier.astype(int)
                    w_adjustor = 20 * (np.abs(ll_int[l] / 3.) - 1)
                    iw_m = iw + w_adjustor.astype(int)
                    prof_avg = np.zeros([iw_m*2 +1])
                    prof_incl = np.ones([iw_m*2 +1])
                    prof_res = np.zeros([iw_m*2+1])
                
                    for h in range(-ih_m, ih_m+1):
                        line_profile = assembled_frame[-iw_m+det_cen[0]:iw_m+det_cen[0]+1, det_cen[1] + ll_int[l]*ll_px + h] 
                        prof_avg += line_profile                    
                        prof_incl[line_profile <= 0] = 0

                    for w in range(0,iw_m + 1):       
                        #print w
                        sym_prod = prof_incl[iw_m - w] * prof_incl[iw_m + w]
                        prof_incl[iw_m - w] = sym_prod
                        prof_incl[iw_m + w] = sym_prod
                        prof_res[iw_m - w] = sym_prod * np.abs( prof_avg[iw_m - w] -  prof_avg[iw_m + w] )
                        prof_res[iw_m + w] = sym_prod * np.abs( prof_avg[iw_m - w] -  prof_avg[iw_m + w] )

                    x= np.arange(-iw_m+det_cen[0],iw_m+det_cen[0]+1)
                    x_res= np.arange(det_cen[0],iw_m+det_cen[0]+1)
                    lp = Layerline_Profile(l, prof_avg)
                    lp.data_incl = prof_incl
                    lp.data_symres = prof_res
                    sym_res = lp.data_symres[len(prof_res)/2:]
                    ll_symres_aves[l] = np.sum(sym_res)/len(sym_res)
                    #p.plot(x, lp.data, label='l='+str(ll_int[l]))
                    #p.plot(x_res, sym_res, label='l_res='+str(ll_int[l]) )            
                    #p.legend()
            
                hits1[i].sym_res_ave = np.max(ll_symres_aves)            
                hits1[i].sym_res_ll = 0#ll_int
                #p.show()
                #assembled_frame[det_cen[0],:] = max_val
            
                #need to save .png of integration regions

                #pq.show(residual_asym)

                #wait = input("PRESS ENTER TO FINISH..")
                #print "DONE"      

            if merge_data == True:
                m.cummul_sum_det(0, assembled_frame)
                #m.merge_RZ(rotated_image, 0, 0, True)
            if verbose_mode == True:                
                print hits1[i].frame, hits1[i].roi_px_mean, hits1[i].phi, '(',hits1[i].phi_peak_to_mean, hits1[i].n_phi_vals,')  ', hits1[i].sym_res_ave
       
    
np.save(file_name, hits1)
print 'Total ', m.frame_count, ' frames merged'

if merge_data == True:
    pq.show(m.cumm_sum/m.frame_count)
    #pq.show(m.merged_data)

#m.merged_data[m.merged_data < 0] = 0 
#print np.min(m.merged_data)
#pq.show(m.merged_data)


#p.figure()
#maxval = np.max(m.cumm_sum)
#thresh = 0.05 * 100 / 200
#if thresh < 0.05:
#    thresh = 0.05
#m.cumm_sum[m.cumm_sum > 0.2 * maxval] = 0.2 * maxval 
#m.cumm_sum[875,887] = thresh *maxval
#p.matshow(m.cumm_sum, cmap='gnuplot', vmin=0, origin="lower")
#p.show()

wait = input("PRESS ENTER TO FINISH..")
print "DONE"

#r1 = [2,2,2,2,3,3,3,3,4,4,4,4]
#r2 = [0.005, 0.01, 0.015, 0.02, 0.005, 0.01, 0.015, 0.02, 0.005, 0.01, 0.015, 0.02]

#r1 = [3]
#r2 = [0.02]

#sys.stdout = open("data_output.txt", "w")

#max_intensity = 80000
#for i in range(0, 1):
    #print "------------------"
    #print "----", r1[i], "  ", r2[i], "----"
    #m=me.merge_events(wavelength, det_dist, e.assemble_frame_w_geom(e.frames[0]).shape)
    #max_frames = len(e.frames)
    #for h in range(0,max_frames):
        #current_frame = np.copy(e.frames[h]) # copy for performing rotation and fft
        #pq.show(e.assemble_frame(current_frame))
        #wait = input("PRESS ENTER TO FINISH..")       
        #mean_px_intensity = np.mean(current_frame)
        #if mean_px_intensity < max_intensity:
            #hit_ft = e.do_fft(current_frame)
            #f = fo.frame_orientate(hit_ft)
            #f.plot_graphs = 0
            #f.determine_orientation(kern_win=5, peak_to_mean_ratio=r1[i], ann_ir=30, ann_or=70)
            #f.calc_autocorrelation() 
            #f.calc_autocorrelation_1d(f.orientation, ray_width=5, peak_threshold=r2[i])
            #pq.show(f.r_map)
            #wait = input("PRESS ENTER TO FINISH..")
            #f.radial_profile_peak_detect(f.orientation, 1, 0.1, r2[i])
            #if f.is_hit_s2 is True and f.n_orientations >= 1: 
            #if f.is_hit_s1 is True :# and f.is_hit_s2 is True :
                #print h, f.orientation[0], f.n_orientations, f.mean_ang_separation, mean_px_intensity, 1, 1
                #print h
                #f.orientation = 0
                #m.merge_RZ( e.assemble_frame_w_geom(e.frames[h]), -f.orientation, 0, False)
                #m.cummul_sum_det(-f.orientation, e.assemble_frame_w_geom(e.frames[h]) )
  

#m.cumm_sum[ 875,:] = np.max(m.cumm_sum) / 2
#pq.show(m.cumm_sum)

#p.figure()
#maxval = np.max(m.cumm_sum)
#thresh = 0.05 * mean_px_intensity / 200
#if thresh < 0.05:
#    thresh = 0.05
#m.cumm_sum[m.cumm_sum > thresh * maxval] = thresh * maxval 
#m.cumm_sum[875,887] = thresh *maxval
#p.imshow(m.cumm_sum, cmap='gnuplot', vmin=0, origin="lower")
#p.show()

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
#wait = input("PRESS ENTER TO FINISH..")
#print "DONE"
