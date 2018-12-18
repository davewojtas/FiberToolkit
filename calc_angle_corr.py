import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

# schema {step, time, tstep,y,x,z,dy,dx,dz,shots}
holeList = np.loadtxt('frames_window_coords.csv',delimiter=',', usecols=(0,1,2,3,4,5))
# schema {frame, orientation,nOrientations, ang_separation,px_intensity
frameList = np.loadtxt('framesList_allStage1.csv',delimiter=',', usecols=(0,1,2,3,4))

frameIndices = frameList[:,0]
frameOrientations = frameList[:,1]

print np.shape(holeList)
res = 3 # number dp
x = np.round(holeList[:,4],res)
y = np.round(holeList[:,3],res)
delta_x = 0.06/1
x= (x/delta_x).astype(int)
y= (y/delta_x).astype(int)
print np.min(x), np.min(y)
print np.max(x), np.max(y)
x = (x - np.min(x))
y = (y - np.min(y))

#print np.min(x), np.min(y)
#print np.max(x), np.max(y)

x_mesh, y_mesh = np.mgrid[slice(0, np.max(y)+1, 1),
                          slice(0, np.max(x)+1, 1)]

frame_coord_index =  np.round(frameIndices * 4206 / 4351).astype(int)

z_mesh = 0*y_mesh + 180
z_mesh[y.astype(int),x.astype(int)] = 0
z_mesh[y[frame_coord_index].astype(int),x[frame_coord_index].astype(int)] = frameOrientations

print np.max(z_mesh)

#for i in range(0, len(x)):
#    z_mesh[y[i],x[i]] = 90

plt.figure()
plt.axis('off')
ax = plt.gca()
im = ax.pcolor(x_mesh, y_mesh, z_mesh, cmap='gnuplot2', vmin = 0, vmax=180)
plt.xlim([0,np.max(x_mesh)])
plt.ylim([0,np.max(y_mesh)])
plt.xlabel('y')
plt.ylabel('x')
plt.axes().set_aspect('equal')
divider = make_axes_locatable(ax)
cax = divider.append_axes("top", size="2.5%", pad=0.2)
#plt.colorbar(im, cax=cax)
#plt.colorbar(im)
cbar = plt.colorbar(im, cax=cax, ticks=[0, 45, 90, 135, 180], orientation='horizontal')
cbar.ax.tick_params(axis='x',direction='in', labeltop='on', labelbottom='off')  # horizontal colorbar
plt.show()

#autoCorr = np.correlate(frameOrientations, frameOrientations, mode='full')
#autoCorr = autoCorr[autoCorr.size/2:]

#x_phys = np.round(holeList[:,4],res)
#y_phys = np.round(holeList[:,3],res)
#r_binSize = 0.1
#r_max = np.sqrt( (np.max(x) - np.min(x))**2 + (np.max(y) - np.min(y))**2)
#r_arraySize = (np.ceil(r_max / r_binSize) + 1).astype(int)
#r_array = np.zeros(r_arraySize)
#r_arrayCount = np.zeros(r_arraySize)

#for i in range(0, len(frame_orientations)-1):
#    xi = x_phys[frame_coord_index[i]]
#    yi = y_phys[frame_coord_index[i]]
#    phii = frame_orientations[i]
#    for j in range(i+1, len(frame_orientations)):
#        xj = x_phys[frame_coord_index[j]]
#        yj = y_phys[frame_coord_index[j]]
#        phij = frame_orientations[j]
#        r = np.sqrt((xi - xj)**2 + ( yi - yj)**2)
#        delta_phi = np.abs(phii - phij)
#        r_array[(r / r_binSize).astype(int)] += delta_phi
#        r_arrayCount +=1

#non_zero_bins = r_array[
        
