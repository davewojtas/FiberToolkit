import get_events as ge
import pyqtgraph as pq
import numpy as np
import sys

e=ge.Events(run=131)
frameslist = np.array([54, 159, 1120, 1227, 1228, 3868])
e.get_frameslist(frameslist, False)

n_frames = len(e.frames)
for i in range(0,n_frames):
    pq.show(e.assemble_frame(e.frames[i]))

wait = input("PRESS ENTER TO CONTINUE..")
    
