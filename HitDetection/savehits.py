#!/usr/bin/env python

import pylab as P
import h5py
import sys
import os
import psana
import numpy as np
from mpi4py import MPI
import argparse

def update_h5file(f, key, data):
    if key in f: 
        del f[key]
    f[key] = data

parser = argparse.ArgumentParser(description='Create powder sum of all frames')
parser.add_argument('run', help='Choose run number', type=int)
#parser.add_argument('-t', '--threshold', help='Lit pixel threshold', type=float)
parser.add_argument('-o', '--output', help='Set output file')
args = parser.parse_args()

'''
if args.threshold is None:
    print 'Missing -t option. Need to specify lit pixel threshold'
    sys.exit(1)
'''

if args.output is None:
	hits_fname = 'r%.4d.h5' % args.run
else:
	hits_fname = args.output

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
num_proc = comm.Get_size()

ds = psana.DataSource('exp=cxils2616:run=%d:idx' % args.run)
run = ds.runs().next()
times = np.array(run.times())
evt = run.event(run.times()[0])
det = psana.Detector('DsaCsPad')
assem_shape = det.image(evt).shape
frame_shape = det.calib(evt).shape
ix, iy = det.indexes_xy(evt)
cx = det.coords_x(evt)
cy = det.coords_y(evt)
cz = det.coords_z(evt)
rad = np.hypot(cx, cy) / det.pixel_size(evt)


with h5py.File('r%.4d.h5'%args.run, 'r') as f:
    litpix = f['scores/litpix'][:] # Total pixels with at least 1 photon
    mcounts = f['scores/mask_counts'][:] # Photons inside mask
    total_counts = f['scores/total_counts'][:] # Total photons/frame

# Hit logic ==============================
ratios = np.zeros(len(mcounts))
ratios[mcounts[:,0]>0] = mcounts[mcounts[:,0]>0, 1].astype('f8') / mcounts[mcounts[:,0]>0, 0]
#hit_indices = np.where(litpix>args.threshold)[0]
hit_indices = np.where((litpix>400000) & (ratios>2.))[0]
# End hit logic ==========================

num_hits = hit_indices.shape[0]
times_p = times[hit_indices][rank::num_proc]
print 'Rank %.3d: %d/%d' % (rank, len(times_p), num_hits)

#f = h5py.File(hits_fname, 'a', driver='mpio', comm=comm)
f = h5py.File(hits_fname, 'a')
update_h5file(f, 'hits/indices', hit_indices)
if 'hits/data' in f: del f['hits/data']
if 'hits/assem' in f: del f['hits/assem']
dset_data = f.create_dataset('hits/data', (num_hits,)+frame_shape, chunks=(1,)+frame_shape, dtype='i4')
dset_assem = f.create_dataset('hits/assem', (num_hits,)+assem_shape, chunks=(1,)+assem_shape, dtype='i4')

for i, t in enumerate(times[hit_indices]):
    if i % num_proc != rank:
        continue
    phot = det.photons(run.event(t), adu_per_photon=23)
    if phot is None:
        nbad_p[0] += 1
        continue
    dset_data[i] = phot
    dset_assem[i] = det.image(run.event(t), nda_in=phot)
    
    if rank == 0:
        sys.stderr.write('\r%d/%d' % (i, num_hits))

if rank == 0:
    sys.stderr.write('\n')
f.close()

MPI.Finalize()
