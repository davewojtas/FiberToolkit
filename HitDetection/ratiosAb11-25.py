#!/usr/bin/env python

import pylab as P
import h5py
import sys
import os
import psana
import numpy as np
from mpi4py import MPI
import argparse

parser = argparse.ArgumentParser(description='Create powder sum of all frames')
parser.add_argument('run', help='Choose run number', type=int)
parser.add_argument('-o', '--output', help='Set output file')
args = parser.parse_args()

if args.output is None:
	hits_fname = 'r%.4d.h5' % args.run
else:
	hits_fname = args.output

mask_fname = '/home/cseuring/LS26-Hitfinding/scripts/mask.h5'

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
num_proc = comm.Get_size()

ds = psana.DataSource('exp=cxils2616:run=%d:idx' % int(sys.argv[1]))
run = ds.runs().next()
times = run.times()
evt = run.event(run.times()[0])
det = psana.Detector('DsaCsPad')
assem_shape = det.image(evt).shape
frame_shape = det.calib(evt).shape
ix, iy = det.indexes_xy(evt)
cx = det.coords_x(evt)
cy = det.coords_y(evt)
cz = det.coords_z(evt)
rad = np.hypot(cx, cy) / det.pixel_size(evt)

times_p = times[rank::num_proc]
powder_p = np.zeros(frame_shape)
print 'Rank %.3d: %.6d/%.6d' % (rank, len(times_p), len(times))

# load user mask from .h5 file
f = h5py.File(mask_fname, 'r')
mask1 = f['amyloidmask'][:].astype(np.bool)
#mask2 = f['amyloidmaskbg1'][:].astype(np.bool)
mask2 = f['amyloidmaskbg2'][:].astype(np.bool)
f.close()

print mask1.shape, mask1.dtype, mask2.shape, mask2.dtype

mask_counts = np.zeros((len(times), 2), dtype='i4')
total_counts = np.zeros((len(times),), dtype='i4')
litpix = np.zeros((len(times),), dtype='i4')
nbad_p = np.zeros((1,), dtype='i4') 

for i, t in enumerate(times_p):
    #if i > 9:
    #    break
    phot = det.photons(run.event(t), adu_per_photon=23)
    if phot is None:
        nbad_p[0] += 1
        continue
    powder_p += phot 
    mask_counts[num_proc*i+rank] = [phot[mask1].sum(), phot[mask2].sum()]
    total_counts[num_proc*i+rank] = phot.sum()
    litpix[num_proc*i+rank] = (phot>0).sum()
    if rank == 0:
        sys.stderr.write('\r%d/%d' % (i, len(times_p)))

def update_h5file(f, key, data):
    if key in f: 
        del f[key]
    f[key] = data

if rank == 0:
    sys.stderr.write('\n')
    red = np.empty(frame_shape)
    nbad = np.empty_like(nbad_p)
    mred = np.empty_like(mask_counts)
    tred = np.empty_like(total_counts)
    lred = np.empty_like(litpix)
    comm.Reduce([powder_p, MPI.DOUBLE], [red, MPI.DOUBLE], op=MPI.SUM, root=0)
    comm.Reduce([nbad_p, MPI.INT], [nbad, MPI.INT], op=MPI.SUM, root=0)
    comm.Reduce([mask_counts, MPI.INT], [mred, MPI.INT], op=MPI.SUM, root=0)
    comm.Reduce([total_counts, MPI.INT], [tred, MPI.INT], op=MPI.SUM, root=0)
    comm.Reduce([litpix, MPI.INT], [lred, MPI.INT], op=MPI.SUM, root=0)
    print 'Found %d frames without data' % nbad
    red /= len(times) - nbad
    assem_red = np.zeros(assem_shape)
    assem_red[ix, iy] = red
    
    f = h5py.File(hits_fname, 'a')
    update_h5file(f, 'geom/cx', cx)
    update_h5file(f, 'geom/cy', cy)
    update_h5file(f, 'geom/cz', cz)
    update_h5file(f, 'geom/ix', ix)
    update_h5file(f, 'geom/iy', iy)
    
    update_h5file(f, 'event/timestamps', np.array([t.time() for t in times], dtype='u8'))
    update_h5file(f, 'event/fiducial', np.array([t.fiducial() for t in times], dtype='u2'))
    
    update_h5file(f, 'powder/calib', red)
    update_h5file(f, 'powder/image', assem_red)
    update_h5file(f, 'powder/nevents', len(times) - nbad)
    
    update_h5file(f, 'scores/mask_counts', mred)
    update_h5file(f, 'scores/masks', np.array([mask1, mask2], dtype='u1'))
    update_h5file(f, 'scores/total_counts', tred)
    update_h5file(f, 'scores/litpix', lred)
    f.close()
else:
    comm.Reduce([powder_p, MPI.DOUBLE], None, op=MPI.SUM, root=0)
    comm.Reduce([nbad_p, MPI.INT], None, op=MPI.SUM, root=0)
    comm.Reduce([mask_counts, MPI.INT], None, op=MPI.SUM, root=0)
    comm.Reduce([total_counts, MPI.INT], None, op=MPI.SUM, root=0)
    comm.Reduce([litpix, MPI.INT], None, op=MPI.SUM, root=0)

MPI.Finalize()
