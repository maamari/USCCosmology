'''
Generates data for ML algo
using pixell
cluster works
'''

import multiprocessing as mp
import healpy as hp
import numpy as np
import pandas as pd
from pixell import enmap, reproject, utils


voids = pd.read_csv('data/voids.csv', header=None)
voids.columns = ['Sample', 'ID', 'RAdeg', 'DEdeg', 'z', 'NGal', 'V', 'Reff',
'nmin', 'delmin', 'r', 'P', 'Dbound']

pl = hp.read_map("data/HFI_SkyMap_353-psb-field-IQU_2048_R3.00_full.fits")
pl[pl < -1e23] = 0

shape, wcs = enmap.fullsky_geometry(res=5.0*utils.arcmin, proj='car')
map_pix = reproject.enmap_from_healpix(pl, shape=shape, wcs=wcs)

signal = []
ras_temp = []
decs_temp = []

ras_fake = np.linspace(0, 360, 100)
decs_fake = np.linspace(-90, 90, 100)

N = len(ras_fake)
# N=10
stack_CMB_kSZ = 0
stack_CMB_kSZ_deproject = 0
c = 0


def worker(i):
    for j in range(N):
        stamp = reproject.thumbnails(map_pix, coords=np.deg2rad([decs_fake[i],
                                     ras_fake[j]]), r=3*utils.arcmin)

        if stamp is None:
            continue
        elif stamp[0][0][0] == 0.0:
            continue

        signal.append(sum(sum(stamp[0]))/9)
        ras_temp.append(ras_fake[j])
        decs_temp.append(decs_fake[i])
    return [signal, ras_temp, decs_temp]

# pool = mp.Pool(processes=8)
# for _ in tqdm(pool.imap_unordered(worker, range(N)), total=len(range(N))):
#     pass

with mp.Pool(processes = 1) as p:
    results = p.map(worker, [x for x in range(N)])

np.savetxt("results353_1.csv", 
           results[:10][0],
           delimiter =", ", 
           fmt ='% s')

np.savetxt("results353_2.csv", 
           results[10:20][0],
           delimiter =", ", 
           fmt ='% s')

np.savetxt("results353_3.csv", 
           results[20:30][0],
           delimiter =", ", 
           fmt ='% s')

np.savetxt("results353_4.csv", 
           results[30:40][0],
           delimiter =", ", 
           fmt ='% s')

np.savetxt("results353_5.csv", 
           results[40:50][0],
           delimiter =", ", 
           fmt ='% s')

np.savetxt("results353_6.csv", 
           results[50:60][0],
           delimiter =", ", 
           fmt ='% s')

np.savetxt("results353_7.csv", 
           results[60:70][0],
           delimiter =", ", 
           fmt ='% s')

np.savetxt("results353_8.csv", 
           results[70:80][0],
           delimiter =", ", 
           fmt ='% s')

np.savetxt("results353_9.csv", 
           results[80:90][0],
           delimiter =", ", 
           fmt ='% s')

np.savetxt("results353_10.csv", 
           results[90:100][0],
           delimiter =", ", 
           fmt ='% s')
