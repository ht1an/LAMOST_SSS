# this is the script for LAMOST-RMS input catalog

import sys
sys.path.append('../../util/')

import numpy as np
import matplotlib.pyplot as plt
import IO_InpCat as II
import csv
from astropy.io import fits
import gcirc

# ----------------------------------------------------------------
# input catalog
finput = '../../../lamost_mrs/Test_SeptOct/Input_source/JiangDengKai/region_plane_new.fits'
fgaia = '../../../catalogue/gaia/base/base.fits'

# --------------------------------------------------
# read input catalog 
h1 = fits.open(finput)
t1 = h1[1].data

ra = t1['ra']
dec = t1['dec']

# center of the region, chosen to be the position of a bright star from Tycho-2 catalogue
center = np.array([89.14070, 28.94226]) # HIP28117
cen_ra = center[0]
cen_dec = center[1]

# max radius of the observed region
max_radius = 2.5

# pre select to reduce the computational time
# lower and upper boundary for dec
dec_l = cen_dec - max_radius
dec_u = cen_dec + max_radius
if cen_dec >= 0:
    if (dec_u > 90.0):
        dec_u = 90.0
    dec_h = dec_u
else:
    if (dec_l < -90.0):
        dec_l = -90.0
    dec_h = -dec_l

dec_h0 = np.rad2deg(np.arccos(max_radius / 180.0)) # 89.204
if (dec_h >= dec_h0):
    ra_l = 0.0
    ra_u = 360.0
else:
    ra_l = cen_ra - max_radius / np.cos(np.deg2rad(dec_h))
    ra_u = cen_ra + max_radius / np.cos(np.deg2rad(dec_h))

if ra_l < 0.0: # 0.0 < ra_u < 360.0
    id1 = (dec >= dec_l) * (dec <= dec_u) * ((ra > 360.0 - ra_l) + (ra < ra_u))
elif ra_u > 360.0: # 0.0 < ra_l < 360.0
    id1 = (dec >= dec_l) * (dec <= dec_u) * ((ra > ra_l) + (ra < ra_u - 360.0))
else:
    id1 = (dec >= dec_l) * (dec <= dec_u) * (ra > ra_l) * (ra < ra_u)

ra_sub1 = ra[id1]
dec_sub1 = dec[id1]

# accurate selection using great circle distance / sky distance
d = gcirc.gcirc(ra_sub1, dec_sub1, cen_ra, cen_dec, u=1)
d = np.rad2deg(d)
id2 = (d <= max_radius)

t1_sub = t1[id1][id2]

# --------------------------------------------------
# read catalog from gaia
h2 = fits.open(fgaia)
t2 = h2[1].data

ra = t2['ra']
dec = t2['dec']

if ra_l < 0.0: # 0.0 < ra_u < 360.0
    id1 = (dec >= dec_l) * (dec <= dec_u) * ((ra > 360.0 - ra_l) + (ra < ra_u))
elif ra_u > 360.0: # 0.0 < ra_l < 360.0
    id1 = (dec >= dec_l) * (dec <= dec_u) * ((ra > ra_l) + (ra < ra_u - 360.0))
else:
    id1 = (dec >= dec_l) * (dec <= dec_u) * (ra > ra_l) * (ra < ra_u)

ra_sub1  = ra[id1]
dec_sub1 = dec[id1]

# accurate selection using great circle distance / sky distance
d = gcirc.gcirc(ra_sub1, dec_sub1, cen_ra, cen_dec, u=1)
d = np.rad2deg(d)
id2 = (d <= max_radius)

t2_sub = t2[id1][id2]

# remove duplicate source (if necessary)
# using source_id
sid1 = t1_sub['source_id']
sid2 = t2_sub['source_id']
id2 = ~np.in1d(sid2, sid1)
t2_sub = t2_sub[id2]

# selection strategy: make g-mag close to uniform distributions
g1 = t1_sub['phot_g_mean_mag']
g2 = t2_sub['phot_g_mean_mag']
idx = (g1 >= 10.0) * (g1 <= 15.0)
g1 = g1[idx]
t1_sub = t1_sub[idx]
idx = (g2 >= 10.0) * (g2 <= 15.0)
g2 = g2[idx]
t2_sub = t2_sub[idx]
print('star number: ',len(t1_sub), len(t2_sub))

# selection strategy: make g-mag close to uniform distributions
N0 = 10000

if (len(t1_sub) + len(t2_sub) > N0):
    print('enough stars, make g-mag selection')
    # enough stars, make g-mag selection
    num_bin = 5 # 10.0 - 15.0
    bins = np.linspace(10.0, 15.0, num_bin+1)
    #bins[0] = 0.0
    h1, b1 = np.histogram(g1, bins)
    num_star = np.floor(N0 / num_bin)

    h2 = num_star - h1
    id5 = np.array([], dtype=int)
    for i in range(num_bin):
        if (h2[i] > 0):
            idx = (g2 >= bins[i]) * (g2 < bins[i+1])
            n = np.sum(idx)
            if (n <= h2[i]):
                id4 = np.nonzero(idx)[0]
            else:
                id3 = np.nonzero(idx)[0]
                id4 = np.random.permutation(id3)[:int(h2[i])]
            id5 = np.hstack([id5, id4])
    t2_sub = t2_sub[id5]
    print('after selection, t2: ', len(t2_sub))

ra1 = t1_sub['ra']
dec1 = t1_sub['dec']
sid1 = t1_sub['source_id']
g1 = t1_sub['phot_g_mean_mag']
pri1 = t1_sub['priority']
idx = (pri1 <= 30)
pri1[idx] = 1
idx = (pri1 > 30) * (pri1 <= 60)
pri1[idx] = 2
idx = (pri1 > 60)
pri1[idx] = 3

ra2 = t2_sub['ra']
dec2 = t2_sub['dec']
sid2 = t2_sub['source_id']
g2 = t2_sub['phot_g_mean_mag']
pri2 = np.array([10]*len(t2_sub))

source_id = np.hstack([sid1, sid2])
ra = np.hstack([ra1, ra2])
dec = np.hstack([dec1, dec2])
g_mag = np.hstack([g1, g2])
pri = np.hstack([pri1, pri2])
pri = pri.astype(int)

# output the catalog for lamost-sss
II.Inp_cat_output(OBJID=source_id, RADEG=ra, DECDEG=dec, MAG0=g_mag, PRI=pri)
