import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units

import sys
sys.path.append('../../util/')
import gcirc
import IO_InpCat as II

# ----------------------------------------------------------------
# input catalogue
# no input catalogue
fgaia  = '../../../catalogue/gaia/base/base.fits'
#fgaia  = 'gaia2.fits'

# --------------------------------------------------
# read input catalog 
# no input catalog

# center of the region, chosen to be the position of a bright star from Tycho-2 catalogue
center = np.array([294.15823, 44.69494]) # HIP96459
cen_ra = center[0]
cen_dec = center[1]

# max radius of the observed region
max_radius = 2.5

# --------------------------------------------------
# read catalog from gaia
h2 = fits.open(fgaia)
t2 = h2[1].data

ra  = t2['ra']
dec = t2['dec']

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

ra_sub1  = ra[id1]
dec_sub1 = dec[id1]

# accurate selection using great circle distance / sky distance
d = gcirc.gcirc(ra_sub1, dec_sub1, cen_ra, cen_dec, u=1)
d = np.rad2deg(d)
id2 = (d <= max_radius)

t2_sub = t2[id1][id2]

# --------------------------------------------------
# select g-mag range [10, 15]
g2 = t2_sub['phot_g_mean_mag']
idx = (g2 >= 10.0) * (g2 <= 15.0)
g2 = g2[idx]
t2_sub = t2_sub[idx]
print('star number: ', len(t2_sub))

# selection strategy: make g-mag close to uniform distributions
N0 = 10000

if (len(t2_sub) > N0):
    print('enough stars, make g-mag selection')
    # enough stars, make g-mag selection
    num_bin = 5 # 10.0 - 15.0
    bins = np.linspace(10.0, 15.0, num_bin+1)
    #bins[0] = 0.0
    num_star = np.floor(N0 / num_bin)

    h2 = np.array([num_star]*num_bin)
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

# ---- ---- ---- ---- ---- ---- ---- ----
ra2 = t2_sub['ra']
dec2 = t2_sub['dec']
sid2 = t2_sub['source_id']
g2 = t2_sub['phot_g_mean_mag']
pri2 = np.array([1]*len(t2_sub))

objid = sid2
ra    = ra2
dec   = dec2
g_mag = g2
pri   = pri2
pri = pri.astype(int)

# output the catalog for lamost-sss
II.Inp_cat_output(OBJID=objid, RADEG=ra, DECDEG=dec, MAG0=g_mag, PRI=pri)
