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
finput = './K1_Xiejiwei_plate.csv'
fgaia  = '../../../catalogue/gaia/base/base.fits'
#fgaia  = 'gaia2.fits'

# --------------------------------------------------
# read input catalog 
dt = np.dtype([('OBJID', 'S19'), ('RADEG', float), ('DECDEG', float), ('MAG0', float), ('PRI', float)])

t1 = np.genfromtxt(finput, dtype=dt, delimiter=',', names=True)
ra  = t1['RADEG']
dec = t1['DECDEG']

# center of the region, chosen to be the position of a bright star from Tycho-2 catalogue
#center = np.array([294.15823, 44.69494]) # HIP96459
#center = np.array([287.54983, 44.55998]) # HIP94174

#center = np.array([285.1426 , 39.84493]) # HIP93329
#center = np.array([283.83376, 43.94609]) # HIP92862
#center = np.array([285.54254, 47.28735]) # HIP93473
#center = np.array([289.21412, 46.9991 ]) # HIP94755
#center = np.array([289.95925, 43,55627]) # HIP94998
#center = np.array([290.11295, 39.9262 ]) # HIP95056
#center = np.array([295.76872, 39.99712]) # HIP97018
#center = np.array([294.79245, 43,8395 ]) # HIP96661
#center = np.array([296.86178, 47.90757]) # HIP97372

#center = np.array([288.47314, 48.34928]) #HIP94487
#center = np.array([286.57098, 41.41378]) #HIP93808
#center = np.array([292.49547, 46.94649]) #HIP95879
#center = np.array([290.98542, 43.38817]) #HIP95352
center = np.array([293.5758 , 41.92722]) #HIP96252
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

ra_sub1  = ra[id1]
dec_sub1 = dec[id1]

# accurate selection using great circle distance / sky distance
d = gcirc.gcirc(ra_sub1, dec_sub1, cen_ra, cen_dec, u=1)
#d = np.rad2deg(d)
id2 = (d <= max_radius)

t1_sub = t1[id1][id2]

# --------------------------------------------------
# read catalog from gaia
h2 = fits.open(fgaia)
t2 = h2[1].data

ra  = t2['ra']
dec = t2['dec']

# pre select to reduce the computational time
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
#d = np.rad2deg(d)
id2 = (d <= max_radius)

t2_sub = t2[id1][id2]

# --------------------------------------------------
# select g-mag range [10, 15]
g1 = t1_sub['MAG0']
g2 = t2_sub['phot_g_mean_mag']
idx = (g1 >= 9.0) * (g1 <= 13.5)
g1 = g1[idx]
t1_sub = t1_sub[idx]
idx = (g2 >= 9.0) * (g2 <= 13.5)
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

# ---- ---- ---- ---- ---- ---- ---- ----
ra1 = t1_sub['RADEG']
dec1 = t1_sub['DECDEG']
objid1 = t1_sub['OBJID']
g1 = t1_sub['MAG0']
pri1 = t1_sub['PRI']

ra2 = t2_sub['ra']
dec2 = t2_sub['dec']
sid2 = t2_sub['source_id']
g2 = t2_sub['phot_g_mean_mag']
pri2 = np.array([10]*len(t2_sub))

objid = np.hstack([objid1, sid2])
ra    = np.hstack([ra1, ra2])
dec   = np.hstack([dec1, dec2])
g_mag = np.hstack([g1, g2])
pri   = np.hstack([pri1, pri2])
pri = pri.astype(int)

# output the catalog for lamost-sss
II.Inp_cat_output(OBJID=objid, RADEG=ra, DECDEG=dec, MAG0=g_mag, PRI=pri)
