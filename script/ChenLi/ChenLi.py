# this is the script for LAMOST-RMS input catalog

import numpy as np
from astropy.io import fits

import sys
sys.path.append('../../util/')
import gcirc
import IO_InpCat as II

# ----------------------------------------------------------------
# read the central coordinates of the fibers 
finput = '../../../lamost_mrs/Test_SeptOct/Input_source/ChenLi/input_9-11_oc_prior5.fits'

# --------------------------------------------------
# read catalog from input catalog 
h1 = fits.open(finput)
t1 = h1[1].data

# all sources within the circle of radius = max_radius
# the following check skipped
'''
ra = t1['RAJ2000']
dec = t1['DEJ2000']

# center of the region, chosen to be the position of a bright star from Tycho-2 catalogue
center = np.array([19.359, 57.632])
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
'''

t1_sub = t1
g1 = t1_sub['GMAG']

id1 = (g1 >= 9.0) * (g1 <= 14.0)
t1_sub = t1_sub[id1]

source_id = t1_sub['SrCIDgaia']
ra = t1_sub['RAJ2000']
dec = t1_sub['DEJ2000']
g_mag = t1_sub['GMAG']
pri = t1_sub['PRIOR']
pri[pri < 3] = 1
pri[pri== 3] = 5

# output the catalog for lamost-sss
II.Inp_cat_output(OBJID=source_id, RADEG=ra, DECDEG=dec, MAG0=g_mag, PRI=pri)
