import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units

import sys
sys.path.append('../../util/')
import gcirc
import IO_InpCat as II
import choose_center

# ----------------------------------------------------------------
# input catalogue
#finput = 'K2_00_plate.csv'
#finput = 'K2_01_plate.csv'
finput = 'K2_04_plate.csv'
#finput = 'K2_05_plate.csv'
#finput = 'K2_06_plate.csv'
#finput = 'K2_08_plate.csv'
#finput = 'K2_13_plate.csv'
#finput = 'K2_14_plate.csv'
#finput = 'K2_16_plate.csv'
finput_supplement = 'K2_04_supplement.txt'
fgaia  = '../../../catalogue/gaia/base/base.fits'
#fgaia  = 'gaia2.fits'
ftycho = '../../../lamost_mrs/bt.fits'

# data format
dt = np.dtype([('OBJID', 'S19'), ('RADEG', float), ('DECDEG', float), ('MAG0', float), ('PRI', float)])
dts = np.dtype([('OBJID', 'S19'), ('RAh', float), ('RAm', float), ('RAs', float), ('DECd', float), ('DECm', float), ('DECs', float), ('MAG0', float), ('Comments', 'S19')])
# --------------------------------------------------
# read supplement catalog 

t1s = np.genfromtxt(finput_supplement, dtype=dts, names=True)
#sid1s = t1s['Epic_Num']
#sid1 = t1_sub['OBJID']
#id1 = np.in1d(sid1, sid1s, invert=True)
#t1_sub = t1_sub[id1]

ra  = (t1s['RAh'] + t1s['RAm']/60.0 + t1s['RAs']/3600.0) * 15.0
dec = t1s['DECd'] + t1s['DECm']/60.0 + t1s['DECs']/3600.0
g1s = t1s['Mag']
idx = (g1s >= 10.0) * (g1s <= 15.0)
t1s = t1s[idx]
ra  = ra[idx]
dec = dec[idx]

# choose 'best' center star from tycho-catalog
(cen_ra, cen_dec, idx, HIP, m) = choose_center.choose_center(ra, dec, ftycho)
print(cen_ra, cen_dec)
print(HIP, m)

t1s_sub = t1s[idx]
ra1s = ra[idx]
dec1s = dec[idx]
sid1s = t1s_sub['Epic_Num']
g1s = t1s_sub['Mag']
comments = t1s_sub['Comments']
pri = np.array([1]*m)
idx = (comments == 'RRLyr')
pri[idx] += 1
print(len(t1s_sub))
t1snew = np.empty(m, dtype=dt)
t1snew['OBJID'] = sid1s
t1snew['RADEG'] = ra1s
t1snew['DECDEG'] = dec1s
t1snew['MAG0'] = g1s
t1snew['PRI'] = pri

# --------------------------------------------------
# read input catalog 
t1 = np.genfromtxt(finput, dtype=dt, delimiter=',', names=True)
sid1 = t1['OBJID']
id1 = np.in1d(sid1, sid1s, invert=True)
t1 = t1[id1]
ra  = t1['RADEG']
dec = t1['DECDEG']

# center of the region, chosen to be the position of a bright star from Tycho-2 catalogue
#center = np.array([ 96.82913, 21.29358]) # plate00, HIP30726
#center = np.array([172.39918,  3.77993]) # plate01, HIP56073
#center = np.array([ 60.27909, 20.19916]) # plate04, HIP18762
#center = np.array([134.57621, 18.30828]) # plate05, HIP44056
#center = np.array([206.80585, -9.70945]) # plate06, HIP67271
#center = np.array([ 14.67684,  4.20579]) # plate08, HIP4582
#center = np.array([ 16.21933,  5.65626]) # plate08, HIP5074

#center = np.array([ 14.58106,  6.84431]) # plate08, HIP4559
#center = np.array([ 16.45516,  4.90837]) # plate08, HIP5141
#center = np.array([ 18.12769,  2.47157]) # plate08, HIP5646

#center = np.array([ 74.02765, 22.57655]) # plate13, HIP22935
#center = np.array([158.75899,  8.65043]) # plate13, HIP51802
#center = np.array([134.39545, 16.23365]) # plate13, HIP44000
#cen_ra = center[0]
#cen_dec = center[1]

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

t1_sub['PRI'] += 2
t1_sub = np.hstack([t1snew, t1_sub])

# --------------------------------------------------
# read catalog from gaia
h2 = fits.open(fgaia)
t2 = h2[1].data

ra  = t2['ra']
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
#d = np.rad2deg(d)
id2 = (d <= max_radius)

t2_sub = t2[id1][id2]

# remove duplicate source (if necessary)
# using cross-matching
ra1  = t1_sub['RADEG']
dec1 = t1_sub['DECDEG']
ra2  = t2_sub['ra']
dec2 = t2_sub['dec']
c1 = SkyCoord(ra=ra1*units.degree, dec=dec1*units.degree)
c2 = SkyCoord(ra=ra2*units.degree, dec=dec2*units.degree)
idx, d2d, d3d = SkyCoord.match_to_catalog_sky(c1, c2) # idx is index, type int
id1 = (d2d <= 1.0 * units.arcsec)
id2 = idx[id1]
id3 = np.array([True]*len(c2))
id3[id2] = False
t2_sub = t2_sub[id3]

# --------------------------------------------------
# select g-mag range [10, 15]
g1 = t1_sub['MAG0']
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
