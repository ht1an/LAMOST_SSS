import numpy as np
from astropy.io import fits

import sys
sys.path.append('../../util/')
import gcirc
import IO_InpCat as II
import read_target

# ----------------------------------------------------------------
# input catalogue
finput = '../../../lamost_mrs/Test_SeptOct/Input_source/ChenLi/input_9-11_oc_prior5.fits'
# observed catalogue
directory = './b3/'
fb0 = directory + 'CATALOG.csv'
fb1 = directory + 'SKY101.TXT'
fb2 = directory + 'SKY102.TXT'

# --------------------------------------------------
# read input catalog 
h0 = fits.open(finput)
t0 = h0[1].data

# all sources within the the circle of radius = max_radius

sid0 = t0['SrCIDgaia']

# ---- ---- ---- ---- ---- ---- ---- ----
# read observed sources
dt1 = np.dtype([('RADEG', float), ('DECDEG', float), ('MAG0', float), ('PRI', int), ('OBJID', 'S19')])
xb0 = np.genfromtxt(fb0, delimiter=',', dtype=dt1, names=True)

xb1 = read_target.read_target(fb1)
xb2 = read_target.read_target(fb2)

sb0 = xb0['OBJID']
sb1 = xb1['OBJID']
sb2 = xb2['OBJID']

sid0_sub = sid0
id1 = np.in1d(sid0_sub, sb1)
sid0_sub = sid0_sub[~id1]
id2 = np.in1d(sid0_sub, sb2)

t0_sub = t0[~id1][~id2]

g = t0_sub['GMAG']
idx = (g >= 10.0) * (g <= 15.0)
t0_sub = t0_sub[idx]

source_id = t0_sub['SrCIDgaia']
ra = t0_sub['RAJ2000']
dec = t0_sub['DEJ2000']
g_mag = t0_sub['GMAG']
pri = t0_sub['PRIOR']
pri[pri < 3] = 1
pri[pri== 3] = 5

# output the catalog for lamost-sss
II.Inp_cat_output(OBJID=source_id, RADEG=ra, DECDEG=dec, MAG0=g_mag, PRI=pri)
