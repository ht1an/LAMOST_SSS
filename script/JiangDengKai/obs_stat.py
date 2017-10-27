import numpy as np
from astropy.io import fits

import sys
sys.path.append('../../util/')
import gcirc
import IO_InpCat as II
import read_target

# ----------------------------------------------------------------
# input catalogue
finput = '../../../lamost_mrs/Test_SeptOct/Input_source/JiangDengKai/region_plane_new.fits'
finput_supplement = '../../../lamost_mrs/Test_SeptOct/Input_source/JiangDengKai/lamost_eb_sb_gaia.fits'

# observed catalogue
directory = './'
fobs = directory + 'SKY101.TXT'

# --------------------------------------------------
# read input catalog 
h0 = fits.open(finput)
t0 = h0[1].data

# read supplement catalog
h1 = fits.open(finput_supplement)
t1 = h1[1].data
g1 = t1['phot_g_mean_mag']
idx = (g1 >= 10.0) * (g1 <= 15.0)
t1 = t1[idx]
sid1 = t1['gaia_id']
sid1_fmt = np.array([int(x[1:]) for x in sid1])
sid1 = sid1_fmt

t2 = read_target.read_target(fobs)
sid2 = t2['OBJID']
nobs = np.sum(np.in1d(sid1, sid2))
print(nobs)
