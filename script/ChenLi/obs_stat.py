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
directory = './f3/'
ff0 = directory + 'CATALOG.csv'
ff1 = directory + 'SKY101.TXT'
ff2 = directory + 'SKY102.TXT'
ff3 = directory + 'SKY103.TXT'
ff4 = directory + 'SKY104.TXT'

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

xf0 = np.genfromtxt(fb0, delimiter=',', dtype=dt1, names=True)
xf1 = read_target.read_target(ff1)
xf2 = read_target.read_target(ff2)
xf3 = read_target.read_target(ff3)
xf4 = read_target.read_target(ff4)

sb0 = xb0['OBJID']
sb1 = xb1['OBJID']
sb2 = xb2['OBJID']
sf0 = xf0['OBJID']
sf1 = xf1['OBJID']
sf2 = xf2['OBJID']
sf3 = xf3['OBJID']
sf4 = xf4['OBJID']

#sb0_sub = sb0
#id_b1 = np.in1d(sb0_sub, sb1)
#sb0_sub = sb0_sub[~id_b1]
#id_b2 = np.in1d(sb0_sub, sb2)
#sf0_sub = sf0
#id_f1 = np.in1d(sf0_sub, sf1)
#sf0_sub = sf0_sub[~id_f1]
#id_f2 = np.in1d(sf0_sub, sf2)
#sf0_sub = sf0_sub[~id_f2]
#id_f3 = np.in1d(sf0_sub, sf3)
#sf0_sub = sf0_sub[~id_f3]
#id_f4 = np.in1d(sf0_sub, sf4)

sid0_sub = sid0
obs_b1 = np.in1d(sid0_sub, sb1)
sid0_sub = sid0_sub[~obs_b1]
obs_b2 = np.in1d(sid0_sub, sb2)
sid0_sub = sid0_sub[~obs_b2]
obs_f1 = np.in1d(sid0_sub, sf1)
sid0_sub = sid0_sub[~obs_f1]
obs_f2 = np.in1d(sid0_sub, sf2)
sid0_sub = sid0_sub[~obs_f2]
obs_f3 = np.in1d(sid0_sub, sf3)
sid0_sub = sid0_sub[~obs_f3]
obs_f4 = np.in1d(sid0_sub, sf4)

p0 = t0['PRIOR']

tb0 = t0
tb1 = tb0[~obs_b1]
tb2 = tb1[~obs_b2]
tf0 = tb2
tf1 = tf0[~obs_f1]
tf2 = tf1[~obs_f2]
tf3 = tf2[~obs_f3]
tf4 = tf3[~obs_f4]

xb1_sub = tb0[obs_b1]
xb2_sub = tb1[obs_b2]
xf1_sub = tf0[obs_f1]
xf2_sub = tf1[obs_f2]
xf3_sub = tf2[obs_f3]
xf4_sub = tf3[obs_f4]

pb0 = xb0['PRI']
pb1 = xb1_sub['PRIOR']
pb2 = xb2_sub['PRIOR']
pf0 = xf0['PRI']
pf1 = xf1_sub['PRIOR']
pf2 = xf2_sub['PRIOR']
pf3 = xf3_sub['PRIOR']
pf4 = xf4_sub['PRIOR']

print('pri==1\ntotal: ',np.sum(p0==1), 'bt input: ', np.sum(pb0==1), 'bt output: ', np.sum(pb1==1), np.sum(pb2==1), 'ft input: ', np.sum(pf0==1), 'ft output: ', np.sum(pf1==1), np.sum(pf2==1), np.sum(pf3==1), np.sum(pf4==1))
print('pri==2\ntotal: ',np.sum(p0==2), 'bt input: ', np.sum(pb0==2), 'bt output: ', np.sum(pb1==2), np.sum(pb2==2), 'ft input: ', np.sum(pf0==2), 'ft output: ', np.sum(pf1==2), np.sum(pf2==2), np.sum(pf3==2), np.sum(pf4==2))
#print('pri==3\ntotal: ',np.sum(p0==3), 'bt input: ', np.sum(pb0==3), 'bt output: ', np.sum(pb1==3), np.sum(pb2==3), 'ft input: ', np.sum(pf0==3), 'ft output: ', np.sum(pf1==3), np.sum(pf2==3), np.sum(pf3==3), np.sum(pf4==3))
print('pri==20\ntotal: ',np.sum(p0==20), 'bt input: ', np.sum(pb0==20), 'bt output: ', np.sum(pb1==20), np.sum(pb2==20), 'ft input: ', np.sum(pf0==20), 'ft output: ', np.sum(pf1==20), np.sum(pf2==20), np.sum(pf3==20), np.sum(pf4==20))
