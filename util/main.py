# this is the script for LAMOST]

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import math
import IO_InpCat as II
import csv
# import pyfits as pf
import random
# ----------------------------------------------------------------
# read the central coordinates of the fibers 
dpath_SSS = "/Users/naoc/Documents/tian/LAMOST/InputCatalog/LAMOST2-MRS/Test_SeptOct/SSS/"
dpath_Inp = "/Users/naoc/Documents/tian/LAMOST/InputCatalog/LAMOST2-MRS/Test_SeptOct/Input_source/"
dpath_Outp = "/Users/naoc/Documents/tian/LAMOST/"
fn_Fib_Cat = "fiberdata.20140115.txt"
fn_gaia = dpath_Inp + "test/GaiaSource_000-004-225.dat"
fn_out = dpath_Outp + "ASS2SSS1.TXT"
fn_inp = dpath_Outp + "CATALOG.csv"
# print (fn_out)
# test_fits = dpath_Inp + "test/GaiaSource_000-004-225.fits"

# # 0 angle in polar frame for each fiber
# # 1 radii in polar frame for each fiber
# # 2 fibercell radii
# # 3 spectroscopy ID
# # 4 fiber id in each spectroscopy
# data = np.loadtxt(dpath_SSS+fn_Fib_Cat,skiprows=4)
# #  convert polar coordinates (theta,r) to rectangular coordinates (x,y)
# x = np.sin(data[:,0]*math.pi/180)*data[:,1]*(-1)
# y = np.cos(data[:,0]*math.pi/180)*data[:,1]*(-1)
# # plot the fibers
# plt.figure(figsize=(4,4))
# plt.scatter(x,y,s=10,linewidths=0.001,c=(data[:,3]*31)%255,cmap='gist_rainbow',edgecolors=None)
# for i_sp in xrange(1,17):
# 	ind = data[:,3]==i_sp
# 	plt.text(np.mean(x[ind])-0.25,np.mean(y[ind])-0.1,str(i_sp),fontsize=15)
# # 	plt.plot(x[ind],y[ind],'.')#,color=i_sp)
# for i_l in xrange(0,12):
# 	xl = np.sin(i_l*30*math.pi/180)*2.5*(-1)
# 	yl = np.cos(i_l*30*math.pi/180)*2.5*(-1)
# 	plt.plot([0,xl],[0,yl],'r--')

# ll = np.linspace(0,360,1000)
# for i_b in xrange(1,6):
# 	xl = np.sin(ll*math.pi/180)*i_b*0.5*(-1)
# 	yl = np.cos(ll*math.pi/180)*i_b*0.5*(-1)
# 	plt.plot(xl,yl,'r--')
# plt.text(-0.06,-2.8,"S")
# plt.text(-0.06,2.6,"N")
# plt.text(-2.8,-0.06,"E")
# plt.text(2.6,-0.06,"W")
# plt.savefig("fiber.eps")
# --------------------------------------------------
# read catalog from gaia
# print(fn_gaia)
Source_ID,ra,dec,MagG = np.loadtxt(fn_gaia,skiprows=1,dtype={'names':('Source_ID','ra','dec','MagG'),
	'formats':('19S',float,float,float)},unpack=True)
Source_ID_S = Source_ID.astype(str)
# Source_ID,ra,dec,MagG = np.loadtxt(fn_gaia,skiprows=1,dtype=None,unpack=True)
ind_mag14 = MagG<=14
# # II.Inp_cat_out(Source_ID[ind_mag14],ra[ind_mag14],dec[ind_mag14],z,MagG[ind_mag14],z,z,z,z,z,'Gaia_MagG14.csv')
# II.Inp_cat_out_NEW(Source_ID[ind_mag14],ra[ind_mag14],dec[ind_mag14],MagG[ind_mag14],'Gaia_MagG14_new.csv')

# II.Inp_cat_output(OBJID=Source_ID[ind_mag14],RADEG=ra[ind_mag14],DECDEG=dec[ind_mag14],MAG0=MagG[ind_mag14],fn='Gaia_MagG14_2.0.csv')
# TEST PRIORITY
nt = len(ra[ind_mag14])
PRI = np.random.randint(3,size=nt)
comment = 'THIS IS USED TO TEST PRIORITY, BUT NO PRIORITY'
II.Inp_cat_output(OBJID=Source_ID_S[ind_mag14],RADEG=ra[ind_mag14],DECDEG=dec[ind_mag14],MAG0=MagG[ind_mag14],fn='Gaia_MagG14_6.0.csv',PRI=PRI)
# PRI[PRI>1]=PRI[PRI>1]*3
# II.Inp_cat_output(OBJID=Source_ID_S[ind_mag14],RADEG=ra[ind_mag14],DECDEG=dec[ind_mag14],MAG0=MagG[ind_mag14],fn='Gaia_MagG14_7.0.csv',PRI=PRI)

# ----------------------------------------------------
# # read the output catalog from SSS
# SID,rao,deco,MagGo,pmrao,pmdeo, LAMBDA,GAMMA,x1,y1,x2,y2,x3,y3,pri= np.loadtxt(fn_out,skiprows = 12,
# 	dtype={'names':('Source_ID','ra','dec','MagG','pmra','pmde','LAMBDA','GAMMA',
# 		'X1','Y1','X2','Y2','X3','Y3','PRI'),
# 	'formats':('S19','f10','f10','f5','f5','S8','f11','f9','f13','f13','f13','f13','f13','f13','i2')},unpack=True)
# print SID[0],rao[0],deco[0],LAMBDA[0],GAMMA[0]
# # convert lambda gamma to x,y
# xo = np.sin(LAMBDA*math.pi/180)*GAMMA*(-1)
# yo = np.cos(LAMBDA*math.pi/180)*GAMMA*(-1)
# plt.figure(figsize=(4,4))
# for i_l in xrange(0,12):
# 	xl = np.sin(i_l*30*math.pi/180)*2.5*(-1)
# 	yl = np.cos(i_l*30*math.pi/180)*2.5*(-1)
# 	plt.plot([0,xl],[0,yl],'r--')

# ll = np.linspace(0,360,1000)
# for i_b in xrange(1,6):
# 	xl = np.sin(ll*math.pi/180)*i_b*0.5*(-1)
# 	yl = np.cos(ll*math.pi/180)*i_b*0.5*(-1)
# 	plt.plot(xl,yl,'r--')
# ind_pri0 = pri==0
# print len(xo[ind_pri0])
# plt.plot(xo[ind_pri0],yo[ind_pri0],'k.',markersize=0.5)
# plt.scatter(x,y,s=10,linewidths=0.001,facecolor=(0,0,0,0),edgecolors="c")

# plt.savefig("star_fiber.eps")

# SIDI, rai, deci,MagGi = np.loadtxt(fn_inp,skiprows=1,dtype={'names':('Source_ID','ra','dec','MagG'),
# 	'formats':('S19','f13','f13','f8')},delimiter=',',unpack=True)
# print SIDI[0],rai[0],deci[0],MagGi[0]
# plt.figure(figsize=(4,4))
# plt.plot(rai,deci,'k.',markersize=0.1)
# plt.plot([115.0875338],[-11.7526143],'r+',markersize=2)
# plt.xlabel("RA")
# plt.ylabel("Dec")
# plt.savefig("datainput.eps")

