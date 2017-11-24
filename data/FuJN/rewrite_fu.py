# this is used to rewrite the data of kepler area from Jianning Fu
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import math as m
import IO_InpCat_3 as II
import csv
import random

dpath = "/Users/naoc/Documents/tian/LAMOST/InputCatalog/LAMOST2-MRS/\
Test_SeptOct/Input_source/FuJianNing/"
fn_ind = ['00','01','04','05','06','08','13','14','16']


#center = np.array([ 96.82913, 21.29358]) # plate00, HIP30726
#center = np.array([172.39918,  3.77993]) # plate01, HIP56073
#center = np.array([ 60.27909, 20.19916]) # plate04, HIP18762
#center = np.array([134.57621, 18.30828]) # plate05, HIP44056
#center = np.array([206.80585, -9.70945]) # plate06, HIP67271
#center = np.array([ 14.67684,  4.20579]) # plate08, HIP4582
#center = np.array([ 74.02765, 22.57655]) # plate13, HIP22935
#center = np.array([158.75899,  8.65043]) # plate14, HIP51802
# center = np.array([134.39545, 16.23365]) # plate15, HIP44000
center_ra = np.array([ 96.82913,172.39918, 60.27909,134.57621,206.80585, 14.67684,\
	 74.02765,158.75899,134.39545])
center_dec = np.array([21.29358,  3.77993, 20.19916, 18.30828, -9.70945,  4.20579,\
	 22.57655,  8.65043, 16.23365])




for i in range(0,9):
	fnB = 'K2_'+fn_ind[i]+'_Bplate.txt'
	fnV = 'K2_'+fn_ind[i]+'_Vplate.txt'
	plate_ind = fn_ind[i]
	objidB, magB, raB, decB, priB, prgB = np.loadtxt(dpath+fnB,skiprows=1, \
		dtype={'names':('ID','mag','ra','dec','pri','prg'),
		'formats':('9S',float,float,float,float,'100S')},unpack=True)
	objidV, magV, raV, decV, priV, prgV = np.loadtxt(dpath+fnV,skiprows=1, \
		dtype={'names':('ID','mag','ra','dec','pri','prg'),
		'formats':('9S',float,float,float,float,'100S')},unpack=True)
	# for checking if there are common stars in both B and V plates
	n_com = 0
	for j in range(0,len(raB)):
		n_com = n_com+len(objidV[objidV==objidB[j]])
	print(i,'there are {n_com} commom in B and V plates'.format(**locals()))

	objid = np.concatenate((objidV,objidB))
	mag = np.concatenate((magV,magB))
	ra = np.concatenate((raV,raB))
	dec = np.concatenate((decV,decB))
	pri = np.concatenate((priV,priB))*0+1

	ind15 = mag<15
	II.Inp_cat_output(OBJID=objid[ind15].astype(str),RADEG=ra[ind15], \
		DECDEG=dec[ind15],MAG0=mag[ind15],\
		fn=dpath+'K2_{plate_ind}_plate.csv'.format(**locals()),PRI=pri[ind15])


	fig = plt.figure(figsize=(6,4))
	plt.plot(raV,decV,'r.',markersize=2,label="V-plate")
	plt.plot(raB,decB,'g.',markersize=2,label="B-plate")
	# px = [134.57621,125.85764]
	# py = [18.30828,18.13644]
	aaa = np.linspace(0,360,10000)*m.pi/180
	rrr = 2.5
	xxx = np.cos(aaa)*rrr
	yyy = np.sin(aaa)*rrr
	plt.plot(center_ra[i],center_dec[i],'k+',markersize=10)
	# plt.plot(px[1],py[1],'k+',markersize=10)
	plt.plot(xxx+center_ra[i],yyy+center_dec[i],'b-')
	# plt.plot(xxx+px[1],yyy+py[1],'b-')

	plt.xlabel("RA")
	plt.ylabel("dec")
	plt.legend(loc=1)
	# plt.show()
	plt.savefig(dpath+'ra_dec_{plate_ind}_plate.eps'.format(**locals()))
	print(type(objidV),'\\\\\\')