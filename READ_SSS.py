# THIS IS USED TO READ THE OUTPUT FILE FROM SSS
# read ASS2SSS1.TXT

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import math
import IO_InpCat as II
import csv
# import pyfits as pf
import random


outpath = "/Users/naoc/Documents/tian/LAMOST/output_fromSSS/TEST/"
fn_A = outpath + "ASS2SSS1_1.TXT"
fn_T = outpath + "test.dat"
fn_S = outpath + "SKY101_TESTPRI_1.TXT"
skip_row_A = 12
skip_row_S = 2
nrow_A = 0
nrow_S = 25
print(nrow_S,\
	nrow_A)
# for testing elements in each line
# nl=0
# with open(fn_A,'r') as FA:
# 	line = FA.readline()
# 	while line:
# 		data = line.split()
# 		print(nl,len(data))
# 		# print(data[])
# 		line = FA.readline()
# 		print("-------------")
# 		nl=nl+1


SourceID_A,ra_A,dec_A,Mag_A,PMARA_A, PMDECLASS_A, LAMBDA_A,GAMMA_A,X1_A,Y1_A,X_A,Y_A,X3_A,Y3_A,PRI_A = \
np.loadtxt(fn_A,skiprows = 12,dtype = {'names':('SourceID','ra','dec','Mag','PMRA','PMDE',\
'LAMBDA','GAMMA','X1','Y1','X','Y','X3','Y3','PRI'),'formats':('19S',float,float,float,float,'8S', \
float,float,float,float,float,float,float,float,int)},unpack=True)
print(SourceID_A[0],ra_A[0],dec_A[0])





