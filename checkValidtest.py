import pickle
import numpy as np
import pylab as pl
import matplotlib as mpl


import matplotlib.cm as cm
from pylab import *
# Since its a sparse matrix
import scipy as sp
import itertools
import scipy.cluster.hierarchy as sph
import os       
from pylab import *     
import knownUnknownParams as p
import paramsearchGA_DopDep_nonlinear as psGA
#import paramsearchGA_DopDep_nonlinearTest as psGA
import checkPDFeaturesStrRed as cPDF
import funcs as fs
import pylab as pl
from scipy.optimize import fsolve
from scipy.integrate import odeint
# These could be any files, all we need is an example matrix
#AllAs = pickle.load(open("AllAsDiffGeatures20free.pickle","r"))

PDTherapyRates=[]
# Checks if A,B follow all the constraints for PD




# Plot the phase potrait of STN-TA,TI in an attempt to understand
#print "A",A1.round(2)
#print "B",B1.round(2)

#A1 = A1*1.0
#B1 = B1*1.0

#print "A",A1.round(2)
#print "B",B1
def dispRates(R):
	t = np.arange(0,len(R['d1']),1)
	pl.figure()
	pl.plot(R['d1'],'b-',label='D1')
	pl.plot(R['d2'],'r-',label='D2')
	pl.plot(R['fsi'],'g-',label='FSI')
	pl.plot(R['ta'],'c-',label='ta')
	pl.plot(R['ti'],'m-',label='ti')
	pl.plot(R['stn'],'y-',label='stn')
	pl.plot(R['gpi'],'-',color='orange',label='gpi')
	pl.plot(R['tha'],'k--',label='tha')
	pl.plot(R['ipctx'],'k-',label='ctx')
	pl.legend()
	pl.show()

pathname = "/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/newParams/output/"
#pathname = "/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/output/"
AllAs = pickle.load(open(pathname+"AllAs5.pickle","r"))
#AllAs = pickle.load(open(pathname+"AllAs0.pickle","r"))
AllBs = pickle.load(open(pathname+"AllBs5.pickle","r"))
#AllBs = pickle.load(open(pathname+"AllBs0.pickle","r"))
clustersPD = pickle.load(open(pathname+"clustersPDMore5.pickle","r"))
labelsPD = pickle.load(open(pathname+"LabeledClustersPD5.pickle","r"))
knownparams = p.params["known"]
#PDFeatures = pickle.load(open(pathname+"PDFeaturesnew_5.pickle","r"))
PDFeatures = pickle.load(open(pathname+"PDFeaturesnew_5.pickle","r"))
ipctx1=dict()
ipctx1["ip"] = np.zeros((1,1001))

dt = p.params["dt"]
timeGrid = np.arange(0.0,len(ipctx1["ip"][0])+dt,dt)
time = np.arange(0,len(timeGrid),1)
ipctx = np.zeros((len(time)))

ids = [(0,3),(0,4),(1,3),(1,4),(2,3),(2,4),(3,1),(3,3),(3,4),(4,1),(4,3),(4,4),(5,3),(5,4),(4,5),(3,5),(0,0),(0,1),(0,2),(0,5)]
Rates = []
# Pick 2 examples from eah cluster
for i,clus1 in enumerate(clustersPD):
	if len(clus1) > 0:
		#randSamps = np.random.randint(0,len(clus1)-1,1)
		if len(clus1) > 1:
			randSamps = np.random.randint(0,len(clus1)-1,1)
			#randSamps = np.random.randint(0,len(AllAs)-1,1)
		else:
			randSamps = [0]
		randSamps=[0]
		print randSamps	
		for j in randSamps:
			print clus1[j]
			A = AllAs[clus1[j]]
			B = AllBs[clus1[j]]	
		print "Check validity"
		A[6][0] = p.params["known"]["gpid1"]
		A[6][5] = p.params["known"]["gpistn"]
		A[6][4] = p.params["known"]["gpiti"]
		A[0][0] = p.params["known"]["d1d1"]
		A[0][1] = p.params["known"]["d1d2"]
		A[1][0] = p.params["known"]["d2d1"]
		A[1][1] = p.params["known"]["d2d2"]
		print A[5,3]	
		A[5,3] = -0.001

		print A[5,3]	
		#A[4][1] = A[4][1]*1.9
		#B[0][0] = B[0][0]*1.5
		'''
		for j,x in enumerate(ids):
			if j < 16:
				A[x[0],x[1]] = (A[x[0],x[1]]/0.8)*0.75
			else:
				B[x[0],x[1]] = (B[x[0],x[1]]/0.8)*0.75
		'''
		SWA,Act,Trans = fs.checkValidity(A,B)
		#PDF = checkPDF(Trans)
		PDFDone =PDFeatures[clus1[j]] 
		PDFNow =fs.checkPDF(Trans) 

		dispRates(Trans)
		'''
		PDF1 = fs.postProcess(PDFDone,1)
		print "GS-Prev",PDF1[0]
		print "SO-Prev",1.-PDF1[1]

		PDF2 = fs.postProcess(PDFNow[0],1)
		print "GS-Now",PDF2[0]
		print "SO-Now",1.-PDF2[1]
		'''
		print labelsPD[i]["label"] 
		print len(clus1)

		#dispRates(SWA)
pl.show()

