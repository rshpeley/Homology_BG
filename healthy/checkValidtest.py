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
import checkPDFeaturesStrRed as cPDF
import funcs as fs
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
	pl.plot(t,R['d1'],'b-',label='D1')
	pl.plot(t,R['d2'],'r-',label='D2')
	pl.plot(t,R['fsi'],'g-',label='FSI')
	pl.plot(t,R['ta'],'c-',label='ta')
	pl.plot(t,R['ti'],'m-',label='ti')
	pl.plot(t,R['stn'],'y-',label='stn')
	pl.plot(t,R['gpi'],'-',color='orange',label='gpi')
	pl.plot(t,R['tha'],'k--',label='tha')
	pl.plot(t,R['ipctx'],'k-',label='ctx')
	pl.legend()
	pl.show()

pathname = "/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/newParams/strictHealthy/output/"
#pathname = "/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/strictHealthy/output/"

AllAs = pickle.load(open(pathname+"AllAs5.pickle","r"))
#AllAs = pickle.load(open(pathname+"AllAs0.pickle","r"))
AllBs = pickle.load(open(pathname+"AllBs5.pickle","r"))
#AllBs = pickle.load(open(pathname+"AllBs0.pickle","r"))
#clustersPD = pickle.load(open(pathname+"clustersH5.pickle","r"))
clustersPD = pickle.load(open(pathname+"clustersH5.pickle","r"))
labelsPD = pickle.load(open(pathname+"LabeledClustersH5.pickle","r"))
knownparams = p.params["known"]
PDFeatures = pickle.load(open(pathname+"PDFeaturesnew_5.pickle","r"))
ipctx1=dict()
ipctx1["ip"] = np.zeros((1,1001))

dt = p.params["dt"]
timeGrid = np.arange(0.0,len(ipctx1["ip"][0])+dt,dt)
time = np.arange(0,len(timeGrid),1)
ipctx = np.zeros((len(time)))
knownids = [(0,0),(0,1),(0,2),(1,0),(1,1),(1,2),(6,0),(6,3),(6,4),(6,5)]
knownterms = ["d1d1","d1d2","d1fsi","d2d1","d2d2","d2fsi","gpid1","gpita","gpiti","gpistn"]
Rates = []
# Pick 2 examples from eah cluster
for i,clus1 in enumerate(clustersPD[::1]):
	if len(clus1) > 0:
		#randSamps = np.random.randint(0,len(clus1)-1,1)
		if len(clus1) > 1:
			randSamps = np.random.randint(0,len(clus1)-1,1)
			#randSamps = [1]
		else:
			randSamps = [0]
		#randSamps = [0]
		print randSamps
		for j in randSamps:
			A = AllAs[clus1[j]]
			B = AllBs[clus1[j]]	
			#for k,x in enumerate(knownids):
			#		A[x[0],x[1]]=knownparams[knownterms[k]]
					
		print "Check validity"
		SWA,Act,Trans = fs.checkValidity(A,B)
		PDFDone =PDFeatures[clus1[j]] 
		PDFNow =fs.checkPDF(Trans) 

		dispRates(Trans)

		PDF1 = fs.postProcess(PDFDone,1)
		print "GS-Prev",PDF1[0]
		print "SO-Prev",1.-PDF1[1]

		PDF2 = fs.postProcess(PDFNow[0],1)
		print "GS-Now",PDF2[0]
		print "SO-Now",1.-PDF2[1]

		print labelsPD[i]["label"] 
		print len(clus1)
		#dispRates(SWA)
pl.show()

