import pickle
import numpy as np
import pylab as pl
import matplotlib as mpl
#mpl.use('Agg')


import matplotlib.cm as cm
from pylab import *
# Since its a sparse matrix
import scipy.sparse.linalg as sp
import scipy.fftpack as spfft
import itertools
import scipy.cluster.hierarchy as sph
import os
from pylab import *
#import paramsearchGA_DopDep as psGA
import paramsearchGA_DopDep_nonlinear as psGA
import knownUnknownParams as p

AllAs=[]
AllASTPs=[]
AllBs=[]
storage_home = '/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/newParams/strictHealthy'
from scipy.optimize import fsolve
ctx = 5.0

def S(x,theta,Qmax):
	sigma=3.8
	funcVal = Qmax*(1./(1.+np.exp(-(x-theta)/sigma)))
	return funcVal 	
	#return x

def calcMean():
	pathname = storage_home+'/output/'
	print "in calcMean"
	knownparams = p.params["known"] 
	leak = -0.05
	Supressed=[]
	Sup_Ind=[]
	k = 0
	#print Allcombs
	#PDFeatures = np.zeros((len(Allcombs),2))
	# 0= akinesia, 1,2 = gpi amplitude spectrum,dominant frequency, 3,4 = TA amplitude,frequency, 5,6=TI amplitude, frequency
#	Samples = np.arange(0,len(Allcombs),1)
#	np.random.shuffle(Samples)
	randSamps = pickle.load(open(pathname+"randSamps.pickle","r"))
	Allcombs = pickle.load(open(pathname+"Allcombs.pickle","r"))

	MeansSWA = []
	MeansAct = []
	MeansTrans=[]
	tempSWA=[]
	#tempSWATI=[]
	tempAct=[]
	tempTrans=[]

	for i,j in enumerate(randSamps):
		print i
		ind = Allcombs[j]
		d1ta = ind[0]
		d2ta = ind[1]
		fsita = ind[2]
		fsiti = ind[3]
		tata = ind[4]
		tati = ind[5]
		tastn = ind[6]
		tita = ind[7]
		titi = ind[8]
		tistn = ind[9]
		stnta = ind[10]
		stnti = ind[11]			
		tid2 = ind[12]
		tad2 = ind[13]	
		d1ti = ind[14]
		d2ti = ind[15]
		jc1 = ind[16]
		jc2 = ind[17]
		jfsictx = ind[18]
		jstnctx = ind[19]

		#A = np.matrix([[knownparams['d1d1'],knownparams['d1d2'],knownparams['d1fsi'],d1ta,knownparams['d1ti'],0.,0.],[knownparams['d2d1'],knownparams['d2d2'],knownparams['d2fsi'],d2ta,knownparams['d2ti'],0.,0.],[0.,0.,knownparams['fsifsi'],fsita,fsiti,0.,0.],[0,tad2,0.,tata,tati,tastn,0],[0.,tid2,0.,tita,titi,tistn,0.],[0.,0.,0.,knownparams['stnta'],stnti,knownparams['stnstn'],0.],[knownparams['gpid1'],0.,0.,knownparams['gpita'],knownparams['gpiti'],knownparams['gpistn'],knownparams['gpigpi']]])
		A = np.matrix([[knownparams['d1d1'],knownparams['d1d2'],knownparams['d1fsi'],d1ta,d1ti,0.,0.],[knownparams['d2d1'],knownparams['d2d2'],knownparams['d2fsi'],d2ta,d2ti,0.,0.],[0.,0.,knownparams['fsifsi'],fsita,fsiti,0.,0.],[0,tad2,0.,tata,tati,tastn,0],[0.,tid2,0.,tita,titi,tistn,0.],[0.,0.,0.,stnta,stnti,knownparams['stnstn'],0.],[knownparams['gpid1'],0.,0.,knownparams['gpita'],knownparams['gpiti'],knownparams['gpistn'],knownparams['gpigpi']]])

		B = np.matrix([jc1,jc2,jfsictx,0,0,jstnctx,0])
		delay = 1.0	
		#Calculate Rates for SWA and Control
		ipctx1=dict()
		#Calculate Rates for SWA and lesion(dopamine depletion)
		Flags = []
		Flags.append("SWA")
		ipctx1["ip"] = np.zeros((1,2001))
		SWADopDepRates = psGA.calcRates(Flags,delay,A,B,False,ipctx1,iptau=p.params["iptau"])
			
		tempSWA.append([np.mean(SWADopDepRates["d1"]),np.mean(SWADopDepRates["d2"]),np.mean(SWADopDepRates["fsi"]),np.mean(SWADopDepRates["ta"]),np.mean(SWADopDepRates["ti"]),np.mean(SWADopDepRates["stn"]),np.mean(SWADopDepRates["gpi"])])


		Flags = []
		Flags.append("Act")
		ipctx1["ip"] = np.zeros((1,2001))
		ActDopDepRates = psGA.calcRates(Flags,delay,A,B,False,ipctx1,iptau=p.params["iptau"])
	
		tempAct.append([ np.mean(ActDopDepRates["d1"]),np.mean(ActDopDepRates["d2"]),np.mean(ActDopDepRates["fsi"]),np.mean(ActDopDepRates["ta"]),np.mean(ActDopDepRates["ti"]),np.mean(ActDopDepRates["stn"]),np.mean(ActDopDepRates["gpi"])])


		Flags = []
		Flags.append("Trans")
		ipctx1["ip"] = np.zeros((1,2001))
		TransRates = psGA.calcRates(Flags,delay,A,B,False,ipctx1,iptau=p.params["iptau"],ipamp=0.0)
		tempTrans.append([ np.mean(TransRates["d1"]),np.mean(TransRates["d2"]),np.mean(TransRates["fsi"]),np.mean(TransRates["ta"]),np.mean(TransRates["ti"]),np.mean(TransRates["stn"]),np.mean(TransRates["gpi"])])

	MeansSWA.append(tempSWA)
	MeansAct.append(tempAct)
	MeansTrans.append(tempTrans)
	#MeansSWA["ti"] = tempSWATI
	#MeansAct["ti"] = tempActTI

	Means = dict()
	Means["SWA"] = MeansSWA
	Means["Act"] = MeansAct
	Means["Trans"] = MeansTrans

	pickle.dump(Means,open(pathname+"MeansTight.pickle","w"))

	
