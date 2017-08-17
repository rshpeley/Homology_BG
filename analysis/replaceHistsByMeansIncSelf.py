import numpy as np
import pickle
import numpy as np
import pylab as pl
import matplotlib as mpl
mpl.use('Agg')
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.fftpack as spfft

import matplotlib.cm as cm
from pylab import *
# Since its a sparse matrix
import scipy.sparse.linalg as sp
import itertools
import scipy.cluster.hierarchy as sph
import os
from pylab import *
from scipy.stats import ks_2samp
import funcs as fs
import paramsearchGA_DopDep_nonlinear as psGA
import knownUnknownParams as p
import checkPDFeaturesStrRed as cPDF
import funcs as fs
import sys

storage_home = os.getcwd() 

path1 = storage_home+'/scripts/' # Path to store the results
path2 = storage_home+'/Dist/' # Path to jdf files
path3 = storage_home+'/jdf/'
path5 = storage_home+'/output/' # Path for output 
pathname = storage_home+'/output/'

AllAsPD = pickle.load(open(storage_home+"../PD/output/AllAs.pickle","r"))
AllBsPD = pickle.load(open(storage_home+"../PD/output/AllBs.pickle","r"))

AllAsH = pickle.load(open(storage_home+"../healthy/output/AllAs.pickle","r"))
AllBsH = pickle.load(open(storage_home+"../healthy/output/AllBs.pickle","r"))

labels=["RobustHealthy","AR","TDMild","TDSevere","RobustPathological","SevereAkinesia","RH<->TDMild","TDMild<->TDSevere","RH<->AR","AR<->SA","SA<->RP","InBetween" ]

KSPos = np.zeros((6,1)) # K-S statistic
KSNeg = np.zeros((14,1))

MeansClus1Neg = np.zeros((14,1)) # Means
MeansClus2Neg = np.zeros((14,1))
MeansClus1Pos =  np.zeros((6,1))
MeansClus2Pos =  np.zeros((6,1))
VarsClus1Neg = np.zeros((14,1)) # Means
VarsClus2Neg = np.zeros((14,1))
VarsClus1Pos =  np.zeros((6,1))
VarsClus2Pos =  np.zeros((6,1))



terms = ["d1ta","d1ti","d2ta","d2ti","fsita","fsiti","tad2","tata","tati","tid2","tita","titi","stnta","stnti","tistn","tastn","d1ctx","d2ctx","fsictx","stnctx"]
# And their respetive ids in the matrix
ids = [(0,3),(0,4),(1,3),(1,4),(2,3),(2,4),(3,1),(3,3),(3,4),(4,1),(4,3),(4,4),(5,3),(5,4),(4,5),(3,5),(0,7),(1,7),(2,7),(5,7)]
posbin = np.arange(0,19,1)*0.8
negbin = np.arange(-8,0.5,0.3)*0.8

cluster1 = []
cluster2 = []

def checkValidityPD(A,B):
	print "A",A.round(2)
	print "B",B.round(2)
	Flags = []
	Flags.append("SWA")
	ipctx=dict()
	ipctx["ip"]=np.zeros((1,2001))
	delay = 1
	SWADopDepRates = psGA.calcRates(Flags,delay,A,B,False,ipctx,iptau=p.params["iptau"])
	Flags = []
	Flags.append("Act")
	ActDopDepRates = psGA.calcRates(Flags,delay,A,B,False,ipctx,iptau=p.params["iptau"])
	Flags = []
	Flags.append("Trans")
	TransDopDepRates = psGA.calcRates(Flags,delay,A,B,False,ipctx,iptau=p.params["iptau"],ipamp=5.)
	# Checks Refer Figure 6I and 6J of Mallet 2008 - All the below firing rates are from Abdi 2015

	tests = fs.checkCondsPD(SWADopDepRates,ActDopDepRates)
	Grades = np.sum(tests)
	if Grades == 13:
		print "Valid"
	else:
		print "Invalid"
		print "Grades",Grades


	return SWADopDepRates, ActDopDepRates, TransDopDepRates, Grades

def checkValidityH(A,B):
	print "A",A.round(2)
	print "B",B.round(2)
	Flags = []
	Flags.append("SWA")
	ipctx=dict()
	ipctx["ip"]=np.zeros((1,2001))
	delay = 1
	SWARates = psGA.calcRates(Flags,delay,A,B,False,ipctx,iptau=p.params["iptau"])

	Flags = []
	Flags.append("Act")
	ActRates = psGA.calcRates(Flags,delay,A,B,False,ipctx,iptau=p.params["iptau"])

	Flags = []
	Flags.append("Trans")
	Rates = psGA.calcRates(Flags,delay,A,B,False,ipctx,iptau=p.params["iptau"],ipamp=5.)


	tests = fs.checkCondsH(SWARates,ActRates)
	Grades = np.sum(tests)
	if Grades == 10:
		print "Valid"
	else:
		print "Invalid"
		print "Grades",Grades


	return SWARates, ActRates, Rates,Grades






def replace(params):
	# Shift from here till dic KSValInfo dict to runrHBM.py and run separate simulation for each parameter

	clusPD = params["clusPD"]
	clusH = params["clusH"]

	MeanPD = params["MeanPD"]
	MeanH = params["MeanH"]
	conn = params["conn"]
	tupleSorted = params["tupleSorted"]
	indexOrd = params["indicesOrder"]
	PDSorted = params["PDSorted"]
	HSorted = params["HSorted"]
	indexOrdPD = params["indicesOrdPD"]
	indexOrdH = params["indicesOrdH"]
	ind = params["ind"]	
	print "MeanPD",MeanPD
	print "MeanH",MeanH 		
	print "conn",conn
	print "tuple",tupleSorted
	print "indexOrd",indexOrd
	samplesPD = []
	samplesChangedPD = []
	samplesH = []
	samplesChangedH = []		
	
	# First all PD samples
	for i,k in enumerate(clusPD):	
		temp = dict()
		A = np.array(AllAsPD)[k]
		B = np.array(AllBsPD)[k]
		for l,x in enumerate(indexOrdPD): # In order of decreasing classification robustness
			if x < 16:
				A[ids[x][0],ids[x][1]] = MeanPD[x][0]  
			else:
				B[0,ids[x][0]] = MeanPD[x][0]
		
		# Now check validity
		SWA,Act,Trans,Grades = checkValidityPD(A,B)	
		

	
		temp["id"] = i
		temp["Grades"] = Grades
		if Grades == 13:
			samplesPD.append(temp)
		else:
			samplesChangedPD.append(temp)
	
	# Then all healthy samples
	for i,k in enumerate(randSampsH):
		temp1 = dict()
		A = np.array(AllAsH)[k]
		B = np.array(AllBsH)[k]
		#print "Before replacing"
		#SWA,Act,Trans,Grades = checkValidityH(A,B)
		for l,x in enumerate(indexOrdH):
			if x < 16:
				A[ids[x][0],ids[x][1]] = MeanH[x][0]  
			else:
				B[0,ids[x][0]] = MeanH[x][0]
		
		# Now check validity
		SWA,Act,Trans,Grades = checkValidityH(A,B)	

	
		temp1["id"] = i
		temp1["Grades"] = Grades
		if Grades == 10:
			samplesH.append(temp1)
		else:
			samplesChangedH.append(temp1)

	
	SamplesAllPD = dict()
	SamplesAllPD["stillPD"] = samplesPD
	SamplesAllPD["changed"] = samplesChangedPD
	SamplesAllPD["otherParams"] = params
	
	SamplesAllH = dict()
	SamplesAllH["stillPD"] = samplesH
	SamplesAllH["changed"] = samplesChangedH
	SamplesAllH["otherParams"] = params
	connStr = ''
	for x in conn:
		connStr+=x+"_"

	pickle.dump(SamplesAllPD,open(pathname+"RHBM_PD_"+connStr+"_Self_Inc_Sep.pickle","w"))
	pickle.dump(SamplesAllH,open(pathname1+"RHBM_H_"+connStr+"_Self_Inc_Sep.pickle","w")) 

		




			
		

