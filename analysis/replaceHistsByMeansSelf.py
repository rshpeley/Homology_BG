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
import datetime
import time
import sys
import funcs as fs

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
	paramComb = dict()
	for i in xrange(len(terms)):
		if i < 16:
			paramComb[terms[i]] = A[ids[i][0]][ids[i][1]]
		else:
			paramComb[terms[i]] = B[0,ids[i][0]]

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

	tests = fs.checkCondsPD(SWADopDepRates,ActDopDepRates)
	print "tests",tests
	Grades = np.sum(tests)
	if Grades == 13:
		print "Valid"
	else:
		print "Invalid"
		print "Grades",Grades


	return SWADopDepRates, ActDopDepRates, TransDopDepRates, Grades,tests

def checkValidityH(A,B):
	print "A",A.round(2)
	print "B",B.round(2)
	paramComb = dict()
	for i in xrange(len(terms)):
		if i < 16:
			paramComb[terms[i]] = A[ids[i][0]][ids[i][1]]
		else:
			paramComb[terms[i]] = B[0,ids[i][0]]
	print "paramComb",paramComb

	Flags = []
	Flags.append("SWA")
	ipctx=dict()
	ipctx["ip"]=np.zeros((1,2001))
	delay = 1
	SWARates = psGA.calcRates(Flags,delay,A,B,False,ipctx,iptau=p.params["iptau"])

	# Calculate Rates for Activation and Lesion
	Flags = []
	Flags.append("Act")
	ActRates = psGA.calcRates(Flags,delay,A,B,False,ipctx,iptau=p.params["iptau"])
	print "SWA:d1,d2,fsi,ta,ti,stn,gpi,tha",np.mean(SWARates["d1"]),np.mean(SWARates["d2"]),np.mean(SWARates["fsi"]),np.mean(SWARates["ta"]),np.mean(SWARates["ti"]),np.mean(SWARates["stn"]),np.mean(SWARates["gpi"]),np.mean(SWARates["tha"])
	print "Act:d1,d2,fsi,ta,ti,stn,gpi,tha",np.mean(ActRates["d1"]),np.mean(ActRates["d2"]),np.mean(ActRates["fsi"]),np.mean(ActRates["ta"]),np.mean(ActRates["ti"]),np.mean(ActRates["stn"]),np.mean(ActRates["gpi"]),np.mean(ActRates["tha"])

	tests = fs.checkCondsH(SWARates,ActRates)
	Grades = np.sum(tests)
	print "tests",tests
	if Grades == 10:
		print "Valid"
	else:
		print "Invalid"
		ind = np.where(tests==0)
		print "condition failed",ind
		print "Grades",Grades
	#return SWARates, ActRates, Rates,Grades, tests
	return SWARates, ActRates, Grades, tests






def replace(params):

	clusPD = params["clusPD"]
	clusH = params["clusH"]

	MeanPD = params["MeanPD"]
	MeanH = params["MeanH"]
	conn = params["conn"]
	CVPD = params["CVPD"]	
	CVH = params["CVH"]	

	print "MeanPD",MeanPD
	print "MeanH",MeanH 		
	print "conn",conn
	print "CVPD",CVPD
	print "CVH",CVH
	print "VarPD",params["VarPD"]
	print "VarH",params["VarH"]
	ind1 = params["ind1"]	
	samplesPD = []
	samplesChangedPD = []
	samplesH = []
	samplesChangedH = []		
	# First all PD samples
	
	for i,x in enumerate(clusPD):	
		temp = dict()
		A = np.array(AllAsPD)[x]
		B = np.array(AllBsPD)[x]
		print "ind1",ind1
		if ind1 < 16:
			A[ids[ind1][0],ids[ind1][1]] = MeanPD  # Replacing by mean
 		else:
			B[0,ids[ind1][0]] = MeanPD
		
		# Now check validity
		SWA,Act,Trans,Grades,tests = checkValidityPD(A,B)	
		

		temp["id"] = i
		temp["Grades"] = Grades
		temp["tests"] = tests
		if Grades == 13:
			samplesPD.append(temp)
		else:
			samplesChangedPD.append(temp)
	
	print "Healthy-----------------------------------------"

	for i,x in enumerate(clusH):
		temp1 = dict()
		A = np.array(AllAsH)[x]
		B = np.array(AllBsH)[x]

		#SWA,Act,Grades,tests = checkValidityH(A,B)	
		#print "Before replacing"
		print "ind1",ind1
		if ind1 < 16:
			A[ids[ind1][0],ids[ind1][1]] = MeanH  
 		else:
			B[0,ids[ind1][0]] = MeanH
		
		# Now check validity
		SWA,Act,Grades,tests = checkValidityH(A,B)	

		temp1["id"] = i
		temp1["Grades"] = Grades
		temp1["tests"] = tests
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

	ts = time.time()
	pickle.dump(SamplesAllPD,open(pathname+"RHBM_PD_Self_"+conn+datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d')+".pickle","w"))
	pickle.dump(SamplesAllH,open(pathname1+"RHBM_H_Self_"+conn+datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d')+".pickle","w")) # Just store whatever is done so far

		




			
		

