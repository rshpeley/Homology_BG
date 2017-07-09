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
import findDistStrRed as fDistPD
import sys


storage_home = '/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/newParams'

path1 = storage_home+'/scripts/' # Path to store the results
path2 = storage_home+'/Dist/' # Path to jdf files
path3 = storage_home+'/jdf/'
path5 = storage_home+'/output/' # Path for output 
pathname = storage_home+'/output/'

ipNo1 = 5
PDFPD = pickle.load(open(pathname+"PDFeaturesnew_"+str(ipNo1)+".pickle","r"))
clusterPD = pickle.load(open(pathname+"clustersPDMore"+str(ipNo1)+".pickle","r"))
labelPD = pickle.load(open(pathname+"LabeledClustersPD"+str(ipNo1)+".pickle","r"))


pathname1 = storage_home+'/strictHealthy/output/'
PDFH = pickle.load(open(pathname1+"PDFeaturesnew_"+str(ipNo1)+".pickle","r"))
clusterH = pickle.load(open(pathname1+"clustersH"+str(ipNo1)+".pickle","r"))
labelH = pickle.load(open(pathname1+"LabeledClustersH"+str(ipNo1)+".pickle","r"))
sys.path.append(storage_home+"/strictHealthy/")
import findDist as fDistH
sys.path.insert(0, storage_home)

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
AllAsPD = pickle.load(open(pathname+"AllAs5.pickle","r"))
AllBsPD = pickle.load(open(pathname+"AllBs5.pickle","r"))

#AllAsH = pickle.load(open(pathname1+"AllAs5.pickle","r"))
#AllAsH = pickle.load(open(pathname1+"AllAsQ.pickle","r"))
AllAsH = pickle.load(open(pathname1+"AllAsQ1.pickle","r"))
#AllBsH = pickle.load(open(pathname1+"AllBs5.pickle","r"))
#AllBsH = pickle.load(open(pathname1+"AllBsQ.pickle","r"))
AllBsH = pickle.load(open(pathname1+"AllBsQ1.pickle","r"))

# RobustPathological
cluster1 = []
# Robust Healthy
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
	if np.abs(paramComb["d1ti"]) > np.abs(paramComb["d1ta"])/1.5  or np.abs(paramComb["d2ti"]) > np.abs(paramComb["d2ta"])/1.5  or np.abs(paramComb["tid2"]) < 1.2 or np.abs(paramComb["tad2"]) < 1.2  or np.abs(paramComb["stnta"]) > np.abs(paramComb["stnti"])/2. or np.abs(paramComb["fsita"]) > np.abs(paramComb["fsiti"])/1.5:	
		return dict(),dict(),dict(),0,np.zeros((13))
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
	#	Flags.append("Act")
	Flags.append("Trans")
	#	Flags.append("DopDep")
	TransDopDepRates = psGA.calcRates(Flags,delay,A,B,False,ipctx,iptau=p.params["iptau"],ipamp=5.)
	# Checks Refer Figure 6I and 6J of Mallet 2008 - All the below firing rates are from Abdi 2015

	tests = fDistPD.checkConds(SWADopDepRates,ActDopDepRates)
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
	if np.abs(paramComb["d1ti"]) > np.abs(paramComb["d1ta"])/1.5  or np.abs(paramComb["d2ti"]) > np.abs(paramComb["d2ta"])/1.5  or np.abs(paramComb["tid2"]) > 1.0 or np.abs(paramComb["tad2"]) > 1.0 or np.abs(paramComb["stnta"]) > np.abs(paramComb["stnti"])/2. or np.abs(paramComb["fsita"]) > np.abs(paramComb["fsiti"])/1.5 or np.abs(paramComb["d1ta"]) > np.abs(paramComb["d2ta"]):	
		return dict(),dict(),0,np.zeros((10))


	Flags = []
	Flags.append("SWA")
	ipctx=dict()
	ipctx["ip"]=np.zeros((1,2001))
	delay = 1
	SWARates = psGA.calcRates(Flags,delay,A,B,False,ipctx,iptau=p.params["iptau"])

	# Calculate Rates for Activation and Lesion
	Flags = []
	Flags.append("Act")
#	Flags.append("Trans")
#	Flags.append("DopDep")
	ActRates = psGA.calcRates(Flags,delay,A,B,False,ipctx,iptau=p.params["iptau"])
	print "SWA:d1,d2,fsi,ta,ti,stn,gpi,tha",np.mean(SWARates["d1"]),np.mean(SWARates["d2"]),np.mean(SWARates["fsi"]),np.mean(SWARates["ta"]),np.mean(SWARates["ti"]),np.mean(SWARates["stn"]),np.mean(SWARates["gpi"]),np.mean(SWARates["tha"])
	print "Act:d1,d2,fsi,ta,ti,stn,gpi,tha",np.mean(ActRates["d1"]),np.mean(ActRates["d2"]),np.mean(ActRates["fsi"]),np.mean(ActRates["ta"]),np.mean(ActRates["ti"]),np.mean(ActRates["stn"]),np.mean(ActRates["gpi"]),np.mean(ActRates["tha"])
	#Flags = []
#	Flags.append("Act")
	#Flags.append("Trans")
#	Flags.append("DopDep")
	#Rates = psGA.calcRates(Flags,delay,A,B,False,ipctx,iptau=p.params["iptau"],ipamp=5.)
	tests = fDistH.checkConds(SWARates,ActRates)
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

def checkPDF(Rates):
	PDFeatures = np.zeros((1,8))
	dt = p.params["dt"]
	i = 0
	howmuch = (np.mean(Rates['gpi'][100/dt:500/dt]) -np.mean(Rates['gpi'][600/dt:1400/dt]))/np.mean(Rates['gpi'][100/dt:500/dt])		# (Orig - Final)/Orig
	PDFeatures[i][0] = howmuch
	noiseLim = 0.5
	timeStart = 510
	#timeEnd = 1480	
	timeEnd1 = 1000 # This took a lot of time to tune, so there is a tradeoff between frequency resolution and time resolution,if frequency resolution increases, longer length of signal, temporal resolution deecreases, which means the peak in specrum becomes wider and shifts towards higher frequencies, ultimately slipping outside the beta band. If shorter length of signal is sent (time_range < len(Rates)), temporal resolution increases peak becomes thinner or more precise but then again shifts towards lower frequencies. So a compromise reached at 1000. For damped oscillations the signal length also needed to be decreased.  	
	timeEnd2 = 650   	
	# An alternatiev idea to SE could be , percentage of power contained in beta-band (periodogram)- because you require some value fot susceptibility to oscillations which varies between 0 and 1
	noise = np.random.uniform(-noiseLim,noiseLim,len(Rates['gpi'][timeStart/dt:timeEnd1/dt]))
        time = np.arange(0,len(noise),1)
        wind = np.exp(-time/10.)
        noise1 = np.convolve(noise,wind,mode='same')

	#Rates['stn'][timeStart/dt:timeEnd/dt] += noise
	#Rates['stn'][timeStart/dt:timeEnd1/dt] += noise1[:len(noise)]

	#se1,dfreq1 = cPDF.spec_entropy(Rates['stn'][timeStart/dt:timeEnd/dt],time_range=Rates["ipctx"][timeStart/dt:timeEnd/dt],freq_range=[0.015,0.025]) 
	se11,dfreq11,maxFreq11,perMax11 =  cPDF.spec_entropy(Rates['stn'][timeStart/dt:timeEnd1/dt],time_range=Rates["ipctx"][timeStart/dt:timeEnd1/dt],freq_range=[0.00009,0.00017]) # 0.00010 = 10 Hz, since resolution = 10^-5 sec, 1000msec, and dt = 0.01 , to better the resolution of frequency captured by fft
	se12,dfreq12,maxFreq12,perMax12 =  cPDF.spec_entropy(Rates['stn'][timeStart/dt:timeEnd2/dt],time_range=Rates["ipctx"][timeStart/dt:timeEnd2/dt],freq_range=[0.00007,0.00017]) # 0.00010 = 10 Hz, since resolution = 10^-5 sec, 1000msec, and dt = 0.01 , to better the resolution of frequency captured by fft
	se13,dfreq13,maxFreq13,perMax13 =  cPDF.spec_entropy(Rates['stn'][timeStart/dt:timeEnd1/dt],time_range=Rates["ipctx"][timeStart/dt:timeEnd1/dt],freq_range=[0.00017,0.00036]) # 0.00010 = 10 Hz, since resolution = 10^-5 sec, 1000msec, and dt = 0.01 , to better the resolution of frequency captured by fft		#se1,dfreq1,maxFreq1,perMax = spec_entropy(TransRates['stn'][timeStart/dt:timeEnd/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd/dt],freq_range=[0.00008,0.00020]) # 0.00010 = 10 Hz, since resolution = 10^-5 sec, 1000msec, and dt = 0.01 , to better the resolution of frequency captured by fft
	se14,dfreq14,maxFreq14,perMax14 = cPDF.spec_entropy(Rates['stn'][timeStart/dt:timeEnd2/dt],time_range=Rates["ipctx"][timeStart/dt:timeEnd2/dt],freq_range=[0.00017,0.00036]) # 0.00010 = 10 Hz, since resolution = 10^-5 sec, 1000msec, and dt = 0.01 , to better the resolution of frequency captured by fft
	ans = np.array([se11,se12,se13,se14])
	maxFreqs = np.array([maxFreq11,maxFreq12,maxFreq13,maxFreq14])
	perMaxs = np.array([perMax11,perMax12,perMax13,perMax14])
	indnan = np.where(np.isnan(ans)==True)
	ans[indnan] = 1
	indmin = np.where(ans==np.min(ans))
	print indmin
	#indmin = np.where(ans==np.min(ans))
	'''
	if se11 < se12:
		PDFeatures[i][2] = se11 
		maxFreq[i][0] = maxFreq11
		perPowerMaxFreq[i][0] = perMax11
	else:
		PDFeatures[i][2] = se12 
		maxFreq[i][0] = maxFreq12
		perPowerMaxFreq[i][0] = perMax12
	'''
	PDFeatures[i][2] = np.min(ans)
	#maxFreq[i][0] = maxFreqs[indmin]
	#perPowerMaxFreq[i][0] = perMaxs[indmin] 

	#Rates['ta'][timeStart/dt:timeEnd1/dt] += noise1[:len(noise)]
	se21,dfreq21,maxFreq21,perMax21 =  cPDF.spec_entropy(Rates['ta'][timeStart/dt:timeEnd1/dt],time_range=Rates["ipctx"][timeStart/dt:timeEnd1/dt],freq_range=[0.00009,0.00017])
	se22,dfreq22,maxFreq22,perMax22 =  cPDF.spec_entropy(Rates['ta'][timeStart/dt:timeEnd2/dt],time_range=Rates["ipctx"][timeStart/dt:timeEnd2/dt],freq_range=[0.00007,0.00017])
	se23,dfreq23,maxFreq23,perMax23 =  cPDF.spec_entropy(Rates['ta'][timeStart/dt:timeEnd1/dt],time_range=Rates["ipctx"][timeStart/dt:timeEnd1/dt],freq_range=[0.00017,0.00036])
	se24,dfreq24,maxFreq24,perMax24 = cPDF.spec_entropy(Rates['ta'][timeStart/dt:timeEnd2/dt],time_range=Rates["ipctx"][timeStart/dt:timeEnd2/dt],freq_range=[0.00017,0.00036])
	ans = np.array([se21,se22,se23,se24])
	maxFreqs = np.array([maxFreq21,maxFreq22,maxFreq23])
	perMaxs = np.array([perMax21,perMax22,perMax23])
	indnan = np.where(np.isnan(ans)==True)
	ans[indnan] = 1
	print ans
	indmin = np.where(ans==np.min(ans))[0]

	'''
	if se21 < se22:
		PDFeatures[i][3] = se21
		PDFeatures[i][4] = dfreq21
		maxFreq[i][1] = maxFreq21
		perPowerMaxFreq[i][1] = perMax21
	else:
		PDFeatures[i][3] = se22
		PDFeatures[i][4] = dfreq22
		maxFreq[i][1] = maxFreq22
		perPowerMaxFreq[i][1] = perMax22
	'''
	PDFeatures[i][3] = np.min(ans)
	#maxFreq[i][1] = maxFreqs[indmin]
	#perPowerMaxFreq[i][1] = perMaxs[indmin] 



	noise = np.random.uniform(-noiseLim,noiseLim,len(Rates['ti'][timeStart/dt:timeEnd1/dt]))
        #time = np.arange(0,len(noise),1)
        #wind = np.exp(-time/10.)
        #noise1 = np.convolve(noise,wind,mode='full')

	#Rates['ti'][timeStart/dt:timeEnd/dt] += noise
	#Rates['ti'][timeStart/dt:timeEnd1/dt] += noise1[:len(noise)]
	#RatesSTP['ti'][502/dt:1400/dt] += noise

	#se3,dfreq3 = cPDF.spec_entropy(Rates['ti'][timeStart/dt:timeEnd/dt],time_range=Rates["ipctx"][timeStart/dt:timeEnd/dt],freq_range=[0.015,0.025])
	se31,dfreq31,maxFreq31,perMax31 =  cPDF.spec_entropy(Rates['ti'][timeStart/dt:timeEnd1/dt],time_range=Rates["ipctx"][timeStart/dt:timeEnd1/dt],freq_range=[0.00009,0.00017])
	se32,dfreq32,maxFreq32,perMax32 =  cPDF.spec_entropy(Rates['ti'][timeStart/dt:timeEnd2/dt],time_range=Rates["ipctx"][timeStart/dt:timeEnd2/dt],freq_range=[0.00007,0.00017])
	se33,dfreq33,maxFreq33,perMax33 =  cPDF.spec_entropy(Rates['ti'][timeStart/dt:timeEnd1/dt],time_range=Rates["ipctx"][timeStart/dt:timeEnd1/dt],freq_range=[0.00017,0.00036])
	se34,dfreq34,maxFreq34,perMax34 = cPDF.spec_entropy(Rates['ti'][timeStart/dt:timeEnd2/dt],time_range=Rates["ipctx"][timeStart/dt:timeEnd2/dt],freq_range=[0.00017,0.00036])
	#se3,dfreq3,maxFreq1,perMax = cPDF.spec_entropy(Rates['ti'][timeStart/dt:timeEnd/dt],time_range=Rates["ipctx"][timeStart/dt:timeEnd/dt],freq_range=[0.00008,0.00020])
	ans = np.array([se31,se32,se33,se34])
	maxFreqs = np.array([maxFreq31,maxFreq32,maxFreq33])
	perMaxs = np.array([perMax31,perMax32,perMax33])
	indnan = np.where(np.isnan(ans)==True)
	ans[indnan] = 1
	print ans
	indmin = np.where(ans==np.min(ans))[0]
	print indmin

	'''
	if se31 < se32:
		PDFeatures[i][5] = se31
		PDFeatures[i][6] = dfreq31
		maxFreq[i][1] = maxFreq31
		perPowerMaxFreq[i][1] = perMax31
	else:
		PDFeatures[i][5] = se32
		PDFeatures[i][6] = dfreq32
		maxFreq[i][1] = maxFreq32
		perPowerMaxFreq[i][1] = perMax32
	'''
	PDFeatures[i][5] = np.min(ans)
	#maxFreq[i][2] = maxFreqs[indmin]
	#perPowerMaxFreq[i][2] = perMaxs[indmin] 

	se41,dfreq41,maxFreq41,perMax41 =  cPDF.spec_entropy(Rates['gpi'][timeStart/dt:timeEnd1/dt],time_range=Rates["ipctx"][timeStart/dt:timeEnd1/dt],freq_range=[0.00009,0.00017])
	se42,dfreq42,maxFreq42,perMax42 =  cPDF.spec_entropy(Rates['gpi'][timeStart/dt:timeEnd2/dt],time_range=Rates["ipctx"][timeStart/dt:timeEnd2/dt],freq_range=[0.00007,0.00017])
	se43,dfreq43,maxFreq43,perMax43 =  cPDF.spec_entropy(Rates['gpi'][timeStart/dt:timeEnd1/dt],time_range=Rates["ipctx"][timeStart/dt:timeEnd1/dt],freq_range=[0.00017,0.00036])
	se44,dfreq44,maxFreq44,perMax44 = cPDF.spec_entropy(Rates['gpi'][timeStart/dt:timeEnd2/dt],time_range=Rates["ipctx"][timeStart/dt:timeEnd2/dt],freq_range=[0.00017,0.00036])
	ans = np.array([se41,se42,se43,se44])
	indnan = np.where(np.isnan(ans)==True)
	ans[indnan] = 1
	print ans
	indmin = np.where(ans==np.min(ans))[0]
	PDFeatures[i][7] = np.min(ans)

	return PDFeatures




def replace(params):
	# Shift from here till dic KSValInfo dict to runrHBM.py and run separate simulation for each parameter

	clusPD = params["clusPD"]
	clusH = params["clusH"]

	idPD = params["idPD"]
	idH = params["idH"]
	MeanPD = params["MeanPD"]
	MeanH = params["MeanH"]
	conn = params["conn"]
	CVPD = params["CVPD"]	
	CVH = params["CVH"]	
	randSampsPD = params["randSampsPD"]
	randSampsH = params["randSampsH"]

	# Lets look at just 1st 5 K-S values
	print "MeanPD",MeanPD
	print "MeanH",MeanH # Confirming indeed they are vrey similar		
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
	#for samp in clusPD:	
	# Only temporarily commented, because wanted healthy results urgently
	
	for i,x in enumerate(randSampsPD):	
		temp = dict()
		A = np.array(AllAsPD)[x]
		B = np.array(AllBsPD)[x]
		print "ind1",ind1
		if ind1 < 16:
			A[ids[ind1][0],ids[ind1][1]] = MeanPD  # Replacing by each other's mean, if truly interchageable, this wont matter. Replace by its own mean, will lead to 500 identical models.
 		else:
			B[0,ids[ind1][0]] = MeanPD
		
		# Now check validity
		SWA,Act,Trans,Grades,tests = checkValidityPD(A,B)	
		
		# Now check which cluster	
		#PDF = checkPDF(Trans)
		#PDF1 = fs.postProcess(PDF[0],1)
	
		#temp["id"] = samp
		temp["id"] = i
		temp["Grades"] = Grades
		temp["tests"] = tests
		#temp["PDFPP"] = PDF1
		#temp["PDF"] = PDF
		if Grades == 13:
			samplesPD.append(temp)
		else:
			samplesChangedPD.append(temp)
	
	print "Healthy-----------------------------------------"
	# Then all healthy samples
	#for i in xrange(len(AllAsH))[:1000]:
	for i,x in enumerate(randSampsH):
		temp1 = dict()
		A = np.array(AllAsH)[x]
		B = np.array(AllBsH)[x]

		SWA,Act,Grades,tests = checkValidityH(A,B)	
		print "Before replacing"
		#SWA,Act,Trans,Grades,tests = checkValidityH(A,B)	
		print "ind1",ind1
		if ind1 < 16:
			A[ids[ind1][0],ids[ind1][1]] = MeanH  # Replacing by 1, coul dbe replaced by either Mean2 , since K-S statistic is very small
 		else:
			B[0,ids[ind1][0]] = MeanH
		
		# Now check validity
		#SWA,Act,Trans,Grades,tests = checkValidityH(A,B)	
		SWA,Act,Grades,tests = checkValidityH(A,B)	

		# Now check which cluster	
		#PDF = checkPDF(Trans)
		#PDF1 = fs.postProcess(PDF[0],1)
	
		#temp1["id"] = samp1
		temp1["id"] = i
		temp1["Grades"] = Grades
		temp1["tests"] = tests
		#temp1["PDFPP"] = PDF1
		#temp1["PDF"] = PDF
		if Grades == 10:
			samplesH.append(temp1)
		else:
			samplesChangedH.append(temp1)

	# Record K-S value wise, index represents which K-S value is this for , 0 = smallest value, 1 = 2nd smallest value 
	
	SamplesAllPD = dict()
	SamplesAllPD["stillPD"] = samplesPD
	SamplesAllPD["changed"] = samplesChangedPD
	SamplesAllPD["otherParams"] = params
	#SamplesAllPD["randSamps"] = randSampsPD
	
	# Record K-S value wise, index represents which K-S value is this for , 0 = smallest value, 1 = 2nd smallest value 
	SamplesAllH = dict()
	SamplesAllH["stillPD"] = samplesH
	SamplesAllH["changed"] = samplesChangedH
	SamplesAllH["otherParams"] = params
	#SamplesAllH["randSamps"] = randSampsH
	ts = time.time()
	pickle.dump(SamplesAllPD,open(pathname+"RHBM_PD_Self_"+conn+datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d')+".pickle","w"))
	pickle.dump(SamplesAllH,open(pathname1+"RHBM_H_Self_"+conn+datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d')+".pickle","w")) # Just store whatever is done so far
	#pickle.dump(SamplesAllH,open(pathname1+"RHBM_H_Self_1_"+conn+".pickle","w")) # Just store whatever is done so far

		




			
		

