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
#AllAsPD = pickle.load(open(pathname+"AllAsPDShuff.pickle","r"))
AllAsPD = pickle.load(open(pathname+"AllAs5.pickle","r"))		# Only consider cluster 4
#AllBsPD_S = pickle.load(open(pathname+"AllBsPDShuff.pickle","r"))
#AllBsPD_S = pickle.load(open(pathname+"AllBsPDShuff1.pickle","r"))
AllBsPD = pickle.load(open(pathname+"AllBs5.pickle","r"))

#AllAsH_S= pickle.load(open(pathname1+"AllAsHShuff.pickle","r"))
#AllAsH_S= pickle.load(open(pathname1+"AllAsHShuff1.pickle","r"))
#AllAsH= pickle.load(open(pathname1+"AllAs5.pickle","r"))		# Only consider cluster 0, since shuffling among the whole array leads to different values
AllAsH= pickle.load(open(pathname1+"AllAsQ.pickle","r"))		# Only consider cluster 0, since shuffling among the whole array leads to different values
#AllBsH_S= pickle.load(open(pathname1+"AllBsHShuff.pickle","r"))
#AllBsH_S= pickle.load(open(pathname1+"AllBsHShuff1.pickle","r"))
#AllBsH= pickle.load(open(pathname1+"AllBs5.pickle","r"))
AllBsH= pickle.load(open(pathname1+"AllBsQ.pickle","r"))

# RobustPathological
cluster1 = []
# Robust Healthy
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
	#	Flags.append("Act")
	Flags.append("Trans")
	#	Flags.append("DopDep")
	TransDopDepRates = psGA.calcRates(Flags,delay,A,B,False,ipctx,iptau=p.params["iptau"],ipamp=5.)
	# Checks Refer Figure 6I and 6J of Mallet 2008 - All the below firing rates are from Abdi 2015

	tests = fDistPD.checkConds(SWADopDepRates,ActDopDepRates)
	Grades = np.sum(tests)
	if Grades == 13:
		print "Valid"
	else:
		print "Invalid"
		print "Grades",Grades

	'''
        tests = np.zeros(13)
        #if np.mean(SWADopDepRates['ti']) > np.mean(SWADopDepRates['ta']):
	if np.round(np.mean(SWADopDepRates['ti'])) >= 19 and np.round(np.mean(SWADopDepRates['ti']))<=35: 				# Mean around 24, so +- 4Hz allowed
		tests[0] = 1
	#if np.mean(ActDopDepRates['ta']) > np.mean(ActDopDepRates['ti']):
	#if np.mean(ActDopDepRates['ti']) >= 10 and np.mean(ActDopDepRates['ti']) <18:		# Mean around 14
	if np.round(np.mean(ActDopDepRates['ti'])) >= 7 and np.round(np.mean(ActDopDepRates['ti'])) <=19:		# Mean around 14
		tests[1] = 1
	#if np.mean(SWADopDepRates['ti']) > np.mean(ActDopDepRates['ti']):
	#if np.mean(ActDopDepRates['ta']) > 16 and np.mean(ActDopDepRates['ta']) <=23: # Mean around 19 # 
	if np.round(np.mean(ActDopDepRates['ta'])) >= 7 and np.round(np.mean(ActDopDepRates['ta'])) <=15: # Mean around 19
		tests[7] = 1
	if np.round(np.mean(SWADopDepRates['ta'])) >=1 and np.round(np.mean(SWADopDepRates['ta'])) <=6: # Mean is 12   < np.mean(ActDopDepRates['ta']):
		tests[8] = 1
	# Refer Fig 5 A,B of Mallet 2008
	GpeSWADopDep = np.mean(SWADopDepRates['ti']+SWADopDepRates['ta'])
	GpeActDopDep = np.mean(ActDopDepRates['ti']+ActDopDepRates['ta'])
	if GpeSWADopDep > GpeActDopDep:
		tests[2] = 1
	# Check if STn is in phase with ctx (in-phase)
	if SWADopDepRates['stn_ctx'] > 0 :
		tests[3] = 1
	#Check if STN is in phase with TA
	if SWADopDepRates['stn_ta'] > 0:# and SWADopDepRates['taFF']>1.5:
		tests[4] = 1
	# Check if STN is anti-phase with TI
	if SWADopDepRates['stn_ti'] < 0:# and SWADopDepRates['tiFF']>1.5:
		tests[5] = 1 	
	# Check if FSI activity is higher than Striatal MSN activity
	#if np.mean(SWADopDepRates['fsi']) > np.mean(SWADopDepRates['d1']) and np.mean(SWADopDepRates['fsi']) > np.mean(SWADopDepRates['d2']):
	if SWADopDepRates['taFF'] > 1 or SWADopDepRates['tiFF'] > 1:
		tests[6] = 1 # Commented because no combination was fulfilling this. Moreover, this is no where specified in mallet.
	# Since this is PD model, check D1 activity < D2 activity
	#if np.mean(SWADopDepRates['d2']) - np.mean(SWADopDepRates['d1']) > 1:
	tests[9] = 1 
	# Sanity test , all rates are fairly above zero
	if np.mean(SWADopDepRates['d1'][1000:]) > 0.5 and np.mean(SWADopDepRates['d2'][1000:]) > 1.0 and np.mean(SWADopDepRates['fsi'][1000:]) > 0.5 and	np.mean(SWADopDepRates['ta'][1000:]) > 1.0 and np.mean(SWADopDepRates['ti'][1000:]) > 1.0 and np.mean(SWADopDepRates['stn'][1000:]) > 1.0 and np.mean(SWADopDepRates['gpi'][1000:]) > 0.1:
		tests[10] = 1
	#if GpeSWADopDep > np.mean(SWADopDepRates['d1']) and GpeSWADopDep > np.mean(SWADopDepRates['d2']):
	tests[11] = 1
	#if  np.mean(SWADopDepRates['stn']) > 50:
	tests[12] = 1
        print "tests",tests
        print "SWA:d1,d2,fsi,ta,ti,stn,gpi,tha",np.mean(SWADopDepRates["d1"]),np.mean(SWADopDepRates["d2"]),np.mean(SWADopDepRates["fsi"]),np.mean(SWADopDepRates["ta"]),np.mean(SWADopDepRates["ti"]),np.mean(SWADopDepRates["stn"]),np.mean(SWADopDepRates["gpi"]),np.mean(SWADopDepRates["tha"])
        print "Act:d1,d2,fsi,ta,ti,stn,gpi,tha",np.mean(ActDopDepRates["d1"]),np.mean(ActDopDepRates["d2"]),np.mean(ActDopDepRates["fsi"]),np.mean(ActDopDepRates["ta"]),np.mean(ActDopDepRates["ti"]),np.mean(ActDopDepRates["stn"]),np.mean(ActDopDepRates["gpi"]),np.mean(ActDopDepRates["tha"])
	Grades = np.sum(tests)	
	if Grades == 13:
		print "Valid"
	else:
		print "Invalid"
		print "Grades",Grades
	'''
	return SWADopDepRates, ActDopDepRates, TransDopDepRates, Grades

def checkValidityH(A,B):
	print "A",A.round(2)
	print "B",B.round(2)
	Flags = []
#	Flags.append("DopDep")
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

	Flags = []
#	Flags.append("Act")
	Flags.append("Trans")
#	Flags.append("DopDep")
	Rates = psGA.calcRates(Flags,delay,A,B,False,ipctx,iptau=p.params["iptau"],ipamp=5.)


	tests = fDistH.checkConds(SWARates,ActRates)
	Grades = np.sum(tests)
	if Grades == 10:
		print "Valid"
	else:
		print "Invalid"
		print "Grades",Grades

	'''
	# Checks Refer Figure 6I and 6J of Mallet 2008
        tests = np.zeros(10)
        # Refer Fig 5 A,B of Mallet 2008
        GpeSWA = np.mean(SWARates['ti']+SWARates['ta'])
        GpeAct = np.mean(ActRates['ti']+ActRates['ta'])

        if GpeSWA < GpeAct:
                tests[0] = 1
        # Check if STn is in phase with ctx (in-phase)
        if SWARates['stn_ctx'] > 0 :
                tests[1] = 1

        #Check if TA is non-modulated 
        #if np.mean(SWARates['ta'][1000:]) > 0 and abs(SWARates['stn_ta']) <= 0.05:
        if np.round(np.mean(SWARates['ta'])) >= 0 and np.round(np.mean(SWARates['ta'])) <= 5. :# and SWARates['taFF'] < 1.5:
        #if np.mean(SWARates['ta']) < np.mean(ActRates['ta']) :# and SWARates['taFF'] < 1.5:
        #if np.mean(SWARates['ta']) - np.mean(ActRates['ta']) <=-1./1. :# and SWARates['taFF'] < 1.5:
                tests[2] = 1

        # Check if TI is non-modulated
        #if np.mean(SWARates['ti'][1000:]) > 0 and abs(SWARates['stn_ti']) <= 0.05:
        if np.round(np.mean(SWARates['ti'])) >= 9.5 and np.round(np.mean(SWARates['ti'])) <= 35:# and SWARates['tiFF'] < 1.5:
        #if np.mean(SWARates['ti']) < np.mean(ActRates['ti']):# and SWARates['tiFF'] < 1.5:
        #if np.mean(SWARates['ti']) - np.mean(ActRates['ti']) <=-1./1.:# and SWARates['tiFF'] < 1.5:
                tests[3] = 1

        # Check if FSI activity is higher than Striatal MSN activity
        #if np.mean(SWARates['fsi']) > np.mean(SWARates['d1']) and np.mean(SWARates['fsi']) > np.mean(SWARates['d2'])
        if SWARates['taFF'] < 1 and SWARates['tiFF'] < 1./1.:
                tests[4] = 1 # Refer to paramSearch_nonlinear.py in ../ directory

        # Check if D1 >= D2, since healthy state
        if np.round(np.mean(ActRates['ti'])) >= 12 and np.round(np.mean(ActRates['ti'])) <=50:             # Mean around 14
        #if np.mean(ActRates['ti']) > np.mean(ActRates['ta']):          # Mean around 14
        #if np.mean(ActRates['ti']) - np.mean(ActRates['ta']) >=1./1.:           # Mean around 14
                tests[5] = 1
        # Sanity test , all rates are fairly above zero
        if np.mean(SWARates['d1'][1000:]) > 1.0 and np.mean(SWARates['d2'][1000:]) > 0.5 and np.mean(SWARates['fsi'][1000:]) > 0.5 and                  np.mean(SWARates['ta'][1000:]) > 0.5 and np.mean(SWARates['ti'][1000:]) > 0.5 and np.mean(SWARates['stn'][1000:]) > 0.5 and np.mean(SWARates['gpi'][1000:]) > 0.1:
                tests[6] = 1
        if GpeSWA > np.mean(SWARates['d1']) and GpeSWA > np.mean(SWARates['d2']):
                tests[7] = 1

        #if abs(np.mean(SWARates['stn']) - np.mean(ActRates['stn']))<=10 and np.mean(SWARates['stn'])<=40: # Additonal constraint since STN firing rates were going very high. Also consult supplementary figure in Mallet 2008
        if np.mean(SWARates['stn'])<=50 and np.mean(SWARates['d1']) > 2*np.mean(SWARates['d2']) and np.mean(ActRates['d1']) >2*np.mean(ActRates['d2']): # Additonal constraint since STN firing rates were going very high. Also consult supplementary figure in Mallet 2008
      		tests[8] = 1
        if np.round(np.mean(ActRates['ta'])) >= 5 and np.round(np.mean(ActRates['ta'])) <=25: # Mean around 19
        #if np.mean(SWARates['ti']) > np.mean(SWARates['ta']): # Mean around 19
        #if np.mean(SWARates['ti']) - np.mean(SWARates['ta']) >=1./1.: # Mean around 19
                tests[9] = 1

        print "tests",tests
        print "SWA:d1,d2,fsi,ta,ti,stn,gpi,tha",np.mean(SWARates["d1"]),np.mean(SWARates["d2"]),np.mean(SWARates["fsi"]),np.mean(SWARates["ta"]),np.mean(SWARates["ti"]),np.mean(SWARates["stn"]),np.mean(SWARates["gpi"]),np.mean(SWARates["tha"])
        print "Act:d1,d2,fsi,ta,ti,stn,gpi,tha",np.mean(ActRates["d1"]),np.mean(ActRates["d2"]),np.mean(ActRates["fsi"]),np.mean(ActRates["ta"]),np.mean(ActRates["ti"]),np.mean(ActRates["stn"]),np.mean(ActRates["gpi"]),np.mean(ActRates["tha"])
	Grades = np.sum(tests)	

	if Grades == 10:
		print "Valid"
	else:
		print "Invalid"
		print "Grades",Grades
	'''
	return SWARates, ActRates, Rates,Grades

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
	sampPD = params["samplesPD"]
	sampH = params["samplesH"]
	#tupleSorted = params["tupleSorted"]
	indexOrd = params["indicesOrder"]
	#PDSorted = params["PDSorted"]
	#HSorted = params["HSorted"]
	indexOrdPD = params["indicesOrdPD"]
	indexOrdH = params["indicesOrdH"]
	randSampsPD = params["randSampsPD"]
	randSampsH = params["randSampsH"]
	seed1 = params["seed"]	
	AllAsPD_S = pickle.load(open(pathname+"AllAsPDShuff_"+str(seed1)+".pickle","r")) # 1,2,3, are different trials, instances of shuffling
	AllBsPD_S = pickle.load(open(pathname+"AllBsPDShuff_"+str(seed1)+".pickle","r"))

	AllAsH_S= pickle.load(open(pathname1+"AllAsHShuff_"+str(seed1)+".pickle","r"))
	AllBsH_S= pickle.load(open(pathname1+"AllBsHShuff_"+str(seed1)+".pickle","r"))

	ind1 = params["ind1"]	
	# print the parameter in unshuffled version
	#x1 = ids[ind1][0]
	#x2 = ids[ind1][1]
	print "parameter:",conn
	# Lets look at just 1st 5 K-S values
	print "MeanPD",MeanPD
	print "MeanH",MeanH # Confirming indeed they are vrey similar		
	print "conn",conn
	#print "tuple",tupleSorted
	print "indexOrd",indexOrd
	samplesPD = []
	samplesChangedPD = []
	samplesH = []
	samplesChangedH = []		
	# First all PD samples
	# Temporary commenting !! Needed healthy data urgently
	'''
	for i,k in enumerate(randSampsPD):	
		temp = dict()
		A = np.array(AllAsPD)[k]
		B = np.array(AllBsPD)[k]
		for l,x in enumerate(indexOrdPD):
			if x < 16:
				A[ids[x][0],ids[x][1]] = AllAsPD_S[k][ids[x][0],ids[x][1]]  # Replacing by 1, coul dbe replaced by either Mean2 , since K-S statistic is very small
			else:
				B[0,ids[x][0]] = AllBsPD_S[k][0,ids[x][0]]
		
		# Now check validity
		SWA,Act,Trans,Grades = checkValidityPD(A,B)	
		
		# Now check which cluster	
		#PDF = checkPDF(Trans)
		#PDF1 = fs.postProcess(PDF,1)
	
		temp["id"] = i
		temp["Grades"] = Grades
		#temp["PDFPP"] = PDF1
		#temp["PDF"] = PDF
		if Grades == 13:
			samplesPD.append(temp)
		else:
			samplesChangedPD.append(temp)
	'''
	# Then all healthy samples
	for i,k in enumerate(randSampsH):
		temp1 = dict()
		A = np.array(AllAsH)[k]
		B = np.array(AllBsH)[k]
		print "Before replacing"

		SWA,Act,Trans,Grades = checkValidityH(A,B)	
		for l,x in enumerate(indexOrdH):
			if x < 16:
				A[ids[x][0],ids[x][1]] = AllAsH_S[k][ids[x][0],ids[x][1]]  # Replacing by 1, coul dbe replaced by either Mean2 , since K-S statistic is very small
			else:
				B[0,ids[x][0]] = AllBsH_S[k][0,ids[x][0]]
		
		# Now check validity
		SWA,Act,Trans,Grades = checkValidityH(A,B)	

		# Now check which cluster	
		#PDF = checkPDF(Trans)
		#PDF1 = fs.postProcess(PDF,1)
	
		temp1["id"] = i
		temp1["Grades"] = Grades
		#temp1["PDFPP"] = PDF1
		#temp1["PDF"] = PDF
		if Grades == 10:
			samplesH.append(temp1)
		else:
			samplesChangedH.append(temp1)

	# Record K-S value wise, index represents which K-S value is this for , 0 = smallest value, 1 = 2nd smallest value 
	'''
	SamplesAllPD = dict()
	SamplesAllPD["stillPD"] = samplesPD
	SamplesAllPD["changed"] = samplesChangedPD
	SamplesAllPD["otherParams"] = params
	'''
	#SamplesAllPD["randSamps"] = randSampsPD

	# Record K-S value wise, index represents which K-S value is this for , 0 = smallest value, 1 = 2nd smallest value 
	SamplesAllH = dict()
	SamplesAllH["stillPD"] = samplesH
	SamplesAllH["changed"] = samplesChangedH
	SamplesAllH["otherParams"] = params
	#SamplesAllH["randSamps"] = randSampsH

	connStr = ''
	for x in conn:
		connStr+=x+"_"

	#pickle.dump(SamplesAllPD,open(pathname+"RHBM_PD_"+connStr+"_Self_Inc_Sep_Shuff_"+str(seed1)+".pickle","w"))
	pickle.dump(SamplesAllH,open(pathname1+"RHBM_H_"+connStr+"_Self_Inc_Sep_Shuff_"+str(seed1)+".pickle","w")) # Just store whatever is done so far

		




			
		

