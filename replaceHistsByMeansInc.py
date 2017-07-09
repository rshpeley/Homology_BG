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
storage_home = '/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams'

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
posbin = np.arange(0,19,1)*0.75
negbin = np.arange(-8,0.5,0.3)*0.75
AllAsPD = pickle.load(open(pathname+"AllAs5.pickle","r"))
AllBsPD = pickle.load(open(pathname+"AllBs5.pickle","r"))

AllAsH = pickle.load(open(pathname1+"AllAs5.pickle","r"))
AllBsH = pickle.load(open(pathname1+"AllBs5.pickle","r"))

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
	SWADopDepRates = psGA.calcRates(Flags,delay,A,B,False,ipctx,iptau=14.)
	Flags = []
	Flags.append("Act")
	ActDopDepRates = psGA.calcRates(Flags,delay,A,B,False,ipctx,iptau=14.)
	Flags = []
	#	Flags.append("Act")
	Flags.append("Trans")
	#	Flags.append("DopDep")
	TransDopDepRates = psGA.calcRates(Flags,delay,A,B,False,ipctx,iptau=14.,ipamp=5.)
	# Checks Refer Figure 6I and 6J of Mallet 2008 - All the below firing rates are from Abdi 2015
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
	SWARates = psGA.calcRates(Flags,delay,A,B,False,ipctx,iptau=10.)

	# Calculate Rates for Activation and Lesion
	Flags = []
	Flags.append("Act")
#	Flags.append("Trans")
#	Flags.append("DopDep")
	ActRates = psGA.calcRates(Flags,delay,A,B,False,ipctx,iptau=10.)

	Flags = []
#	Flags.append("Act")
	Flags.append("Trans")
#	Flags.append("DopDep")
	Rates = psGA.calcRates(Flags,delay,A,B,False,ipctx,iptau=10.,ipamp=5.)

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

	return SWARates, ActRates, Rates,Grades

def checkPDF(Rates):
	PDFeatures = np.zeros((1,7))
	dt = p.params["dt"]
	i = 0
	howmuch = (np.mean(Rates['gpi'][100/dt:500/dt]) -np.mean(Rates['gpi'][600/dt:1400/dt]))/np.mean(Rates['gpi'][100/dt:500/dt])		# (Orig - Final)/Orig
	PDFeatures[i][0] = howmuch
	start = 520
	noise = np.random.uniform(0,2,len(Rates['gpi'][start/dt:1400/dt]))
	#noise = 0
	Rates['stn'][start/dt:1400/dt] += noise

	se1,dfreq1,m1,n1 = cPDF.spec_entropy(Rates['stn'][start/dt:1400/dt],time_range=Rates["ipctx"][start/dt:1400/dt],freq_range=[0.0001,0.0003]) 
	PDFeatures[i][2] = se1 

	#PDFeatures[i][2] = dfreq		
	
	#slice = Rates['ta'][900:1600]
	#x = 100
	#checkZeroRates = [ 1 for i in xrange(len(slice)) if np.mean(slice[i:i+x]) < 1]
	#if np.sum(checkZeroRates) < 0.5*len(checkZeroRates):
	noise = np.random.uniform(0,2,len(Rates['ta'][start/dt:1400/dt]))
	#noise = 0
	Rates['ta'][start/dt:1400/dt] += noise

	se2,dfreq2,m2,n2 = cPDF.spec_entropy(Rates['ta'][start/dt:1400/dt],time_range=Rates["ipctx"][start/dt:1400/dt],freq_range=[0.0001,0.0003])
	PDFeatures[i][3] = se2
	PDFeatures[i][4] = dfreq2



	#slice = Rates['ti'][900:1600]
	#x = 100
	#checkZeroRates = [ 1 for i in xrange(len(slice)) if np.mean(slice[i:i+x]) < 1]
	#if np.sum(checkZeroRates) < 0.5*len(checkZeroRates):
	noise = np.random.uniform(0,2,len(Rates['ti'][start/dt:1400/dt]))
	#noise = 0
	Rates['ti'][start/dt:1400/dt] += noise

	se3,dfreq3,m3,n3 = cPDF.spec_entropy(Rates['ti'][start/dt:1400/dt],time_range=Rates["ipctx"][start/dt:1400/dt],freq_range=[0.0001,0.0003])

	PDFeatures[i][5] = se3
	PDFeatures[i][6] = dfreq3

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
	KSval = params["KS-val"]	
	# Lets look at just 1st 5 K-S values
	print "MeanPD",MeanPD
	print "MeanH",MeanH # Confirming indeed they are vrey similar		
	print "conn",conn
	print "KS-val",KSval
	print "VarPD",params["VarPD"]
	print "VarH",params["VarH"]
	ind1 = params["ind1"]	
	samplesPD = []
	samplesChangedPD = []
	samplesH = []
	samplesChangedH = []		
	# First all PD samples
	for samp in clusPD:	
		temp = dict()
		A = np.array(AllAsPD)[samp]
		B = np.array(AllBsPD)[samp]
		print "ind1",ind1
		for l,x in enumerate(ind1):
			if x < 16:
				A[ids[x][0],ids[x][1]] = MeanH[l]  # Replacing by 1, coul dbe replaced by either Mean2 , since K-S statistic is very small
			else:
				B[0,ids[x][0]] = MeanH[l]
		
		# Now check validity
		SWA,Act,Trans,Grades = checkValidityPD(A,B)	
		
		# Now check which cluster	
		PDF = checkPDF(Trans)
		PDF1 = fs.postProcess(PDF,1)
	
		temp["id"] = samp
		temp["Grades"] = Grades
		temp["PDFPP"] = PDF1
		temp["PDF"] = PDF
		if Grades == 13:
			samplesPD.append(temp)
		else:
			samplesChangedPD.append(temp)

	# Then all healthy samples
	for samp1 in clusH:
		temp1 = dict()
		A = np.array(AllAsH)[samp1]
		B = np.array(AllBsH)[samp1]

		print "ind1",ind1
		for l,x in enumerate(ind1):
			if x < 16:
				A[ids[x][0],ids[x][1]] = MeanPD[l]  # Replacing by 1, coul dbe replaced by either Mean2 , since K-S statistic is very small
			else:
				B[0,ids[x][0]] = MeanPD[l]
		
		# Now check validity
		SWA,Act,Trans,Grades = checkValidityH(A,B)	

		# Now check which cluster	
		PDF = checkPDF(Trans)
		PDF1 = fs.postProcess(PDF,1)
	
		temp1["id"] = samp1
		temp1["Grades"] = Grades
		temp1["PDFPP"] = PDF1
		temp1["PDF"] = PDF
		if Grades == 10:
			samplesH.append(temp1)
		else:
			samplesChangedH.append(temp1)

	# Record K-S value wise, index represents which K-S value is this for , 0 = smallest value, 1 = 2nd smallest value 
	SamplesAllPD = dict()
	SamplesAllPD["stillPD"] = samplesPD
	SamplesAllPD["changed"] = samplesChangedPD
	SamplesAllPD["otherParams"] = params

	# Record K-S value wise, index represents which K-S value is this for , 0 = smallest value, 1 = 2nd smallest value 
	SamplesAllH = dict()
	SamplesAllH["stillPD"] = samplesH
	SamplesAllH["changed"] = samplesChangedH
	SamplesAllH["otherParams"] = params

	connStr = ''
	for x in conn:
		connStr+=x+"_"

	pickle.dump(SamplesAllPD,open(pathname+"RHBM_PD_"+connStr+"_Inc.pickle","w"))
	pickle.dump(SamplesAllH,open(pathname1+"RHBM_H_"+connStr+"_Inc.pickle","w")) # Just store whatever is done so far

		




			
		

