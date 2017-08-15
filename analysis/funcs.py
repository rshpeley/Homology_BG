import pickle
import numpy as np
import pylab as pl
import matplotlib as mpl
#mpl.use('Agg')
#from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.fftpack as spfft

import matplotlib.cm as cm
from pylab import *
# Since its a sparse matrix
import scipy.sparse.linalg as sp
import itertools
import scipy.cluster.hierarchy as sph
import os
from pylab import *
#import tree_anal
#import paramsearchGA_DopDep as psGA
import paramsearchGA_DopDep_nonlinear as psGA
import knownUnknownParams as p
import checkPDFeaturesStrRed as cPDF
#from sklearn.cluster import KMeans
knownparams = p.params["known"]


def postProcess(PDFeatures1,which):
        PDFeatures2 = np.copy(PDFeatures1)
        cols = np.shape(PDFeatures2)[1]
        for x in xrange(cols):
                indinf = np.where(np.isinf(PDFeatures2[:,x])==True)[0]
                if len(indinf) > 0:
                        PDFeatures2[indinf,x] = 0
                indnan = np.where(np.isnan(PDFeatures2[:,x])==True)[0]
                if len(indnan) > 0:
                        PDFeatures2[indnan,x] = 0

        indLN = np.where(PDFeatures2[:,0] < -1)[0]
        for x in indLN:
                PDFeatures2[x,0] = np.random.uniform(-0.75,-0.99,1)

        term1 = (PDFeatures2[:,2]+PDFeatures2[:,3]+PDFeatures2[:,5])/3.
        print term1

        PDFeatures2[:,1] = 1.-term1
        print "PDF",PDFeatures2[:,1]
        return PDFeatures2




def checkCondsPD(SWADopDepRates, ActDopDepRates):
	# Checks Refer Figure 6I and 6J of Mallet 2008
	tests = np.zeros(13)
	if np.round(np.mean(SWADopDepRates['ti'])) >= 19 and np.round(np.mean(SWADopDepRates['ti']))<=60: 				# Mean around 24, so +- 4Hz allowed
		tests[0] = 1
	
	if np.round(np.mean(ActDopDepRates['ti'])) >= 7 and np.round(np.mean(ActDopDepRates['ti'])) <=19:		# Mean around 14
		tests[1] = 1

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
	if SWADopDepRates['taFF'] > 1 or SWADopDepRates['tiFF'] > 1:
		tests[6] = 1 # Commented because no combination was fulfilling this. Moreover, this is no where specified in mallet.
	
	if np.mean(ActDopDepRates["ti"]) < np.mean(SWADopDepRates["ti"]) and np.mean(ActDopDepRates["ta"]) > np.mean(SWADopDepRates["ta"]):
		tests[9] = 1 
	# Sanity test , all rates are fairly above zero
	if np.mean(SWADopDepRates['d1'][100./p.params["dt"]:]) > 0.5 and np.mean(SWADopDepRates['d2'][100./p.params["dt"]:]) > 0.5 and np.mean(SWADopDepRates['fsi'][100./p.params["dt"]:]) > 0.5 and	np.mean(SWADopDepRates['ta'][100./p.params["dt"]:]) > 1.0 and np.mean(SWADopDepRates['ti'][100./p.params["dt"]:]) > 1.0 and np.mean(SWADopDepRates['stn'][100./p.params["dt"]:]) > 1.0 and np.mean(SWADopDepRates['gpi'][100./p.params["dt"]:]) > 0.1:
		tests[10] = 1
	tests[11] = 1
	tests[12] = 1


	return tests





def checkCondsH(SWARates,ActRates):
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

	if np.round(np.mean(SWARates['ta'])) >= 0 and np.round(np.mean(SWARates['ta'])) <= 5 :# and SWARates['taFF'] < 1.5:
		tests[2] = 1
	
	if np.round(np.mean(SWARates['ti'])) >= 9.5 and np.round(np.mean(SWARates['ti'])) <= 35:# and SWARates['tiFF'] < 1.5:
		tests[3] = 1 	
	
	if SWARates['taFF'] < 0.9 and SWARates['tiFF'] < 0.9:
		tests[4] = 1 # Refer to paramSearch_nonlinear.py in ../ directory
	
	if np.round(np.mean(ActRates['ti'])) >= 12 and np.round(np.mean(ActRates['ti'])) <=50:		# Mean around 14
		tests[5] = 1
	# Sanity test , all rates are fairly above zero
	if np.mean(SWARates['d1'][1000:]) > 1.0 and np.mean(SWARates['d2'][100./p.params["dt"]:]) > 0.5 and np.mean(SWARates['fsi'][1000:]) > 0.1 and 			np.mean(SWARates['ta'][1000:]) > 0.5 and np.mean(SWARates['ti'][1000:]) > 0.5 and np.mean(SWARates['stn'][1000:]) > 0.1 and np.mean(SWARates['gpi'][100./p.params["dt"]:]) > 0.1:
		tests[6] = 1

	if np.mean(ActRates["ti"]) >np.mean(SWARates["ti"]) and np.mean(ActRates["ta"]) > np.mean(SWARates["ta"]):
		tests[7] = 1

        if np.mean(SWARates['stn'])<=50: 
      		tests[8] = 1	
	if np.round(np.mean(ActRates['ta'])) >= 5 and np.round(np.mean(ActRates['ta'])) <=25: 
		tests[9] = 1

	return tests

