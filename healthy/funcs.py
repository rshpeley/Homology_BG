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
                #PDFeatures2[x,0] = (PDFeatures2[x,0]/(np.abs(PDFeatures2[x,0])*2.)+np.random.uniform(-0.5,0.5,1))
                PDFeatures2[x,0] = np.random.uniform(-0.75,-0.99,1)

#       for x in indLN:
                #PDFeatures2[x,0] = (PDFeatures2[x,0]/(np.abs(PDFeatures2[x,0])*2.)+np.random.uniform(-0.5,0.5,1))  
#               PDFeatures2[x,0] = (PDFeatures2[x,0]/maxGS)  
        # The idea is to not truncuate the values > -1 to exactly -1, but instead somewhere near -1
        # The idea of simply dropping the value > -1, leads to a discrepancy in indices in clusters - leads to mysteriosu appearnces of 
        # seeminly 'AR' points in TDMild clusters.
        #PDFeatures2[:,1] = (PDFeatures2[:,2]+PDFeatures2[:,3]+PDFeatures2[:,5])/3.
        term1 = (PDFeatures2[:,2]+PDFeatures2[:,3]+PDFeatures2[:,5])/3.
        #if len(PDFeatures2[0])>7:
        #       term1 = np.min((PDFeatures2[:,2],PDFeatures2[:,3],PDFeatures2[:,5],PDFeatures2[:,7]),axis=0)
        #else:
        #term1 = np.min((PDFeatures2[:,2],PDFeatures2[:,3],PDFeatures2[:,5]),axis=0)
        print term1
        if which == "PD":
                if len(PDFeatures2[0])>7:
                        term2 = (PDFeatures2[:,2]+PDFeatures2[:,3]+PDFeatures2[:,5] + PDFeatures2[:,7])/4.
                        print term2
                        PDFeatures2[:,1] = np.min( (term1,term2),0)
                        #print term1
                else:
                        PDFeatures2[:,1] = term1
        else:
                PDFeatures2[:,1] = term1
        #PDFeatures2[:,1] = (PDFeatures2[:,2]+PDFeatures2[:,3])/2.
        print "PDF",PDFeatures2[:,1]
        return PDFeatures2






def S(x,theta,Qmax):
	sigma=4.
	funcVal = Qmax*(1./(1.+np.exp(-(x-theta)/sigma)))
	return funcVal 	
	#return x


def mypsd(Rates,time_range,bin_w = 5., nmax = 4000):

      #psd, max_value, freq,h = misc2.psd_sp(spikes[:,1],nr_bins,nr_neurons)
      #print time_range
      bins = np.arange(0,len(time_range),1)
      #print bins
      a,b = np.histogram(Rates, bins)
      ff = abs(np.fft.fft(Rates- np.mean(Rates)))**2
      Fs = 1./(1*0.001)
      freq2 = np.fft.fftfreq(len(bins))[0:len(bins/2)+1]
      freq = np.linspace(0,Fs/2,len(ff)/2+1)
      px = ff[0:len(ff)/2+1]
      max_px = np.max(px[1:])
      idx = px == max_px
      corr_freq = freq[pl.find(idx)]
      new_px = px
      max_pow = new_px[pl.find(idx)]
      return new_px,freq,corr_freq[0],freq2, max_pow


def spec_entropy(Rates,time_range=[],bin_w = 5.,freq_range = []):
	'''Function to calculate the spectral entropy'''

        power,freq,dfreq,dummy,dummy = mypsd(Rates,time_range,bin_w = bin_w)
        if freq_range != []:
                power = power[(freq>freq_range[0]) & (freq < freq_range[1])]
                freq = freq[(freq>freq_range[0]) & (freq < freq_range[1])]
        k = len(freq)
        power = power/sum(power)
        sum_power = 0
        for ii in range(k):
                sum_power += (power[ii]*np.log(power[ii]))
        spec_ent = -(sum_power/np.log(k))
        return spec_ent,dfreq



def amp_freq(timeSeries):
	freq = spfft.rfftfreq(len(timeSeries),1)
	fft = spfft.rfft(timeSeries)
        amp = np.abs(fft)
        indmax = np.where(amp==np.max(amp))[0][0]
        oi = (np.sum(amp[indmax-1:indmax+1]))/np.sum(amp)

	indmax = np.where(amp==np.max(amp))[0][0]
	dfreq = freq[indmax]*1000         # Sampling frequency = 1/1ms = 1/0.001= 1000
	
	return oi,dfreq

def calcPDFeatures(A,B):
	PDFeatures = np.zeros((1,7))
#	PDFeaturesSTP = np.zeros((1,7))

	delay = 1.0	
	#Calculate Rates for SWA and Control
	ipctx1=dict()
	#Calculate Rates for SWA and lesion(dopamine depletion)
	Flags = []
	Flags.append("Trans")
	i = 0	
	ipctx1["ip"] = np.zeros((1,2001))
	SWADopDepRates = psGA.calcRates(Flags,delay,A,B,False,ipctx1)
#	SWADopDepRatesSTP = psGA.calcRates(Flags,2,A,B,False,ipctx1)

	dt = p.params["dt"]
	howmuch = (np.mean(SWADopDepRates['gpi'][100/dt:500/dt]) -np.mean(SWADopDepRates['gpi'][600/dt:1400/dt]))/np.mean(SWADopDepRates['gpi'][100/dt:500/dt])		# (Orig - Final)/Orig
	PDFeatures[i][0] = howmuch


#	howmuchSTP = (np.mean(SWADopDepRatesSTP['gpi'][100/dt:500/dt]) -np.mean(SWADopDepRatesSTP['gpi'][600/dt:1400/dt]))/np.mean(SWADopDepRatesSTP['gpi'][100/dt:500/dt])		# (Orig - Final)/Orig
#	PDFeaturesSTP[i][0] = howmuchSTP


	# New features
	'''
	Flags=[]
	Flags.append("Trans")
	ipctx1["ip"] = np.zeros((1,2001))
	scale=1
	TransRates = psGA.calcRates(Flags,delay,A*scale,B*scale,False,ipctx1)
	'''
	#We dont need to call it twice
	TransRates = SWADopDepRates.copy()
	#import pdb
	#pdb.set_trace()
#	TransRatesSTP = psGA.calcRates(Flags,2,A,B*scale,False,ipctx1)
	#Add noise to remove the DC peak at zero. SE will be calculated correctly then
	# However adding noise at the input stage disturbs frmation of oscillation

	# GPe mostly dont oscillate - hence contributes to increasing SE
	'''
		noise = np.random.uniform(-4,4,len(TransRates['gpi'][520/dt:1400/dt]))
		TransRates['gpi'][550/dt:1400/dt] += noise
		TransRatesSTP['gpi'][550/dt:1400/dt] += noise
		#oi,dfreq = amp_freq(TransRates['gpi'][100:]-np.mean(TransRates['gpi'][100:]))	
		se,dfreq = spec_entropy(TransRates['gpi'][550/dt:1400/dt],time_range=TransRates["ipctx"][550/dt:1400/dt]) 
		seSTP,dfreqSTP = spec_entropy(TransRatesSTP['gpi'][550/dt:1400/dt],time_range=TransRatesSTP["ipctx"][550/dt:1400/dt]) 
			
		#if np.sum(checkZeroRates) < 0.5*len(checkZeroRates):
#		PDFeatures[i][1] = np.var(TransRates['gpi'][1100:2400])/np.mean(TransRates['gpi'][1100:2400]) 
		PDFeatures[i][1] = se
		PDFeaturesSTP[i][1] = seSTP
	'''
	# To catch damped oscillations	
	noise = np.random.uniform(0,0.05,len(TransRates['gpi'][502/dt:1400/dt]))
	TransRates['stn'][502/dt:1400/dt] += noise
#	TransRatesSTP['stn'][502/dt:1400/dt] += noise

	se1,dfreq1 = spec_entropy(TransRates['stn'][502/dt:1400/dt],time_range=TransRates["ipctx"][502/dt:1400/dt]) 
	PDFeatures[i][2] = se1 
#	seSTP1,dfreqSTP1 = spec_entropy(TransRatesSTP['stn'][502/dt:1400/dt],time_range=TransRatesSTP["ipctx"][502/dt:1400/dt]) 
#	PDFeaturesSTP[i][2] = seSTP1 

	noise = np.random.uniform(0,0.05,len(TransRates['ta'][502/dt:1400/dt]))
	TransRates['ta'][502/dt:1400/dt] += noise
#	TransRatesSTP['ta'][502/dt:1400/dt] += noise

	se2,dfreq2 = spec_entropy(TransRates['ta'][502/dt:1400/dt],time_range=TransRates["ipctx"][502/dt:1400/dt])
	PDFeatures[i][3] = se2
	PDFeatures[i][4] = dfreq2

#	seSTP2,dfreqSTP2 = spec_entropy(TransRatesSTP['ta'][502/dt:1400/dt],time_range=TransRatesSTP["ipctx"][502/dt:1400/dt])
#	PDFeaturesSTP[i][3] = seSTP2
#	PDFeaturesSTP[i][4] = dfreqSTP2


	noise = np.random.uniform(0,0.05,len(TransRates['ti'][502/dt:1400/dt]))
	TransRates['ti'][502/dt:1400/dt] += noise
#	TransRatesSTP['ti'][502/dt:1400/dt] += noise

	se3,dfreq3 = spec_entropy(TransRates['ti'][502/dt:1400/dt],time_range=TransRates["ipctx"][502/dt:1400/dt])

	PDFeatures[i][5] = se3
	PDFeatures[i][6] = dfreq3
	 
#	seSTP3,dfreqSTP3 = spec_entropy(TransRatesSTP['ti'][502/dt:1400/dt],time_range=TransRatesSTP["ipctx"][502/dt:1400/dt])
#
#	PDFeaturesSTP[i][5] = seSTP3
#	PDFeaturesSTP[i][6] = dfreqSTP3

	#return PDFeatures, PDFeaturesSTP,SWADopDepRates,TransRates
	return PDFeatures #, PDFeaturesSTP

def findMus(cluster,temp1,index):
	el = dict()
        if len(cluster) > 0:
                muGS = np.mean(temp1[cluster,0])
                muSE = np.mean(temp1[cluster,1])
        else:
                muGS = np.mean(temp1[0])
                muSE = np.mean(temp1[1])

	el["Ids"] = cluster
	el["index"] = index
	if muGS >0.9 and muSE >0.9 :
		el["label"] = "RobustHealthy"
		el["muGpiSup"] = muGS
		el["muSE"] = muSE
	elif muGS < 0.1 and muSE >0.9:
		el["label"] = "AR"
		el["muGpiSup"] = muGS
		el["muSE"] = muSE

	elif muGS > 0.7 and muSE < 0.8 and muSE > 0.5:
		el["label"] = "TDMild"
		el["muGpiSup"] = muGS
		el["muSE"] = muSE

	elif muGS > 0.7 and muSE < 0.5:
		el["label"] = "TDSevere"
		el["muGpiSup"] = muGS
		el["muSE"] = muSE

	elif muGS < 0.5 and muSE < 0.5:
		el["label"] = "RobustPathological" 
		el["muGpiSup"] = muGS
		el["muSE"] = muSE
	elif muGS < -0.1:
		el["label"] = "SevereAkinesia"
		el["muGpiSup"] = muGS
		el["muSE"] = muSE
	elif muGS > 0.8 and muSE > 0.8:
		el["label"] = "TowardsHealthy"
		el["muGpiSup"] = muGS
		el["muSE"] = muSE
	elif muGS < 0.1 and muGS > -0.1 and muSE > 0.5 and muSE < 0.85:
		el["label"] = "TowardsAkinesia"
		el["muGpiSup"] = muGS
		el["muSE"] = muSE
	
	else:
		el["label"] = "InBetween" 
		el["muGpiSup"] = muGS
		el["muSE"] = muSE
	return el



def findOnlyLabel(tempPDF):
        Allpoints=[]
        tempPDF1 = np.copy(tempPDF)
        for x in tempPDF1:
                el = dict()
                muGS = x[0]
                muSE = x[1]
                if muGS >0.9 and muSE >0.9 :
                        el["label"] = "RobustHealthy"
                        el["GpiSup"] = muGS
                        el["SE"] = muSE
                elif muGS < 0.1 and muSE >0.9:
                        el["label"] = "AR"
                        el["GpiSup"] = muGS
                        el["SE"] = muSE
                elif muGS > 0.7 and muSE < 0.8 and muSE > 0.5:
                        el["label"] = "TDMild"
                        el["GpiSup"] = muGS
                        el["SE"] = muSE
                elif muGS > 0.7 and muSE < 0.5:
                        el["label"] = "TDSevere"
                        el["GpiSup"] = muGS
                        el["SE"] = muSE
                elif muGS < 0.5 and muSE < 0.5:
                        el["label"] = "RobustPathological"
                        el["GpiSup"] = muGS
                        el["SE"] = muSE
                elif muGS < -0.1:
                        el["label"] = "SevereAkinesia"
                        el["GpiSup"] = muGS
                        el["SE"] = muSE

                elif muGS > 0.8 and muSE > 0.8:
                        el["label"] = "TowardsHealthy"
                        el["muGpiSup"] = muGS
                        el["muSE"] = muSE
                elif muGS < 0.1 and muGS > -0.1 and muSE > 0.5 and muSE < 0.85:
                        el["label"] = "TowardsAkinesia"
                        el["muGpiSup"] = muGS
                        el["muSE"] = muSE

                else:
                        el["label"] = "InBetween"
                        el["GpiSup"] = muGS
                        el["SE"] = muSE
                Allpoints.append(el)
        return Allpoints



def checkValidity(A,B):
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
        if np.mean(SWARates['d1'][1000:]) > 1.0 and np.mean(SWARates['d2'][1000:]) > 0.5 and np.mean(SWARates['fsi'][1000:]) > 0.1 and                  np.mean(SWARates['ta'][1000:]) > 0.5 and np.mean(SWARates['ti'][1000:]) > 0.5 and np.mean(SWARates['stn'][1000:]) > 0.1 and np.mean(SWARates['gpi'][1000:]) > 0.1:
                tests[6] = 1
        #if GpeSWA > np.mean(SWARates['d1']) and GpeSWA > np.mean(SWARates['d2']):
        tests[7] = 1

        #if abs(np.mean(SWARates['stn']) - np.mean(ActRates['stn']))<=10 and np.mean(SWARates['stn'])<=40: # Additonal constraint since STN firing rates were going very high. Also consult supplementary figure in Mallet 2008
        #if np.mean(SWARates['stn'])<=50 and np.mean(SWARates['d1']) > 1*np.mean(SWARates['d2']) and np.mean(ActRates['d1']) >1*np.mean(ActRates['d2']): # Additonal constraint since STN firing rates were going very high. Also consult supplementary figure in Mallet 2008
        #if np.mean(SWARates['stn'])<=50  and np.mean(SWARates['d1']) /np.mean(SWARates['d2']) > 1.2 and np.mean(ActRates['d1']) /np.mean(ActRates['d2']) >1.2: # and B[0][0] > 1.2*B[0][1] : # Additonal constraint since STN firing rates were going very high. Also consult supplementary figure in Mallet 2008 and (np.mean(ActRates['d1']) /np.mean(ActRates['d2'])) >  (np.mean(SWARates['d1']) /np.mean(SWARates['d2']))
        if np.mean(SWARates['stn'])<=50  and np.mean(SWARates['d1']) /np.mean(SWARates['d2']) > 1.5 and np.mean(ActRates['d1']) /np.mean(ActRates['d2']) >1.5: # and B[0][0] > 1.2*B[0][1] : # Additonal constraint since STN firing rates were going very high. Also consult supplementary figure in Mallet 2008 and (np.mean(ActRates['d1']) /np.mean(ActRates['d2'])) >  (np.mean(SWARates['d1']) /np.mean(SWARates['d2']))
      		tests[8] = 1
        #if np.round(np.mean(ActRates['ta'])) >= 5 and np.round(np.mean(ActRates['ta'])) <=25 and np.mean(SWARates['fsi'])*2.<np.mean(SWARates['d1']) : # Mean around 19
        if np.round(np.mean(ActRates['ta'])) >= 5 and np.round(np.mean(ActRates['ta'])) <=25 and np.mean(SWARates['fsi'])*1.5<np.mean(SWARates['d1']) : # Mean around 19
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

	return SWARates, ActRates, Rates


def checkPDF(TransRates):
	PDFeatures = np.zeros((1,8))
	dt = p.params["dt"]
	i = 0
	howmuch = (np.mean(TransRates['gpi'][100/dt:500/dt]) -np.mean(TransRates['gpi'][600/dt:1400/dt]))/np.mean(TransRates['gpi'][100/dt:500/dt])		# (Orig - Final)/Orig
	PDFeatures[i][0] = howmuch
	#start = 520
	noiseLim = 0.5
	timeStart = 510
	timeStart2 = 520
	#timeEnd = 1480	
	timeEnd1 = 1000 # This took a lot of time to tune, so there is a tradeoff between frequency resolution and time resolution,if frequency resolution increases, longer length of signal, temporal resolution deecreases, which means the peak in spectrum becomes wider and shifts towards higher frequencies, ultimately slipping outside the beta band. If shorter length of signal is sent (time_range < len(TransRates)), temporal resolution increases peak becomes thinner or more precise but then again shifts towards lower frequencies. So a compromise reached at 1000. For damped oscillations the signal length also needed to be decreased.  	
	timeEnd2 = 650   	
	timeEnd3 = 900   	
	# An alternatiev idea to SE could be , percentage of power contained in beta-band (periodogram)- because you require some value fot susceptibility to oscillations which varies between 0 and 1
	noise = np.random.uniform(-noiseLim,noiseLim,len(TransRates['gpi'][timeStart/dt:timeEnd1/dt]))
        time = np.arange(0,len(noise),1)
        wind = np.exp(-time/10.)
        noise1 = np.convolve(noise,wind,mode='same')

	#TransRates['stn'][timeStart/dt:timeEnd/dt] += noise
	#TransRates['stn'][timeStart/dt:timeEnd1/dt] += noise1[:len(noise)]

	#se1,dfreq1 = spec_entropy(TransRates['stn'][timeStart/dt:timeEnd/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd/dt],freq_range=[0.015,0.025]) 
	se11,dfreq11,maxFreq11,perMax11 = cPDF.spec_entropy(TransRates['stn'][timeStart/dt:timeEnd1/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd1/dt],freq_range=[0.00009,0.00017]) # 0.00010 = 10 Hz, since resolution = 10^-5 sec, 1000msec, and dt = 0.01 , to better the resolution of frequency captured by fft
	se12,dfreq12,maxFreq12,perMax12 = cPDF.spec_entropy(TransRates['stn'][timeStart/dt:timeEnd2/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd2/dt],freq_range=[0.00007,0.00017]) # 0.00010 = 10 Hz, since resolution = 10^-5 sec, 1000msec, and dt = 0.01 , to better the resolution of frequency captured by fft
	se13,dfreq13,maxFreq13,perMax13 = cPDF.spec_entropy(TransRates['stn'][timeStart/dt:timeEnd1/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd1/dt],freq_range=[0.00017,0.00036]) # 0.00010 = 10 Hz, since resolution = 10^-5 sec, 1000msec, and dt = 0.01 , to better the resolution of frequency captured by fft
	se14,dfreq14,maxFreq14,perMax14 = cPDF.spec_entropy(TransRates['stn'][timeStart/dt:timeEnd2/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd2/dt],freq_range=[0.00017,0.00036]) # 0.00010 = 10 Hz, since resolution = 10^-5 sec, 1000msec, and dt = 0.01 , to better the resolution of frequency captured by fft
	se15,dfreq15,maxFreq15,perMax15 = cPDF.spec_entropy(TransRates['stn'][timeStart2/dt:timeEnd3/dt],time_range=TransRates["ipctx"][timeStart2/dt:timeEnd3/dt],freq_range=[0.00017,0.00036]) # 0.00010 = 10 Hz, since resolution = 10^-5 sec, 1000msec, and dt = 0.01 , to better the resolution of frequency captured by fft

	#se1,dfreq1,maxFreq1,perMax = spec_entropy(TransRates['stn'][timeStart/dt:timeEnd/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd/dt],freq_range=[0.00008,0.00020]) # 0.00010 = 10 Hz, since resolution = 10^-5 sec, 1000msec, and dt = 0.01 , to better the resolution of frequency captured by fft
	ans = np.array([se11,se12,se13,se14,se15])
	indnan = np.where(np.isnan(ans)==True)
	ans[indnan] = 1
	print ans
	indmin = np.where(ans==np.min(ans))[0]
	print indmin
	#PDFeatures[i][2] = ans[indmin]
	PDFeatures[i][2] = np.min(ans)
#	if se11 < se12:
#		PDFeatures[i][2] = se11 
		#maxFreq[i][0] = maxFreq11
		#perPowerMaxFreq[i][0] = perMax11
#	else:
#		PDFeatures[i][2] = se12 
		#maxFreq[i][0] = maxFreq12
		#perPowerMaxFreq[i][0] = perMax12
	
	#TransRates['ta'][timeStart/dt:timeEnd/dt] += noise
	#TransRates['ta'][timeStart/dt:timeEnd1/dt] += noise1[:len(noise)]

	#se2,dfreq2 = spec_entropy(TransRates['ta'][timeStart/dt:timeEnd/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd/dt],freq_range=[0.015,0.025])
	se21,dfreq21,maxFreq21,perMax21 = cPDF.spec_entropy(TransRates['ta'][timeStart/dt:timeEnd1/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd1/dt],freq_range=[0.00009,0.00017])
	se22,dfreq22,maxFreq22,perMax22 = cPDF.spec_entropy(TransRates['ta'][timeStart/dt:timeEnd2/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd2/dt],freq_range=[0.00007,0.00017])
	se23,dfreq23,maxFreq23,perMax23 = cPDF.spec_entropy(TransRates['ta'][timeStart/dt:timeEnd1/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd1/dt],freq_range=[0.00017,0.00036])
	se24,dfreq24,maxFreq24,perMax24 = cPDF.spec_entropy(TransRates['ta'][timeStart/dt:timeEnd2/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd2/dt],freq_range=[0.00017,0.00036])
	se25,dfreq25,maxFreq25,perMax25 = cPDF.spec_entropy(TransRates['ta'][timeStart2/dt:timeEnd3/dt],time_range=TransRates["ipctx"][timeStart2/dt:timeEnd3/dt],freq_range=[0.00017,0.00036])
	#se2,dfreq2,maxFreq1,perMax = spec_entropy(TransRates['ta'][timeStart/dt:timeEnd/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd/dt],freq_range=[0.00008,0.00020])
	ans = np.array([se21,se22,se23,se24,se25])
	indnan = np.where(np.isnan(ans)==True)
	ans[indnan] = 1
	print ans
	indmin = np.where(ans==np.min(ans))[0]
	print indmin
	#PDFeatures[i][3] = ans[indmin]
	PDFeatures[i][3] = np.min(ans)

	'''
	if se21 < se22:
		PDFeatures[i][3] = se21
		PDFeatures[i][4] = dfreq21
		#maxFreq[i][1] = maxFreq21
		#perPowerMaxFreq[i][1] = perMax21
	else:
		PDFeatures[i][3] = se22
		PDFeatures[i][4] = dfreq22
		#maxFreq[i][1] = maxFreq22
		#perPowerMaxFreq[i][1] = perMax22
	'''	

	noise = np.random.uniform(-noiseLim,noiseLim,len(TransRates['ti'][timeStart/dt:timeEnd1/dt]))

	#TransRates['ti'][timeStart/dt:timeEnd/dt] += noise
	#TransRates['ti'][timeStart/dt:timeEnd1/dt] += noise1[:len(noise)]

	#se3,dfreq3 = spec_entropy(TransRates['ti'][timeStart/dt:timeEnd/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd/dt],freq_range=[0.015,0.025])
	se31,dfreq31,maxFreq31,perMax31 = cPDF.spec_entropy(TransRates['ti'][timeStart/dt:timeEnd1/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd1/dt],freq_range=[0.00009,0.00017])
	se32,dfreq32,maxFreq32,perMax32 = cPDF.spec_entropy(TransRates['ti'][timeStart/dt:timeEnd2/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd2/dt],freq_range=[0.00007,0.00017])
	se33,dfreq33,maxFreq33,perMax33 = cPDF.spec_entropy(TransRates['ti'][timeStart/dt:timeEnd1/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd1/dt],freq_range=[0.00017,0.00036])
	se34,dfreq34,maxFreq34,perMax34 = cPDF.spec_entropy(TransRates['ti'][timeStart/dt:timeEnd2/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd2/dt],freq_range=[0.00017,0.00036])
	se35,dfreq35,maxFreq35,perMax35 = cPDF.spec_entropy(TransRates['ti'][timeStart2/dt:timeEnd3/dt],time_range=TransRates["ipctx"][timeStart2/dt:timeEnd3/dt],freq_range=[0.00017,0.00036])
	#se3,dfreq3,maxFreq1,perMax = spec_entropy(TransRates['ti'][timeStart/dt:timeEnd/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd/dt],freq_range=[0.00008,0.00020])

	ans = np.array([se31,se32,se33,se34,se35])
	indnan = np.where(np.isnan(ans)==True)
	ans[indnan] = 1
	print ans
	indmin = np.where(ans==np.min(ans))[0]
	print indmin

	#PDFeatures[i][5] = ans[indmin]
	PDFeatures[i][5] = np.min(ans)

	'''
	if se31 < se32:
		PDFeatures[i][5] = se31
		PDFeatures[i][6] = dfreq31
		#maxFreq[i][1] = maxFreq31
		#perPowerMaxFreq[i][1] = perMax31
	else:
		PDFeatures[i][5] = se32
		PDFeatures[i][6] = dfreq32
		#maxFreq[i][1] = maxFreq32
		#perPowerMaxFreq[i][1] = perMax32
	'''

	se41,dfreq41,maxFreq41,perMax41 = cPDF.spec_entropy(TransRates['gpi'][timeStart/dt:timeEnd1/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd1/dt],freq_range=[0.00009,0.00017])
	se42,dfreq42,maxFreq42,perMax42 = cPDF.spec_entropy(TransRates['gpi'][timeStart/dt:timeEnd2/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd2/dt],freq_range=[0.00007,0.00017])
	se43,dfreq43,maxFreq43,perMax43 = cPDF.spec_entropy(TransRates['gpi'][timeStart/dt:timeEnd1/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd1/dt],freq_range=[0.00017,0.00036])
	se44,dfreq44,maxFreq44,perMax44 = cPDF.spec_entropy(TransRates['gpi'][timeStart/dt:timeEnd2/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd2/dt],freq_range=[0.00017,0.00036])
	se45,dfreq45,maxFreq45,perMax45 = cPDF.spec_entropy(TransRates['gpi'][timeStart2/dt:timeEnd3/dt],time_range=TransRates["ipctx"][timeStart2/dt:timeEnd3/dt],freq_range=[0.00017,0.00036])
	'''
	if se41 < se42:
		PDFeatures[i][7] = se41
	else:
		PDFeatures[i][7] = se42
	'''
	ans = np.array([se41,se42,se43,se44,se45])
	indnan = np.where(np.isnan(ans)==True)
	ans[indnan] = 1
	print ans
	indmin = np.where(ans==np.min(ans))[0]
	print indmin
	#PDFeatures[i][7] = ans[indmin]
	PDFeatures[i][7] = np.min(ans)

	return PDFeatures
