
# Right now 2 features have been identified. 1) How much GPi is supressed when the striatal activity is high
# 2) How much is the network suceptible to oscillilations. 
# Every matrix is tagged with a 2-feature vector, showing how the matrix has scored on both the features.



import pickle
import numpy as np
import pylab as pl
import matplotlib as mpl
mpl.use('Agg')


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
storage_home = os.getcwd()+"../PD" # Or healthy 
from scipy.optimize import fsolve
ctx = 5.0
terms1 = ["d1","d2","fsi","ta","ti","stn","gpi","ipctx"]






def mypsd(Rates,time_range,bin_w = 5., nmax = 4000):

      bins = np.arange(0,len(time_range),1)
      #print bins
      a,b = np.histogram(Rates, bins)
      ff = (1./len(bins))*abs(np.fft.fft(Rates- np.mean(Rates)))**2
      Fs = 1./(1*0.001)
      freq2 = np.fft.fftfreq(len(bins))[0:len(bins/2)+1] # d= dt
      freq = np.fft.fftfreq(len(bins))[:len(ff)/2+1]
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
                power = power[(freq>=freq_range[0]) & (freq <= freq_range[1])]
                freq = freq[(freq>=freq_range[0]) & (freq <= freq_range[1])]
		maxFreq = freq[np.where(power==np.max(power))]*1000*100
		perMax = (np.max(power)/np.sum(power))*100
        k = len(freq)
        power = power/sum(power)
        sum_power = 0
        for ii in range(k):
                sum_power += (power[ii]*np.log(power[ii]))
        spec_ent = -(sum_power/np.log(k))
        return spec_ent,dfreq,maxFreq,perMax



def amp_freq(timeSeries):
	freq = spfft.rfftfreq(len(timeSeries),1)
	fft = spfft.rfft(timeSeries)
        amp = np.abs(fft)
        indmax = np.where(amp==np.max(amp))[0][0]
        oi = (np.sum(amp[indmax-1:indmax+1]))/np.sum(amp)

	indmax = np.where(amp==np.max(amp))[0][0]
	dfreq = freq[indmax]*1000         # Sampling frequency = 1/1ms = 1/0.001= 1000
	
	return oi,dfreq



def classification(params):
	pathname = storage_home+'/output/'
	print "in classification"
	knownparams = p.params["known"] 
	leak = -0.05
	Supressed=[]
	Sup_Ind=[]
	k = 0
	
	AllAs = pickle.load(open(pathname+"AllAs.pickle","r"))
	AllBs = pickle.load(open(pathname+"AllBs.pickle","r"))
	randSamps = np.random.randint(0,len(AllAs),1000) 

	# 0= GS, 1,2 = gpi SO,dominant frequency, 3,4 = TA SO,frequency, 5,6=TI SO, frequency
	PDFeatures = np.ones((len(randSamps),8))
	maxFreq = np.zeros((len(randSamps),3))
	perPowerMaxFreq = np.zeros((len(randSamps),3))
	ipAmp = params["ipAmp"]
	print ipAmp
	baseLineFR = []
	for i,j in enumerate(randSamps):
		print i
		A = AllAs[j]
		B = AllBs[j]

		delay = 1.0	
		ipctx1=dict()
		Flags = []
		Flags.append("Trans")
		
		ipctx1["ip"] = np.zeros((1,2001))
		SWADopDepRates = psGA.calcRates(Flags,delay,A,B,False,ipctx1,ipAmp,iptau=p.params["iptau"])

		# Record the ratio of Gpi rates between 500-1000/Rates before 500
		dt = p.params["dt"]
		howmuch = (np.mean(SWADopDepRates['gpi'][100/dt:500/dt]) -np.mean(SWADopDepRates['gpi'][600/dt:1400/dt]))/np.mean(SWADopDepRates['gpi'][100/dt:500/dt])		# (Orig - Final)/Orig
		PDFeatures[i][0] = howmuch




		# New features
		Flags=[]
		Flags.append("Trans")
		ipctx1["ip"] = np.zeros((1,2001))
		scale=1
		TransRates = psGA.calcRates(Flags,delay,A*scale,B*scale,False,ipctx1,ipAmp,iptau=p.params["iptau"])
		if ipAmp == 0:
			temp = dict()
			for x in terms1:
				temp[x] = np.mean(TransRates[x])
			baseLineFR.append(temp) 
		noiseLim = 0.5
		timeStart = 510
		timeStart2 = 520
		timeEnd1 = 1000   	
		timeEnd2 = 650   	
		timeEnd3 = 900   	
		# An alternatiev idea to SE could be , percentage of power contained in beta-band (periodogram)- because you require some value for susceptibility to oscillations which varies between 0 and 1
		noise = np.random.uniform(-noiseLim,noiseLim,len(TransRates['gpi'][timeStart/dt:timeEnd1/dt]))
		time = np.arange(0,len(noise),1)
		wind = np.exp(-time/10.)
		noise1 = np.convolve(noise,wind,mode='same')
		ind = np.where(noise1<0)
		noise1[ind] = 0.

		se11,dfreq11,maxFreq11,perMax11 =   spec_entropy(TransRates['stn'][timeStart/dt:timeEnd1/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd1/dt],freq_range=[0.00009,0.00017]) # 0.00010 = 10 Hz, since resolution = 10^-5 sec, 1000msec, and dt = 0.01 , to better the resolution of frequency captured by fft
		se12,dfreq12,maxFreq12,perMax12 =   spec_entropy(TransRates['stn'][timeStart/dt:timeEnd2/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd2/dt],freq_range=[0.00007,0.00017]) # 0.00010 = 10 Hz, since resolution = 10^-5 sec, 1000msec, and dt = 0.01 , to better the resolution of frequency captured by fft
		se13,dfreq13,maxFreq13,perMax13 =   spec_entropy(TransRates['stn'][timeStart/dt:timeEnd1/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd1/dt],freq_range=[0.00017,0.00036]) # 0.00010 = 10 Hz, since resolution = 10^-5 sec, 1000msec, and dt = 0.01 , to better the resolution of frequency captured by fft
		se14,dfreq14,maxFreq14,perMax14 =   spec_entropy(TransRates['stn'][timeStart/dt:timeEnd2/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd2/dt],freq_range=[0.00017,0.00036]) # 0.00010 = 10 Hz, since resolution = 10^-5 sec, 1000msec, and dt = 0.01 , to better the resolution of frequency captured by fft
		se15,dfreq15,maxFreq15,perMax15 =   spec_entropy(TransRates['stn'][timeStart2/dt:timeEnd3/dt],time_range=TransRates["ipctx"][timeStart2/dt:timeEnd3/dt],freq_range=[0.00017,0.00036]) # 0.00010 = 10 Hz, since resolution = 10^-5 sec, 1000msec, and dt = 0.01 , to better the resolution of frequency captured by fft
		se16,dfreq16,maxFreq16,perMax16 = spec_entropy(TransRates['stn'][timeEnd1/dt:(timeEnd1+480.)/dt],time_range=TransRates["ipctx"][timeEnd1/dt:(timeEnd1+500.)/dt],freq_range=[0.0001,0.00025]) # 0.00010 = 10 Hz, since resolution = 10^-5 sec, 1000msec, and dt = 0.01 , to better the resolution of frequency captured by fft , If the oscillation occurs at the end of the stimulation period

		ans = np.array([se11,se12,se13,se14,se15,se16])
		indnan = np.where(np.isnan(ans)==True)
		ans[indnan] = 1
		print ans
		indmin = np.where(ans==np.min(ans))[0]
		print indmin
		PDFeatures[i][2] = np.min(ans)

		se21,dfreq21,maxFreq21,perMax21 =   spec_entropy(TransRates['ta'][timeStart/dt:timeEnd1/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd1/dt],freq_range=[0.00009,0.00017])
		se22,dfreq22,maxFreq22,perMax22 =   spec_entropy(TransRates['ta'][timeStart/dt:timeEnd2/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd2/dt],freq_range=[0.00007,0.00017])
		se23,dfreq23,maxFreq23,perMax23 =   spec_entropy(TransRates['ta'][timeStart/dt:timeEnd1/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd1/dt],freq_range=[0.00017,0.00036])
		se24,dfreq24,maxFreq24,perMax24 =   spec_entropy(TransRates['ta'][timeStart/dt:timeEnd2/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd2/dt],freq_range=[0.00017,0.00036])
		se25,dfreq25,maxFreq25,perMax25 =   spec_entropy(TransRates['ta'][timeStart2/dt:timeEnd3/dt],time_range=TransRates["ipctx"][timeStart2/dt:timeEnd3/dt],freq_range=[0.00017,0.00036])
		se26,dfreq26,maxFreq26,perMax26 = spec_entropy(TransRates['ta'][timeEnd1/dt:(timeEnd1+480.)/dt],time_range=TransRates["ipctx"][timeEnd1/dt:(timeEnd1+500.)/dt],freq_range=[0.0001,0.00025]) # 0.00010 = 10 Hz, since resolution = 10^-5 sec, 1000msec, and dt = 0.01 , to better the resolution of frequency captured by fft , If the oscillation occurs at the end of the stimulation period
		ans = np.array([se21,se22,se23,se24,se25,se26])
		indnan = np.where(np.isnan(ans)==True)
		ans[indnan] = 1
		print ans
		indmin = np.where(ans==np.min(ans))[0]
		print indmin
		PDFeatures[i][3] = np.min(ans)

		noise = np.random.uniform(-noiseLim,noiseLim,len(TransRates['ti'][timeStart/dt:timeEnd1/dt]))

		se31,dfreq31,maxFreq31,perMax31 =   spec_entropy(TransRates['ti'][timeStart/dt:timeEnd1/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd1/dt],freq_range=[0.00009,0.00017])
		se32,dfreq32,maxFreq32,perMax32 =   spec_entropy(TransRates['ti'][timeStart/dt:timeEnd2/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd2/dt],freq_range=[0.00007,0.00017])
		se33,dfreq33,maxFreq33,perMax33 =   spec_entropy(TransRates['ti'][timeStart/dt:timeEnd1/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd1/dt],freq_range=[0.00017,0.00036])
		se34,dfreq34,maxFreq34,perMax34 =   spec_entropy(TransRates['ti'][timeStart/dt:timeEnd2/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd2/dt],freq_range=[0.00017,0.00036])
		se35,dfreq35,maxFreq35,perMax35 =   spec_entropy(TransRates['ti'][timeStart2/dt:timeEnd3/dt],time_range=TransRates["ipctx"][timeStart2/dt:timeEnd3/dt],freq_range=[0.00017,0.00036])
		se36,dfreq36,maxFreq36,perMax36 = spec_entropy(TransRates['ti'][timeEnd1/dt:(timeEnd1+480.)/dt],time_range=TransRates["ipctx"][timeEnd1/dt:(timeEnd1+500.)/dt],freq_range=[0.0001,0.00025])

		ans = np.array([se31,se32,se33,se34,se35,se36])
		indnan = np.where(np.isnan(ans)==True)
		ans[indnan] = 1
		print ans
		indmin = np.where(ans==np.min(ans))[0]
		print indmin

		PDFeatures[i][5] = np.min(ans)


		se41,dfreq41,maxFreq41,perMax41 =   spec_entropy(TransRates['gpi'][timeStart/dt:timeEnd1/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd1/dt],freq_range=[0.00009,0.00017])
		se42,dfreq42,maxFreq42,perMax42 =   spec_entropy(TransRates['gpi'][timeStart/dt:timeEnd2/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd2/dt],freq_range=[0.00007,0.00017])
		se43,dfreq43,maxFreq43,perMax43 =   spec_entropy(TransRates['gpi'][timeStart/dt:timeEnd1/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd1/dt],freq_range=[0.00017,0.00036])
		se44,dfreq44,maxFreq44,perMax44 =   spec_entropy(TransRates['gpi'][timeStart/dt:timeEnd2/dt],time_range=TransRates["ipctx"][timeStart/dt:timeEnd2/dt],freq_range=[0.00017,0.00036])
		se45,dfreq45,maxFreq45,perMax45 =   spec_entropy(TransRates['gpi'][timeStart2/dt:timeEnd3/dt],time_range=TransRates["ipctx"][timeStart2/dt:timeEnd3/dt],freq_range=[0.00017,0.00036])
		se46,dfreq46,maxFreq46,perMax46 = spec_entropy(TransRates['gpi'][timeEnd1/dt:(timeEnd1+480.)/dt],time_range=TransRates["ipctx"][timeEnd1/dt:(timeEnd1+500.)/dt],freq_range=[0.0001,0.00025])

		ans = np.array([se41,se42,se43,se44,se45,se46])
		indnan = np.where(np.isnan(ans)==True)
		ans[indnan] = 1
		print ans
		indmin = np.where(ans==np.min(ans))[0]
		print indmin
		PDFeatures[i][7] = np.min(ans)
 
		if i%250==0:
			pickle.dump(PDFeatures,open(pathname1+"PDFeaturesnew_"+str(ipAmp)+"_"+str(i)+".pickle","w"))			

	pickle.dump(baseLineFR,open(pathname1+"baseLineFR_"+str(ipAmp)+".pickle","w"))
	pickle.dump(PDFeatures,open(pathname1+"PDFeaturesnew_"+str(ipAmp)+".pickle","w"))
	pickle.dump(perPowerMaxFreq,open(pathname1+"PercentageInMaxPower_"+str(ipAmp)+".pickle","w"))
	pickle.dump(maxFreq,open(pathname1+"maxFreq_"+str(ipAmp)+".pickle","w"))

