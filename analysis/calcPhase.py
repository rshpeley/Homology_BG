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

from scipy.optimize import curve_fit
import scipy.signal as sciSig

knownparams = p.params["known"]


storage_home = os.getcwd()+"../PD"  # Or healthy

path1 = storage_home+'/scripts/' # Path to store the results
path2 = storage_home+'/Dist/' # Path to jdf files
path3 = storage_home+'/jdf/'
path5 = storage_home+'/output/' # Path for output 

AllAs = pickle.load(open(path5+"AllAs.pickle","r"))
AllBs = pickle.load(open(path5+"AllBs.pickle","r"))

terms = ["d1ta","d1ti","d2ta","d2ti","fsita","fsiti","tad2","tata","tati","tid2","tita","titi","stnta","stnti","tistn","tastn","jc1","jc2","jfsictx","jstnctx"]
#terms1 = ["D1","D2","FSI","TA","TI","STN","GPi","Ctx"]
terms1 = ["d1","d2","fsi","ta","ti","stn","gpi","ipctx"]

# This is only to measure the phase difference for SWA 
def func(x,amp,base,pha):
	return amp*np.sin(2.*np.pi*0.002*x+pha)+base



def calcPhase():
	ipctx1 = dict()
	ipctx1["ip"] = np.zeros((1,2001))	
	dt = p.params["dt"]	

	randSamp = np.random.randint(0,len(AllAs),1000) 
	# Phases of three nuclei
	phasMus = dict()
	Flags = []
	Flags.append("SWA")	
	phasTAvsTA = []
	phasTIvsSTN = []
	phasTAvsSTN = []
	for i,samp in enumerate(randSamp):
		#temp = np.zeros((8,len(freqs),len(fftfreq)/10+1))
		A = AllAs[samp]
		B = AllBs[samp]
		Rates = psGA.calcRates(Flags,1.0,A,B,False,ipctx1,iptau=p.params["iptau"])
		
			

		# Filter the signals first at SWA frequency, or the phase difference between cortical input signal and TA,TI not calculated correctly due to multiple frequencies in the signal
		#B, A = sciSig.butter(2,np.array([0.0001,0.0005]),btype='low')
		B, A = sciSig.butter(2,np.array([0.00005]),btype='low')
		taFilt = sciSig.filtfilt(B, A, Rates['ta'])#, padlen=150)
		tiFilt = sciSig.filtfilt(B, A, Rates['ti'])#, padlen=150)
		stnFilt = sciSig.filtfilt(B, A, Rates['stn'])#, padlen=150)
		tG = np.arange(0,len(ipctx1["ip"][0])+dt,dt)
		
		fftfreq = np.fft.fftfreq(len(taFilt),d=dt)
		fftta = np.fft.rfft(taFilt-np.mean(taFilt))
		fftti = np.fft.rfft(tiFilt-np.mean(tiFilt))
		fftstn = np.fft.rfft(stnFilt-np.mean(stnFilt)) 
		
		fftip = np.fft.rfft(Rates['ipctx'] -np.mean(Rates['ipctx']))
		maxta = np.where(np.abs(fftta)==np.max(np.abs(fftta)))[0]
		maxti = np.where(np.abs(fftti)==np.max(np.abs(fftti)))[0]
		maxstn = np.where(np.abs(fftstn) == np.max(np.abs(fftstn)))[0]
		maxip = np.where(np.abs(fftip) == np.max(np.abs(fftip)))[0]

		print "maxta,maxti,maxstn,maxip",maxta,maxti,maxstn,maxip


		phasTAvsTA.append(np.angle(fftta[maxta]/fftta[maxta]))
		
		phasTIvsSTN.append(np.mean([np.angle(fftti[maxti]/fftstn[maxti]),np.angle(fftti[maxstn]/fftstn[maxstn])] ))
		phasTAvsSTN.append(np.mean([np.angle(fftta[maxta]/fftstn[maxta]),np.angle(fftta[maxstn]/fftstn[maxstn])] ))
	phasMus["ta_ta"] = phasTAvsTA
	phasMus["stn_ti"] = phasTIvsSTN
	phasMus["stn_ta"] = phasTAvsSTN
			 
	pickle.dump(phasMus,open(path5+"Phases.pickle","w"))
	
