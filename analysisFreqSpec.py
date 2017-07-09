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
#from sklearn.cluster import KMeans
import funcs as fs 

from scipy.optimize import curve_fit

knownparams = p.params["known"]


terms = ["d1ta","d1ti","d2ta","d2ti","fsita","fsiti","tad2","tata","tati","tid2","tita","titi","stnta","stnti","tistn","tastn","jc1","jc2","jfsictx","jstnctx"]
#terms1 = ["D1","D2","FSI","TA","TI","STN","GPi","Ctx"]
terms1 = ["d1","d2","fsi","ta","ti","stn","gpi","ipctx"]



storage_home = '/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams'

path1 = storage_home+'/scripts/' # Path to store the results
path2 = storage_home+'/Dist/' # Path to jdf files
path3 = storage_home+'/jdf/'
path5 = storage_home+'/output/' # Path for output 

AllAs = pickle.load(open(path5+"AllAs5.pickle","r"))
AllBs = pickle.load(open(path5+"AllBs5.pickle","r"))



def func(x,amp,base):
	return amp*np.sin(2.*np.pi*func.f*x)+base


func.f=0
def analFreqSpec(params1):
	clusPD = params1["clusPD"]
	clusNum = params1["i"]

	freqSpec = dict()	
	freqs = np.arange(0,400,5)
		
	ipctx1 = dict()
	ipctx1["ip"] = np.zeros((1,1001))	
	dt = p.params["dt"]	
	clusPDTemp = np.copy(clusPD)
	np.random.shuffle(clusPDTemp)
	randSamp = clusPDTemp[:30] # Take 30 random samples from each cluster
	freqSpec["samps"] = randSamp
	# Record the amplitude spectrums of all nuclei of all samples
	# Make a dummy array for frequencies 
	#tG = np.arange(0,len(ipctx1["ip"]),dt)
	#fftfreq = np.fft.rfftfreq(len(tG),0.01)
	fftfreq = pickle.load(open(path5+"rfftfreq.pickle","r")) # Had to calculate on local machine, since this version of numpy (1.7.2) doesnt have np.fft.rfftfreq
	inds = np.where(fftfreq*1000.<=400)[0]
	# Use np.fft.rfft for calculating only +ve side of spectrum. Also 0,400 is an interesting range. to plot use freq*1000
	ampSpec = []


	ampMus = np.zeros((len(randSamp),8,len(freqs))) # Record the mean amplitude of all nuclei for all frequencies
	Flags = []
	Flags.append("DiffFreq")	
	for i,samp in enumerate(randSamp):
		#temp = np.zeros((8,len(freqs),len(fftfreq)/10+1))
		temp = np.zeros((8,len(freqs),len(inds)))
		A = AllAs[samp]
		B = AllBs[samp]
		for j,fre in enumerate(freqs):
			Rates = psGA.calcRates(Flags,1.0,A,B,False,ipctx1,ipfreq=fre)
			for k,nuc in enumerate(terms1):
				ffttemp = np.fft.rfft(Rates[nuc]-np.mean(Rates[nuc]))
				fretemp = fftfreq[np.where(np.abs(ffttemp)==np.max(np.abs(ffttemp)))]
				func.f = fretemp
				tG = np.arange(0,len(ipctx1["ip"][0])+dt,dt)
				#print "np.shape(Rates[nuc])",np.shape(Rates[nuc])
				#print np.shape(tG) 
				#print "samp",samp
				#print "freq",fre
				#print "nuc",nuc
				if fre == 0 and nuc == 'ipctx': # Some error for this condition
					continue
				amp,b = curve_fit(func,tG[100./dt:],Rates[nuc][100./dt:])[0]	# Discard forst 100msecs to avoid transient effects
				#ampMus[i][k][j] = np.mean(Rates[nuc])
				ampMus[i][k][j] = np.abs(amp)
				temp[k][j][:] = (np.abs(ffttemp)/np.max(np.abs(ffttemp)))[inds]
		ampSpec.append(temp)

	freqSpec["ampSpec"] = ampSpec
	freqSpec["ampsMus"] = ampMus

	pickle.dump(freqSpec,open(path5+"freqSpec_"+str(clusNum)+".pickle","w"))
				
		
			

