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
import scipy.signal as sciSig

knownparams = p.params["known"]

ipctx1 = dict()
ipctx1["ip"] = np.zeros((1,2001))	
dt = p.params["dt"]	

Rates = pickle.load(open("TestRates.pickle","r"))


fig = pl.figure()
t1 = fig.add_subplot(211)
t1.plot(Rates['ipctx'],'k-',label='Ctx')
t1.plot(Rates['ti'],'m--',label='TI')
t1.plot(Rates['ta'],'c--',label='TA')
t1.plot(Rates['stn'],'y--',label='STN')
t1.legend()

B, A = sciSig.butter(2,[0.0001,0.0005],btype='low')
taFilt = sciSig.filtfilt(B, A, Rates['ta'])#, padlen=150)
tiFilt = sciSig.filtfilt(B, A, Rates['ti'])#, padlen=150)
stnFilt = sciSig.filtfilt(B, A, Rates['stn'])#, padlen=150)
t2 = fig.add_subplot(212)
t2.plot(taFilt,'c-',label='TA filtered')
t2.plot(tiFilt,'m-',label='TI filtered')
t2.plot(stnFilt,'y-',label='STN filtered')
t2.legend()
tG = np.arange(0,len(ipctx1["ip"][0])+dt,dt)
start = 0
#paramsTA = curve_fit(func,tG,Rates['ta'],p0=(np.std(Rates['ta']),np.mean(Rates['ta']),0))[0]
#paramsTA = sciSig.hilbert(Rates['ta'])
#paramsTA = sciSig.hilbert(taFilt[start/dt:])
#paramsSTN = curve_fit(func,tG,Rates['stn'],p0=(np.std(Rates['stn']),np.mean(Rates['stn']),0))[0]	
#paramsSTN = sciSig.hilbert(Rates['stn'])	
#paramsSTN = sciSig.hilbert(stnFilt[start/dt:])	
#paramsTI = curve_fit(func,tG,Rates['ti'],p0=(np.std(Rates['ti']),np.mean(Rates['ti']),0))[0]
#paramsTI = sciSig.hilbert(Rates['ti'])
#paramsTI = sciSig.hilbert(tiFilt[start/dt:])
#paramsCtx = sciSig.hilbert(Rates['ipctx'][start/dt:])
#TAvsTA = np.angle(np.inner( paramsTA, np.conj(paramsTA) ) / sqrt( np.inner(paramsTA,np.conj(paramsTA)) * np.inner(paramsTA,np.conj(paramsTA)) ))

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
TAvsTA = np.angle(np.angle(fftta[maxta]/fftta[maxta]))

#TIvsSTN = np.angle(np.inner( paramsTI, np.conj(paramsSTN) ) / sqrt( np.inner(paramsTI,np.conj(paramsTI)) * np.inner(paramsSTN,np.conj(paramsSTN)) ) )
TIvsSTN = np.mean([np.angle(fftti[maxti]/fftstn[maxti]),np.angle(fftti[maxstn]/fftstn[maxstn])] )
#TAvsSTN = np.angle(np.inner( paramsTA, np.conj(paramsSTN) ) / sqrt( np.inner(paramsTA,np.conj(paramsTA)) * np.inner(paramsSTN,np.conj(paramsSTN)) ) )
TAvsSTN = np.mean([np.angle(fftta[maxta]/fftstn[maxta]),np.angle(fftta[maxstn]/fftstn[maxstn])] )

print "TAvsTA",(180./np.pi)*TAvsTA
print "TIvsSTN",(180./np.pi)*TIvsSTN
print "TAvsSTN",(180./np.pi)*TAvsSTN

pl.show()
