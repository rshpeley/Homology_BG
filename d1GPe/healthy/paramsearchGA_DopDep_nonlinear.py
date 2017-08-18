
import matplotlib as mpl
mpl.use('Agg')

import pylab as pl
import numpy as np
import matplotlib.cm as cm
from pylab import *
# Since its a sparse matrix
import scipy.sparse.linalg as sp
import itertools
import marshal
import pickle
# Plotting the nullclines, how they vary with lambda-CTX
import logging
import gc
import knownUnknownParams as p
from scipy.integrate import odeint
import scipy.stats as sc_st
import os
def paramSearch(d):
	path = os.getcwd()+"/output/" 


	logging.basicConfig(level=logging.DEBUG)
	# Time trajectories
	leak = -0.1
	knownparams = p.params["known"] 

	unknownparams = p.params["unknown"] 

	ValidParams=[]
	InValidParams=[]
	numEpochs = 300
	lenList = len(unknownparams['d1ta'])
	lenUnknowns = len(unknownparams)	
	terms = ["d1ta","d2ta","fsita","fsiti","tata","tati","tastn","tita","titi","tistn","stnta","stnti","tid2","tad2","tid1","tad1","d1ti","d2ti","jc1","jc2","jfsictx","jstnctx"]		
	posInd = [6,9,18,19,20,21]
	negInd = np.array(list(set(np.arange(0,len(terms)))-set(posInd)))

	# Generate 300 uniform distributions for each unknown parameter
	
	elmt = np.zeros((numEpochs,lenUnknowns))
	for i in xrange(lenUnknowns):
		
		temp = np.random.uniform(np.min(unknownparams[terms[i]]),np.max(unknownparams[terms[i]]),numEpochs) 	
		elmt[:,i] = temp.T # rows are no of elements in GA and columns is number of uknown parameters
	

	newrandom = 0
	#Generate numEpochs combinations
	for epoch in xrange(4000):
		print "epoch",epoch
		if newrandom == 1:
			elmt = np.zeros((numEpochs,lenUnknowns))
			for i in xrange(lenUnknowns):
				temp = np.random.uniform(np.min(unknownparams[terms[i]]),np.max(unknownparams[terms[i]]),numEpochs) 	
				elmt[:,i] = temp.T # rows are no of elements in GA and columns is number of uknown parameters

		Grades = np.zeros(len(elmt))
		if epoch > 5: # Save the last comination so that we have a seed for next GA, stupid cluster has a limit of 24 hours
			pickle.dump(elmt,open(path+"LastComb.pickle","w"))

		for i,ind in enumerate(elmt):
			print ind
			d1ta = ind[0]
			d2ta = ind[1]
			fsita = ind[2]
			fsiti = ind[3]
			tata = ind[4]
			tati = ind[5]
			tastn = ind[6]
			tita = ind[7]
			titi = ind[8]
			tistn = ind[9]
			stnta = ind[10]
			stnti = ind[11]			
			tid2 = ind[12]
			tad2 = ind[13]	
			tid1 = ind[14]
			tad1 = ind[15]
			d1ti = ind[16]
			d2ti = ind[17]
			jc1 = ind[18]
			jc2 = ind[19]
			#jc2 = jc1 # Prevent trivial solutions
			jfsictx = ind[20]
			jstnctx = ind[21]
			
		
			A = np.matrix([[knownparams['d1d1'],knownparams['d1d2'],knownparams['d1fsi'],d1ta,d1ti,0.,0.],[knownparams['d2d1'],knownparams['d2d2'],knownparams['d2fsi'],d2ta,d2ti,0.,0.],[0.,0.,knownparams['fsifsi'],fsita,fsiti,0.,0.],[tad1,tad2,0.,tata,tati,tastn,0],[tid1,tid2,0.,tita,titi,tistn,0.],[0.,0.,0.,knownparams["stnta"],stnti,knownparams['stnstn'],0.],[knownparams['gpid1'],0.,0.,knownparams['gpita'],knownparams['gpiti'],knownparams['gpistn'],knownparams['gpigpi']]])
			#print A
			B = np.matrix([jc1,jc2,jfsictx,0,0,jstnctx,0])
			delay = 1.0	
			ipctx = dict()
			ipctx["ip"]=np.zeros((1,2001))	
		
			Flags = []
			Flags.append("SWA")

			SWARates = calcRates(Flags,delay,A,B,False,ipctx,iptau=p.params["iptau"])

			Flags = []
			Flags.append("Act")
			ActRates = calcRates(Flags,delay,A,B,False,ipctx,iptau=p.params["iptau"])

			tests = np.zeros(10)
			GpeSWA = np.mean(SWARates['ti']+SWARates['ta'])
			GpeAct = np.mean(ActRates['ti']+ActRates['ta'])
			
			if GpeSWA < GpeAct:
				tests[0] = 1
			# Check if STn is in phase with ctx (in-phase)
			if SWARates['stn_ctx'] > 0 :
				tests[1] = 1

			if np.round(np.mean(SWARates['ta'])) >= 0 and np.round(np.mean(SWARates['ta'])) <= 5 :
				tests[2] = 1
			
			# Check if TI is non-modulated
			if np.round(np.mean(SWARates['ti'])) >= 9.5 and np.round(np.mean(SWARates['ti'])) <= 45:
				tests[3] = 1 	
			
			if SWARates['taFF'] < 0.9 and SWARates['tiFF'] < 0.9:
				tests[4] = 1 
			
			if np.round(np.mean(ActRates['ti'])) >= 12 and np.round(np.mean(ActRates['ti'])) <=50:		
				tests[5] = 1
			# Sanity test , all rates are fairly above zero
			if np.mean(SWARates['d1'][1000:]) > 1.0 and np.mean(SWARates['d2'][100./p.params["dt"]:]) > 0.5 and np.mean(SWARates['fsi'][1000:]) > 0.1 and 			np.mean(SWARates['ta'][1000:]) > 0.5 and np.mean(SWARates['ti'][1000:]) > 0.5 and np.mean(SWARates['stn'][1000:]) > 0.5 and np.mean(SWARates['gpi'][100./p.params["dt"]:]) > 0.1:
				tests[6] = 1

			if np.mean(ActRates["ti"]) >np.mean(SWARates["ti"]) and np.mean(ActRates["ta"]) > np.mean(SWARates["ta"]):
				tests[7] = 1


		        if np.mean(SWARates['stn'])<=50 : 
      				tests[8] = 1
			if np.round(np.mean(ActRates['ta'])) >= 5 and np.round(np.mean(ActRates['ta'])) <=25  : 
				tests[9] = 1
			
			
			print "tests",tests	
			print ind
			print "SWA:d1,d2,fsi,ta,ti,stn,gpi,tha",np.mean(SWARates["d1"]),np.mean(SWARates["d2"]),np.mean(SWARates["fsi"]),np.mean(SWARates["ta"]),np.mean(SWARates["ti"]),np.mean(SWARates["stn"]),np.mean(SWARates["gpi"]),np.mean(SWARates["tha"])
			print "Act:d1,d2,fsi,ta,ti,stn,gpi,tha",np.mean(ActRates["d1"]),np.mean(ActRates["d2"]),np.mean(ActRates["fsi"]),np.mean(ActRates["ta"]),np.mean(ActRates["ti"]),np.mean(ActRates["stn"]),np.mean(ActRates["gpi"]),np.mean(ActRates["tha"])
			Grades[i] = np.sum(tests)	
				
		# Select the ones with th5 higest grades. Select here with 7,8 grades, but only 8 in findDist. To enable more parameter combs.
		cutoff = 9
		sorted_list = sorted(zip(Grades,elmt),reverse=True,key=lambda x:x[0])
		nextgen = [sorted_list[i] for i in xrange(len(sorted_list)) if sorted_list[i][0]>cutoff]
		no9s = len([e[0] for e in nextgen if e[0] == 9.])
		no10s = len([e[0] for e in nextgen if e[0] == 10.])
		indices = np.array([i for i,e in enumerate(nextgen) if e[0] >9.]) 
		print indices
		print "no9s",no9s
		print "no10s",no10s

		#print "nextgen",nextgen
		if epoch%1 == 0 and len(nextgen) > 0:
			pickle.dump(nextgen,open(path+"Combinations_new"+str(epoch+0)+".pickle","w"))


		#Generate numEpochs combinations from these by random crossovers
		if len(nextgen) > 0:
			# Changed this to get a new set of combinations, since this grid has a limit of 24 hours 
			print "len(nextgen)",len(nextgen)
			howmany = 1	
			#howmany = 3	
			which1 = [np.random.randint(0,len(nextgen)) for i in xrange(numEpochs/2)]	 	
			which2 = [np.random.randint(0,len(nextgen)) for i in xrange(numEpochs/2)]	 	
			startsN = [negInd[np.random.randint(0,len(negInd),1)] for i in xrange(numEpochs/2)] # separately for negative and positive, you dont want to end up with +ve value in inhibitory connection
			startsP = [posInd[np.random.randint(0,len(posInd),1)] for i in xrange(numEpochs/2)] # separately for negative and positive, you dont want to end up with +ve value in inhibitory connection
			newSet = []

			if len(indices) > 2:
				for i in xrange(0,numEpochs/2,1):
					#print nextgen[which1[i]]
					temp1 = np.copy(nextgen[which1[i]][1])		
					temp2 = np.copy(nextgen[which2[i]][1])
					#print "temp1",temp1
					#print "temp2",temp2
					#print "starts[i]",starts[i]
					temp3N = np.copy(temp1[startsN[i]:startsN[i]+howmany])
					temp3P = np.copy(temp1[startsP[i]:startsP[i]+howmany])
					if startsN[i]+howmany < len(temp1):
					#print "temp3",temp3
						temp1[startsN[i]:startsN[i]+howmany] = np.copy(temp2[startsN[i]:startsN[i]+howmany])
						temp2[startsN[i]:startsN[i]+howmany] = temp3N
					if startsP[i]+howmany < len(temp1):
					#print "temp3",temp3
						temp1[startsP[i]:startsP[i]+howmany] = np.copy(temp2[startsP[i]:startsP[i]+howmany])
						temp2[startsP[i]:startsP[i]+howmany] = temp3P	

					#print "newtemp1",temp1
					#print "newtemp2",temp2
					newSet.append(temp1)
					newSet.append(temp2)
				randNo = 100
				for x in xrange(randNo):
					mutset = np.random.randint(0,len(newSet))
					mutbit = np.random.randint(0,lenUnknowns)
					newSet[mutset][mutbit] = np.random.uniform(np.min(unknownparams[terms[mutbit]]),np.max(unknownparams[terms[mutbit]]),1)


			else:
				morerand = 100
				newSet.append(nextgen[0][1])
				for i in xrange(morerand):
					mutbit = np.random.randint(0,lenUnknowns)
					temp = np.copy(nextgen[0][1])
					temp[mutbit] = np.random.uniform(np.min(unknownparams[terms[mutbit]]),np.max(unknownparams[terms[mutbit]]),1)
					newSet.append(temp)
	 
			prev = np.copy(elmt)
			elmt = np.copy(newSet)
		
			gc.collect()
			newrandom=0
		else:
			newrandom = 1
			continue
	
	return elmt												
tempipc = np.zeros((1,1001))
dtt = p.params["dt"]
tGt = np.arange(0.0,len(tempipc[0])+dtt,dtt)
timet = np.arange(0,len(tGt),1)

ipFC = np.zeros((len(timet)))
noise = np.random.uniform(low=-0.1,high=0.1,size=len(ipFC))
#tau = 5.
#tau = 20.
tau = 10.
#expKern = np.exp(-np.arange(0,50,dtt)/tau)
expKern = np.exp(-np.arange(0,100,dtt)/tau)
LFPnoise = np.convolve(noise,expKern,mode='same')
ipFC = np.copy(LFPnoise)+1.

ipFC[np.where(ipFC<0)]=0.0

#--------------------------------------------------------
timeD = np.arange(0,1+dtt,dtt)
ipD = np.zeros(len(timeD))
ipD[20:80] = 99999999


#--------------------------------------------------------

def S(x,theta,Qmax):
	#sigma=4.
	sigma=4.
	funcVal = Qmax*(1./(1.+np.exp(-(x-theta)/sigma)))
	return funcVal 	

def getipFC(tG):
        dt = p.params["dt"]
        return ipFC[tG/dt]

def ipDiffFreqs(f,t):
        ipctx =2*np.sin(2.0*np.pi*f*t) +2.0
        return ipctx
def ipTFreq(t,ipamp):
	return ipamp


def ipSWA(t):

	ipctx =2*np.sin(2.0*np.pi*0.002*t) +2.0
	return ipctx

def ipT(t,ipamp):
	if t >=500 and t <=1500:
		x = ipamp
	else:
		x = 0.0
	return x
def ipAct(t):
	ipctx = 2*np.sin(2.0*np.pi*0.02*t)+2.5
	return ipctx
def getipD(tG):
	dt = p.params["dt"]
	return ipD[tG/dt]

def func(Rates,timeGrid,knownparams,unknownparams,time,Flag,ipamp,ipfreq,iptau):
	leak1 = -1
	tau = iptau
	ip = 0
	if Flag == "SWA":
		ip = ipSWA(timeGrid)
	if Flag == "Trans":
		ip = ipT(timeGrid,ipamp)
	if Flag == "TransFreq":
		ip = ipTFreq(timeGrid,ipamp)

	if Flag == "Act":
		ip = ipAct(timeGrid)
        if Flag == 'funcConn':
                #ip = ipFC[timeGrid]
                ip = getipFC(timeGrid)
        if Flag == 'DiffFreq':
                ip = ipDiffFreqs(ipfreq,timeGrid)

	if Flag == "delta":
		ip = getipD(timeGrid)	

#	print "ip",ip
	z1 = Rates[0]*knownparams['d1d1'] + Rates[1]*knownparams['d1d2'] + knownparams['d1fsi']*Rates[2] + unknownparams['jc1']*ip + Rates[3]*unknownparams['d1ta']+ Rates[4]*unknownparams['d1ti']
	d1bydt = (1./tau)*leak1*Rates[0] +(1./tau)* S(z1,0.1,65)

	z2 = Rates[1]*knownparams['d2d2'] + Rates[0]*knownparams['d2d1'] + knownparams['d2fsi']*Rates[2] +unknownparams['jc2']*ip + Rates[3]*unknownparams['d2ta']+Rates[4]*unknownparams['d2ti']	
	d2bydt = (1./tau)*leak1*Rates[1]+(1./tau)* S(z2,0.1,65)

	z3 = Rates[3]*unknownparams['fsita'] + Rates[4]*unknownparams['fsiti']+unknownparams['jfsictx']*ip + Rates[2]*knownparams['fsifsi']	
	dfsibydt = (1./tau)*leak1*Rates[2]+(1./tau)*S(z3,0.1,80)

	z4 = Rates[3]*unknownparams['tata'] + Rates[4]*unknownparams['tati'] + Rates[5]*unknownparams['tastn'] + Rates[1]*unknownparams['tad2']+ Rates[0]*unknownparams['tad1'] 
	dtabydt = (1./tau)*leak1*Rates[3]+(1./tau)* S(z4,0.4,75)

	z5 = Rates[3]*unknownparams['tita'] + Rates[4]*unknownparams['titi'] + Rates[5]*unknownparams['tistn'] + Rates[1]*unknownparams['tid2']+ Rates[0]*unknownparams['tid1']
	dtibydt = (1./tau)*leak1*Rates[4]+(1./tau)*S(z5,0.4,125)

 	z6 = unknownparams['jstnctx']*ip + Rates[4]*unknownparams['stnti'] + Rates[3]*unknownparams['stnta'] +Rates[5]*knownparams['stnstn']
	dstnbydt = (1./tau)*leak1*Rates[5]+(1./tau)*S(z6,0.4,500)

	z7 = Rates[0]*knownparams['gpid1']+Rates[3]*knownparams['gpita']+ Rates[4]*knownparams['gpiti'] + Rates[5]*knownparams['gpistn'] + Rates[6]*knownparams['gpigpi']
	dgpibydt = (1./tau)*leak1*Rates[6]+(1./tau)*S(z7,0.1,250) 

        z8 = 1
        dthabydt = (1./tau)*leak1*Rates[7]+(1./tau)*S(z8,3,10)


	return [d1bydt,d2bydt,dfsibydt,dtabydt,dtibydt,dstnbydt,dgpibydt,dthabydt]

def calcRates(Flags,d,A,B,fourier,ipctx1,ipamp=5.,ipfreq=2.,iptau=1):
	leak = -0.15
	if d == 1:
		knownparams = p.params["known"].copy() # Normal 
	if d == 2:
		knownparams = p.params["known"].copy() # Short term plasticity
		knownparams['gpid1'] = knownparams['gpid1']*2
	
	unknownparams = dict()
	unknownparams['d1ta'] = A[0,3]
	unknownparams['d2ta'] = A[1,3]
	unknownparams['fsita'] = A[2,3]
	unknownparams['fsiti'] = A[2,4]
	unknownparams['tata'] = A[3,3]
	unknownparams['tati'] = A[3,4]
	unknownparams['titi'] = A[4,4]
	unknownparams['tita'] = A[4,3]
	unknownparams['stnti'] = A[5,4]
	unknownparams['stnta'] = A[5,3]
	unknownparams['tastn'] = A[3,5]
	unknownparams['tistn'] = A[4,5]
	unknownparams['tid2'] = A[4,1]
	unknownparams['tad2'] = A[3,1]
	unknownparams['tid1'] = A[4,0]
	unknownparams['tad1'] = A[3,0]

	unknownparams['d1ti'] = A[0,4]
	unknownparams['d2ti'] = A[1,4]
	unknownparams['jc1'] = B[0,0]
	unknownparams['jc2'] = B[0,1]
	unknownparams['jfsictx'] = B[0,2]
	unknownparams['jstnctx'] = B[0,5]


	# Leak isnt considered for d1,d2,fsi,stn,gpi here, so include for them in equations!!!

	dt = p.params["dt"]
        timeGrid = np.arange(0.0,len(ipctx1["ip"][0])+dt,dt)
	time = np.arange(0,len(timeGrid),1)
	ipctx = np.zeros((len(time)))
	frate = 0.0
	for which in Flags:
			# Simulating dopamine depletion
		
		if which == 'SWA':
			ipctx[:] =2*np.sin(2.0*np.pi*0.002*timeGrid) +2.0
			frate = 250.	
			print "SWA"
		if which == 'Both':
			continue
		if which == 'STRint':
			# By default STR is non-PD, so set STN to PD
			knownparams['jstnctx'] = 5.5
		if which == 'Trans':
			#ipctx[0] = 5.0	
			ipctx[500/dt:1500/dt]= ipamp	
         	if which == "NoisyTrans":
                        noise = np.random.uniform(low=-0.1,high=0.1,size=len(ipctx[500/dt:1500/dt]))
                        tau = 5.
                        expKern = np.exp(-np.arange(0,30,1)/tau)
                        LFPnoise = np.convolve(noise,expKern,mode='same')
                        ipctx[500/dt:1500/dt] = 5.0+noise
	
		if which == 'funcConn':
                       # Apparently the method to generate pink noise is:, amplitude spectrum varies as 1/f
                        freq = np.fft.fftfreq(len(ipctx))
			mult= 1./len(freq)
                        ampSpec = (mult*np.random.uniform(0,2,len(freq)))/freq
                        timeSig = (1/mult)*np.fft.ifft(ampSpec[1:])
                        ipctx[1:] = np.abs(timeSig)

		if which == 'STNint':
			# Set STR to PD
			knownparams['jc1'] = 7.5 # LTP breaks down in D1
			knownparams['jc2'] = 8.5 # LTD breaks down in D2

		if which == 'Act':
			ipctx[:] = 2*np.sin(2.0*np.pi*0.02*timeGrid)+2.5
			frate = 100.
                if which == "DiffFreq":
                        ipctx[:] = 2*np.sin(2.0*np.pi*(ipfreq/1000.)*timeGrid)+2.0

	# These values were tuned such that no firing rates go to zero
	lamd1 = 0
	lamd2 = 0
	lamfsi =0
	lamgpi =0 #All these nuclei have tonic firing    # initial conditions - fixed point
	lamstn =0
	lamta = 0 
	lamti = 0#lamta = 0.0 

	initRates = np.ones((8))*1
	print "which",which
	finalRates = odeint(func,initRates,timeGrid,hmax=0.1,args = (knownparams,unknownparams,time,which,ipamp,ipfreq/1000.,iptau))	
	
	#print "finalRates",finalRates
	Rates = dict()
	Rates['d1'] = finalRates[:,0]
	Rates['d2'] = finalRates[:,1]
	Rates['fsi'] = finalRates[:,2]
	Rates['ta'] = finalRates[:,3]
	Rates['ti'] = finalRates[:,4]
	Rates['stn'] = finalRates[:,5]
	Rates['gpi'] = finalRates[:,6]	
	Rates['tha'] = finalRates[:,7]	
	Rates['A'] = A
	if which == 'funcConn':
		Rates['ipctx'] = ipFC
	else:
		Rates['ipctx'] = ipctx
	#print "Rates",Rates
	# Phase difference found by cross correlogram
	# Ignore first 100 msecs/timesteps to ignore transients
	#print np.shape(ipctx[1000:])
	#print np.shape(finalRates[:,5][1000:])
	if which != 'funcConn':
		Rates['stn_ctx'] = np.corrcoef(ipctx[1000/dt:],finalRates[:,5][1000/dt:])[0][1]	
		Rates['stn_ta'] = np.corrcoef(finalRates[:,3][1000/dt:],finalRates[:,5][1000/dt:])[0][1]	
		Rates['stn_ti'] = np.corrcoef(finalRates[:,4][1000/dt:],finalRates[:,5][1000/dt:])[0][1]	
		Rates['taFF'] = np.var(Rates['ta'][1000/dt:])/np.mean(Rates['ta'][1000/dt:])
		Rates['tiFF'] = np.var(Rates['ti'][1000/dt:])/np.mean(Rates['ti'][1000/dt:])	
	
	
	'''
	if fourier == True:
		fft_ti = np.fft.fft(Rates['ti'][:len(Rates['ipctx'][0])])
		fft_ta = np.fft.fft(Rates['ta'][:len(Rates['ipctx'][0])])
		fft_stn = np.fft.fft(Rates['stn'][:len(Rates['ipctx'][0])])
		fft_ctx = np.fft.fft(Rates['ipctx'][0])
		#dt = Rates['ipctx'][0][1]-Rates['ipctx'][0][0]
			
		freqs = np.fft.fftfreq(len(Rates['ipctx'][0]))
		indc = np.argmax(np.abs(fft_ctx[1:]))
		inds = np.argmax(np.abs(fft_stn[1:-1]))
		freq = freqs[indc]
		freq_hz = abs(freq*frate)
		# Return the phase at the index of cortical input	
		Rates['fft_ti_c'] = np.angle(fft_ti[indc])
		Rates['fft_ta_c'] = np.angle(fft_ta[indc])
		Rates['fft_stn_c'] = np.angle(fft_stn[indc])
		Rates['fft_ctx_c'] = np.angle(fft_ctx[indc])

		Rates['fft_ti_s'] = np.angle(fft_ti[inds])
		Rates['fft_ta_s'] = np.angle(fft_ta[inds])
		Rates['fft_stn_s'] = np.angle(fft_stn[inds])
		Rates['fft_ctx_s'] = np.angle(fft_ctx[inds])
	'''
	return Rates 





def plotTime(Rates,tit,Fourier): 
	from matplotlib import *
	time = np.arange(0,2000,1.)
	fig = pl.figure()
	pl.title('Time dependent behavior- '+tit)
	Ax = pl.gca()
	Ax.set_xticklabels(Ax.get_xticklabels(),visible=False)
	Ax.set_yticklabels(Ax.get_yticklabels(),visible=False)

	t1 = fig.add_subplot(8,1,1)
	t1.plot(time,Rates['d1'][:len(time)],'b-',label='D1',linewidth=1.5)
	setp(t1.get_xticklabels(),visible=False)
	for x in t1.get_yticklabels()[::2]:
		x.set_visible(False)
	t1.legend()

	t2 = fig.add_subplot(8,1,2,sharex=t1)
	setp(t2.get_xticklabels(),visible=False)
	t2.plot(time,Rates['d2'][:len(time)],'r-',label='D2',linewidth=1.5)
	#setp(t.get_yticklabels(),visible=False)
	for x in t2.get_yticklabels()[::2]:
		x.set_visible(False)

	t2.legend()

	t3 =fig.add_subplot(8,1,3,sharex=t1)
	t3.plot(time,Rates['fsi'][:len(time)],'g-',label='FSI',linewidth=1.5)
	setp(t3.get_xticklabels(),visible=False)
	for x in t3.get_yticklabels()[::2]:
		x.set_visible(False)

	#setp(t.get_yticklabels(),visible=False)
	t3.legend()

	t4 = fig.add_subplot(8,1,4,sharex=t1)
	t4.plot(time,Rates['ta'][:len(time)],'b-',label='GPe-TA',linewidth=1.5)
	setp(t4.get_xticklabels(),visible=False)
	for x in t4.get_yticklabels()[::2]:
		x.set_visible(False)

	#setp(t.get_yticklabels(),visible=False)
	t4.legend()

	t5 = fig.add_subplot(8,1,5,sharex=t1)
	t5.plot(time,Rates['ti'][:len(time)],'b-',label='GPe-TI',linewidth=1.5)
	setp(t5.get_xticklabels(),visible=False)
	for x in t5.get_yticklabels()[::2]:
		x.set_visible(False)

	#setp(t.get_yticklabels(),visible=False)
	t5.legend()

	t6 = fig.add_subplot(8,1,6,sharex=t1)
	t6.plot(time,Rates['stn'][:len(time)],'k-',label='STN',linewidth=1.5)
	setp(t6.get_xticklabels(),visible=False)
	for x in t6.get_yticklabels()[::2]:
		x.set_visible(False)

	#setp(t.get_yticklabels(),visible=False)
	t6.legend()

	#pl.savefig('Recovery_DBS_test.pdf')
	t7 = fig.add_subplot(8,1,7,sharex=t1)
	t7.plot(time,Rates['gpi'][:len(time)],'k-',label='Gpi',linewidth=1.5)
	setp(t7.get_xticklabels(),visible=False)
	for x in t7.get_yticklabels()[::2]:
		x.set_visible(False)

	#setp(t.get_yticklabels(),visible=False)
	t7.legend()

	t8 = fig.add_subplot(8,1,8,sharex=t1)
	t8.plot(time,Rates['ipctx'][0],'k-',label='Ctx',linewidth=1.5)
	#setp(t8.get_xticklabels(),visible=False)
	for x in t8.get_yticklabels()[::2]:
		x.set_visible(False)

	#setp(t.get_yticklabels(),visible=False)
	t8.legend()
	pl.savefig(tit+"TimeDep"+".pdf")	
	
	if Fourier == True:
		
		fig1 = pl.figure()
		# Fourier components analysis
		fft_ti = np.fft.fft(Rates['ti'][:len(Rates['ipctx'][0])])
		fft_ta = np.fft.fft(Rates['ta'][:len(Rates['ipctx'][0])])
		fft_stn = np.fft.fft(Rates['stn'][:len(Rates['ipctx'][0])])
		fft_ctx = np.fft.fft(Rates['ipctx'][0])

		freqs = np.fft.fftfreq(len(Rates['ipctx'][0]),Rates['ipctx'][0][1]-Rates['ipctx'][0][0])
		t9 = fig1.add_subplot(211)
		t9.set_title(tit+' Fourier Magnitude')
		t9.plot(freqs,np.abs(fft_ti),'b',label='TI')
		t9.plot(freqs,np.abs(fft_ta),'r',label='TA')
		t9.plot(freqs,np.abs(fft_stn),'g',label='STN')
		#t9.plot(freqs,np.abs(fft_ctx),'k',label='Ctx')
		t9.set_xlim(-0.15,0.15)
		t9.legend()

		t10 = fig1.add_subplot(212)
		t10.set_title(tit+' Fourier Phase') 
		t10.plot(freqs,np.angle(fft_ti),color='b',marker='.',label='TI')
		t10.plot(freqs,np.angle(fft_ta),color='r',marker='.',label='TA')
		t10.plot(freqs,np.angle(fft_stn),color='g',marker='.',label='STN')	
		#t10.plot(freqs,np.angle(fft_ctx),color='k',marker='^',label='Ctx')
#		t10.set_xlim(-0.15,0.15)
		t10.legend()
		pl.savefig("Fourier"+tit+".pdf")	
		




def diffD1D2(Flags):
	# To check if the difference between D1 and D2 amplifies downstream
	# First decide which model to use
	delay = 1
	if Flags == "Allsym":
		(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,params) = calcRates(Flags,delay)
		# D = Direct pathway, ID = Indirect pathway, HD = Hyperdirect pathway
		# Reducing a full recurrent matrix leads to postive and negative contributions in ID and HD instead of pure just positive contributions
		D = params['gpid1']
		print "Direct",D
		
		de1 = 1. + (params['d1d1']*params['fsifsi']) - (params['stnti']*params['stnstn'])
		Ex1 = (params['stnti']*params['d1d1']*params['fsifsi']*params['gpistn'])/de1  		
		Ex2 = (params['gpid1']*params['stnstn']*params['d1ta']*params['fsifsi'])/de1		
		IDpos = params['stnti']*params['tid2']*(params['d1ta']*params['gpid1']*params['tistn']+ Ex1 + Ex2)
		
		print "IDpos",IDpos

		Ex3 = (params['d1ta']+params['d1ta']*params['tata']+((params['stnta']*params['fsiti']*params['d1fsi'])/de1))
		IDneg = params['gpid1']+params['gpiti']*params['tid2']+params['gpid1']*params['stnstn']*params['tid2']*Ex3	
		print "IDneg",IDneg

		HDpos = (params['jstnctx']*params['gpid1']*params['stnstn']*params['fsiti']*params['d1fsi'])/de1
		print "HDpos",HDpos
		Ex4 = params['jstnctx']*params['d1d1']*params['fsifsi']*params['stnti']*params['gpistn']
		Ex5 = params['jstnctx']*params['gpid1']*params['stnstn']*params['d1ta']*params['fsifsi']
		HDneg = (Ex4 + Ex5)/de1	
		print "HDneg",HDneg

		d1fix = np.mean(d1[:-10])
		d2fix = np.mean(d2[:-10])
		DelMSN = d1fix - d2fix
		DelGpi = (D*d1fix) + ((IDpos - IDneg)*d2fix)

	return (DelMSN,DelGpi)			

def stabilityDiffData():
	delay=1
	Flags = []
	Flags.append('none')
	(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,params) = calcRates(Flags,delay)
	time = np.arange(0.0,2000.0,1.0)
	plotTime(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,time,'Normal Conditions',False)

	#(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B) = calcRates('gpeAsymm')
	#plotTime(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,time)

	#(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B) = calcRates('stngpe')
	#plotTime(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,time)
	Flags = []
	Flags.append('mallet')
	(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,params) = calcRates(Flags,delay)
	plotTime(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,time,'Bias from Mallet',False)

	Flags = []
	Flags.append('mastro')
	(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,params) = calcRates(Flags,delay)
	plotTime(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,time,'Norm- Mastro',False)

	Flags=[]
	Flags.append('Allsym')
	(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,params) = calcRates(Flags,delay)
	plotTime(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,time,'Norm - Allsym',False)

	Flags = []
	Flags.append('mastro')
	Flags.append('dopdep')
	(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,params) = calcRates(Flags,delay)
	plotTime(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,time,'Dop Dep-Mastro',False)

	Flags=[]
	Flags.append('dopdep')
	(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,params) = calcRates(Flags,delay)
	plotTime(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,time,'Dop Dep - Norm',False)

	Flags=[]
	Flags.append('mallet')
	Flags.append('dopdep')
	(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,params) = calcRates(Flags,delay)
	plotTime(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,time,'Dop Dep - Mallet',False)

	Flags=[]
	Flags.append('Allsym')
	Flags.append('dopdep')
	(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,params) = calcRates(Flags,delay)
	plotTime(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,time,'Dop Dep - Allsym',False)

	pl.show()

def checkPhaseRel():
	delay=1
	time = np.arange(0.0,2000.0,1.0)
	Flags = []
	Flags.append('none')
	(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,params) = calcRates(Flags,delay)
	plotTime(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,time,'Original tid2 in Normal Model',True)
	Flags.append('tid2')
	(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,params) = calcRates(Flags,delay)
	plotTime(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,time,'Depleted tid2 in Normal Model',True)


	Flags = []
	Flags.append('mallet')
	(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,params) = calcRates(Flags,delay)
	plotTime(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,time,'Original tid2 in Mallet',True)
	Flags.append('tid2')
	(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,params) = calcRates(Flags,delay)
	plotTime(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,time,'Depleted tid2 in Mallet',True)

	Flags = []
	Flags.append('mastro')
	(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,params) = calcRates(Flags,delay)
	plotTime(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,time,'Normal tid2 in Mastro',True)
	Flags.append('tid2')
	(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,params) = calcRates(Flags,delay)
	plotTime(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,time,'Depleted tid2 in Mastro',True)

	Flags = []
	Flags.append('Allsym')
	(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,params) = calcRates(Flags,delay)
	plotTime(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,time,'Normal tid2 in Allsym Model',True)
	Flags.append('tid2')
	(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,params) = calcRates(Flags,delay)
	plotTime(d1,d2,fsi,ti,ta,stn,gpi,ipctx,A,B,time,'Depleted tid2 in Allsym Model',True)

	pl.show()


delay = 1.0
#elmt = paramSearch(delay)
#pickle.dump(elmt,open("Combinations_DopDep.pickle","w"))
'''
elmt = pickle.load(open("Combinations.pickle","r"))
	knownparams = {
	 
	'd1d1' : -0.7,
#	'd1d1' : -0.0,
	'd1d2' : -0.73,
#	'd1d2' : -0.0,
	'd2d1' : -0.16,
#	'd2d1' : -0.0,
	'd2d2' : -0.97,
#	'd2d2' : -0.0,
	'd2fsi' : -0.22,
#	'd2fsi' : -0.0,
	'd1fsi' : -0.324,
#	'd1fsi' : -0.0,
	'gpid1' : -8.9,
	'gpiti' : -0.17,
	'tid2' : -3.27,
	'gpistn' : 1.98,
	'tad2' : 0.0,
	'stnstn' : 0.001,
	'gpita' : 0.0,
	'gpigpi' : 0.0,
	'd1ti' : 0.0,
	'd2ti' : 0.0,
	'fsifsi':-0.0012,
	'jc1' : 7.5,
	'jc2' : 6.5,
	'jstnctx' : 4.5, # Trying to tri-phasic response in GPi
	'jfsictx' : 1.8
	}

	unknownparams = {
	'stnta' : np.arange(-0.01,-0.65,-0.15),
	#'stnta' : np.arange(-0.01,-0.05,-0.05),
	'stnti' : np.arange(-0.01,-0.65,-0.15),
	#'stnti' : np.arange(-0.01,-0.05,-0.05),
	'tistn' : np.arange(0.01,0.82,0.2),
	#'tistn' : np.arange(0.01,0.05,0.05),
	'tastn' : np.arange(0.01,0.82,0.2),
	#'tastn' : np.arange(0.01,0.05,0.05),
	'tata': np.arange(-0.01,-0.65,-0.15),
	#'tata': np.arange(-0.01,-0.05,-0.05),
	'tati': np.arange(-0.01,-0.65,-0.15),
	#'tati': np.arange(-0.01,-0.05,-0.05),
	'tita': np.arange(-0.01,-0.65,-0.15),
	#'tita': np.arange(-0.01,-0.05,-0.05),
	'titi': np.arange(-0.01,-0.65,-0.15),
	#'titi': np.arange(-0.01,-0.05,-0.05),
	'd1ta' : np.arange(-0.01,-0.65,-0.15),
	#'d1ta' : np.arange(-0.01,-0.05,-0.05),
	'd2ta' : np.arange(-0.01,-0.65,-0.15),
	#'d2ta' : np.arange(-0.01,-0.05,-0.05),
	'fsita' : np.arange(-0.01,-0.65,-0.15),
	#'fsita' : np.arange(-0.01,-0.05,-0.05),
	'fsiti' : np.arange(-0.01,-0.65,-0.15),
	'tid2' : np.arange(0,-5,-1),
	'tad2' : np.arange(0,-5,-1)

	#'fsiti' : np.arange(-0.01,-0.05,-0.05),
	}

for i,ind in enumerate(elmt):
	d1ta = unknownparams['d1ta'][ind[0]]
	d2ta = unknownparams['d2ta'][ind[1]]
	fsita = unknownparams['fsita'][ind[2]]
	fsiti = unknownparams['fsiti'][ind[3]]
	tata = unknownparams['tata'][ind[4]]
	tita = unknownparams['tita'][ind[5]]
	tastn = unknownparams['tastn'][ind[6]]
	tita = unknownparams['tita'][ind[7]]
	titi = unknownparams['titi'][ind[8]]
	tistn = unknownparams['tistn'][ind[9]]
	stnta = unknownparams['stnta'][ind[10]]
	stnti = unknownparams['stnti'][ind[11]]			
	A = np.matrix([[leak+knownparams['d1d1'],knownparams['d1d2'],knownparams['d1fsi'],d1ta,knownparams['d1ti'],0.,0.],[knownparams['d2d1'],leak+knownparams['d2d2'],knownparams['d2fsi'],d2ta,knownparams['d2ti'],0.,0.],[0.,0.,leak+knownparams['fsifsi'],fsita,fsiti,0.,0.],[0,knownparams['tad2'],0.,leak+tata,tita,tastn,0],[0.,knownparams['tid2'],0.,tita,leak+titi,tistn,0.],[0.,0.,0.,stnta,stnti,leak+knownparams['stnstn'],0.],[knownparams['gpid1'],0.,0.,knownparams['gpita'],knownparams['gpiti'],knownparams['gpistn'],leak+knownparams['gpigpi']]])

Flags = []
Flags.append("SWA")
Flags.append("Ctrl")
Rates = calcRates(Flags,delay,A)
plotTime(Rates,"SWA - Ctrl",True) 

Flags = []
Flags.append("SWA")
Flags.append("DopDep")
Rates = calcRates(Flags,delay,I[0])
plotTime(Rates,"SWA - DopDep",True) 

Flags = []
Flags.append("Act")
Flags.append("Ctrl")
Rates = calcRates(Flags,delay,I[0])
plotTime(Rates,"Act - Ctrl",True) 

Flags = []
Flags.append("Act")
Flags.append("DopDep")
Rates = calcRates(Flags,delay,I[0])
plotTime(Rates,"Act - DopDep",True) 
'''





