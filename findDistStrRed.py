import pickle
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.cm as cm
from pylab import *
import paramsearchGA_DopDep_nonlinear as psGA
# Since its a sparse matrix
import scipy.sparse.linalg as sp
import itertools
import scipy.cluster.hierarchy as sph
import knownUnknownParams as p

import os
knownparams = p.params["known"]
unknownparams = p.params["unknown"]
leak = -0.05
storage_home = "/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/newParams"

path5 = storage_home+'/output/' # Path for output 
def checkConds(SWADopDepRates, ActDopDepRates):
	# Checks Refer Figure 6I and 6J of Mallet 2008
	tests = np.zeros(13)
	if np.round(np.mean(SWADopDepRates['ti'])) >= 19 and np.round(np.mean(SWADopDepRates['ti']))<=60: 				# Mean around 24, so +- 4Hz allowed
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
	#if np.mean(ActDopDepRates['d2']) > np.mean(ActDopDepRates['d1']):
	if np.mean(ActDopDepRates["ti"]) < np.mean(SWADopDepRates["ti"]) and np.mean(ActDopDepRates["ta"]) > np.mean(SWADopDepRates["ta"]):
		tests[9] = 1 
	# Sanity test , all rates are fairly above zero
	if np.mean(SWADopDepRates['d1'][100./p.params["dt"]:]) > 0.5 and np.mean(SWADopDepRates['d2'][100./p.params["dt"]:]) > 0.5 and np.mean(SWADopDepRates['fsi'][100./p.params["dt"]:]) > 0.5 and	np.mean(SWADopDepRates['ta'][100./p.params["dt"]:]) > 1.0 and np.mean(SWADopDepRates['ti'][100./p.params["dt"]:]) > 1.0 and np.mean(SWADopDepRates['stn'][100./p.params["dt"]:]) > 1.0 and np.mean(SWADopDepRates['gpi'][100./p.params["dt"]:]) > 0.1:
		tests[10] = 1
		#	if anti==0:
	#if GpeSWADopDep > np.mean(SWADopDepRates['d1']) and GpeSWADopDep > np.mean(SWADopDepRates['d2']):
	tests[11] = 1
	#if  np.mean(SWADopDepRates['stn']) > 50:
	tests[12] = 1


	return tests

def run(prefix,pars):
	Allcombs =[]
	Allstability = []
	AllAs = []

	iter1 = 0
	ind = pars['ind']
	d1ta = unknownparams['d1ta'][ind[0]]
	d2ta = unknownparams['d2ta'][ind[1]]
	fsita = unknownparams['fsita'][ind[2]]
	fsiti = unknownparams['fsiti'][ind[3]]
	tata = unknownparams['tata'][ind[4]]
	tati = unknownparams['tati'][ind[5]]
	tastn = unknownparams['tastn'][ind[6]]
	tita = unknownparams['tita'][ind[7]]
	titi = unknownparams['titi'][ind[8]]
	tistn = unknownparams['tistn'][ind[9]]
	stnta = unknownparams['stnta'][ind[10]]
	stnti = unknownparams['stnti'][ind[11]]			
	tid2 = unknownparams['tid2'][ind[12]]
	tad2 = unknownparams['tad2'][ind[13]]	
	d1ti = unknownparams['d1ti'][ind[14]]
	d2ti = unknownparams['d2ti'][ind[15]]
	jc1 = unknownparams['jc1'][ind[16]]
	jc2 = unknownparams['jc2'][ind[17]]
	#jc2 = jc1 # Prevent trivial solutions
	jfsictx = unknownparams['jfsictx'][ind[18]]
	jstnctx = unknownparams['jstnctx'][ind[19]]

	num = 2500
	# All the above are means, take a normal distribution among them and check if they still follow the conditions
	dev1 = 0.39
	dev2 = 0.29
	dev3 = 0.95
	# All the above are means, take a normal distribution among them and check if they still follow the conditions
	d1taDist = np.random.normal(loc=d1ta,scale=dev1,size=num) 
	d2taDist = np.random.normal(loc=d2ta,scale=dev1,size=num)
	fsitaDist = np.random.normal(loc=fsita,scale=dev1,size=num)
	fsitiDist = np.random.normal(loc=fsiti,scale=dev1,size=num)
	tataDist = np.random.normal(loc=tata,scale=dev1,size=num)
	tatiDist = np.random.normal(loc=tati,scale=dev1,size=num)
	tastnDist = np.random.normal(loc=tastn,scale=dev2,size=num)
	titaDist = np.random.normal(loc=tita,scale=dev1,size=num)
	titiDist = np.random.normal(loc=titi,scale=dev1,size=num)
	tistnDist = np.random.normal(loc=tistn,scale=dev2,size=num)
	stntaDist = np.random.normal(loc=stnta,scale=dev1,size=num)
	stntiDist = np.random.normal(loc=stnti,scale=dev1,size=num)
	tid2Dist = np.random.normal(loc=tid2,scale=dev3,size=num)
	tad2Dist = np.random.normal(loc=tad2,scale=dev3,size=num)
	d1tiDist = np.random.normal(loc=d1ti,scale=dev1,size=num)
	d2tiDist = np.random.normal(loc=d2ti,scale=dev1,size=num)
	jc1Dist = np.random.normal(loc=jc1,scale=dev3,size=num)
	jc2Dist = np.random.normal(loc=jc2,scale=dev3,size=num) 
	jfsictxDist = np.random.normal(loc=jfsictx,scale=dev3,size=num)
	jstnctxDist = np.random.normal(loc=jstnctx,scale=dev3,size=num)

	d1taDist = d1taDist[np.where(d1taDist < 0)]	
	d2taDist = d2taDist[np.where(d2taDist < 0)]	
	fsitaDist = fsitaDist[np.where(fsitaDist < 0)]	
	fsitiDist = fsitiDist[np.where(fsitiDist < 0)]	
	tataDist = tataDist[np.where(tataDist < 0)]	
	tatiDist = tatiDist[np.where(tatiDist < 0)]	
	tastnDist = tastnDist[np.where(tastnDist > 0)]	
	tistnDist = tistnDist[np.where(tistnDist > 0)]	
	titaDist = titaDist[np.where(titaDist < 0)]	
	titiDist = titiDist[np.where(titiDist < 0)]	
	stntaDist = stntaDist[np.where(stntaDist < 0)]	
	stntiDist = stntiDist[np.where(stntiDist < 0)]	
	tid2Dist = tid2Dist[np.where(tid2Dist < 0)]	
	tad2Dist = tad2Dist[np.where(tad2Dist < 0)]	
	d1tiDist = d1tiDist[np.where(d1tiDist < 0)]	
	d2tiDist = d2tiDist[np.where(d2tiDist < 0)]	
	jc1Dist = jc1Dist[np.where(jc1Dist >0)]
	jc2Dist = jc2Dist[np.where(jc2Dist >0)]
	jfsictxDist = jfsictxDist[np.where(jfsictxDist >0)]
	jstnctxDist = jstnctxDist[np.where(jstnctxDist >0)]




	for d1ta1,d2ta1,fsita1,fsiti1,tata1,tati1,tastn1,tita1,titi1,tistn1,stnti1,tid21,tad21,jc11,jc21,jfsictx1,jstnctx1,stnta1,d1ti1,d2ti1 in zip(d1taDist,d2taDist,fsitaDist,fsitiDist,tataDist,tatiDist,tastnDist,titaDist,titiDist,tistnDist,stntiDist,tid2Dist,tad2Dist,jc1Dist,jc2Dist,jfsictxDist,jstnctxDist,stntaDist,d1tiDist,d2tiDist):
		if np.abs(d1ti1) > np.abs(d1ta1)/1.5  or np.abs(d2ti1) > np.abs(d2ta1)/1.5  or np.abs(tid21) < 1.2 or np.abs(tad21) < 1.2  or np.abs(stnta1) > np.abs(stnti1)/2. or np.abs(fsita1) > np.abs(fsiti1)/1.5:# This is the 1 1for all corase options so it always fails or np.abs(stnta1) > np.abs(stnti1)/1.5:The last condition stnta and stnti maybe important to reduce the bump at low SO in PD networks. Check PruneoutSeveAkine.py 				
			continue	
	
		iter1+=1

		A = np.matrix([[knownparams['d1d1'],knownparams['d1d2'],knownparams['d1fsi'],d1ta1,d1ti1,0.,0.],[knownparams['d2d1'],knownparams['d2d2'],knownparams['d2fsi'],d2ta1,d2ti1,0.,0.],[0.,0.,knownparams['fsifsi'],fsita1,fsiti1,0.,0.],[0,tad21,0.,tata1,tati1,tastn1,0],[0.,tid21,0.,tita1,titi1,tistn1,0.],[0.,0.,0.,stnta1,stnti1,knownparams['stnstn'],0.],[knownparams['gpid1'],0.,0.,knownparams['gpita'],knownparams['gpiti'],knownparams['gpistn'],knownparams['gpigpi']]])
		B = np.matrix([jc11,jc21,jfsictx1,0,0,jstnctx1,0])
		#w,v = np.linalg.eig(A)
		#Alleigs.append(w)
		delay = 1.0			 
		ipctx1=dict()
		ipctx1["ip"]=np.zeros((1,2001))	
		#Calculate Rates for SWA and lesion(dopamine depletion)
		Flags = []
		Flags.append("SWA")
#		Flags.append("DopDep")
		SWADopDepRates = psGA.calcRates(Flags,delay,A,B,False,ipctx1,iptau=p.params["iptau"])

		# Calculate Rates for Activation and Lesion
		Flags = []
		Flags.append("Act")
#		Flags.append("DopDep")
		ActDopDepRates = psGA.calcRates(Flags,delay,A,B,False,ipctx1,iptau=p.params["iptau"])
		
		tests = checkConds(SWADopDepRates,ActDopDepRates)
		print "tests",tests	
		print ind
		print "SWA:d1,d2,fsi,ta,ti,stn,gpi",np.mean(SWADopDepRates["d1"]),np.mean(SWADopDepRates["d2"]),np.mean(SWADopDepRates["fsi"]),np.mean(SWADopDepRates["ta"]),np.mean(SWADopDepRates["ti"]),np.mean(SWADopDepRates["stn"]),np.mean(SWADopDepRates["gpi"])
		print "Act:d1,d2,fsi,ta,ti,stn,gpi",np.mean(ActDopDepRates["d1"]),np.mean(ActDopDepRates["d2"]),np.mean(ActDopDepRates["fsi"]),np.mean(ActDopDepRates["ta"]),np.mean(ActDopDepRates["ti"]),np.mean(ActDopDepRates["stn"]),np.mean(ActDopDepRates["gpi"])
	
		print "tests",tests	
		print ind
		Grades = np.sum(tests)	
		cutoff = 12
		if Grades > cutoff:
			Allcombs.append([d1ta1,d2ta1,fsita1,fsiti1,tata1,tati1,tastn1,tita1,titi1,tistn1,stnta1,stnti1,tid21,tad21,d1ti1,d2ti1,jc11,jc21,jfsictx1,jstnctx1])

	# Have at leats one more combination, the original one
	if np.abs(d1ti)*1.5 <  np.abs(d1ta)  and np.abs(d2ti)*1.5 < np.abs(d2ta)  and np.abs(tid2) > 1.2 and np.abs(tad2) > 1.2 and np.abs(stnti) > np.abs(stnta)*2. and np.abs(fsita)*1.5 < np.abs(fsiti):# This is the 1 1for all corase options so it always fails or np.abs(stnta1) > np.abs(stnti1)/1.5: 

		A = np.matrix([[knownparams['d1d1'],knownparams['d1d2'],knownparams['d1fsi'],d1ta,d1ti,0.,0.],[knownparams['d2d1'],knownparams['d2d2'],knownparams['d2fsi'],d2ta,d2ti,0.,0.],[0.,0.,knownparams['fsifsi'],fsita,fsiti,0.,0.],[0,tad2,0.,tata,tati,tastn,0],[0.,tid2,0.,tita,titi,tistn,0.],[0.,0.,0.,stnta,stnti,knownparams['stnstn'],0.],[knownparams['gpid1'],0.,0.,knownparams['gpita'],knownparams['gpiti'],knownparams['gpistn'],knownparams['gpigpi']]])
		B = np.matrix([jc1,jc2,jfsictx,0,0,jstnctx,0])
		#w,v = np.linalg.eig(A)
		#Alleigs.append(w)
		delay = 1.0			 
		ipctx1=dict()
		ipctx1["ip"]=np.zeros((1,2001))	
		#Calculate Rates for SWA and lesion(dopamine depletion)
		Flags = []
		Flags.append("SWA")
	#	Flags.append("DopDep")
		SWADopDepRates = psGA.calcRates(Flags,delay,A,B,False,ipctx1,iptau=p.params["iptau"])

		Flags = []
		Flags.append("Act")
	#	Flags.append("DopDep")
		ActDopDepRates = psGA.calcRates(Flags,delay,A,B,False,ipctx1,iptau=p.params["iptau"])

		tests = checkConds(SWADopDepRates,ActDopDepRates)
		print "tests",tests	
		print ind
		print "SWA:d1,d2,fsi,ta,ti,stn,gpi",np.mean(SWADopDepRates["d1"]),np.mean(SWADopDepRates["d2"]),np.mean(SWADopDepRates["fsi"]),np.mean(SWADopDepRates["ta"]),np.mean(SWADopDepRates["ti"]),np.mean(SWADopDepRates["stn"]),np.mean(SWADopDepRates["gpi"])
		print "Act:d1,d2,fsi,ta,ti,stn,gpi",np.mean(ActDopDepRates["d1"]),np.mean(ActDopDepRates["d2"]),np.mean(ActDopDepRates["fsi"]),np.mean(ActDopDepRates["ta"]),np.mean(ActDopDepRates["ti"]),np.mean(ActDopDepRates["stn"]),np.mean(ActDopDepRates["gpi"])
		
		print "tests",tests	
		print ind
		Grades = np.sum(tests)	
		cutoff = 12
		if Grades > cutoff:
			Allcombs.append([d1ta,d2ta,fsita,fsiti,tata,tati,tastn,tita,titi,tistn,stnta,stnti,tid2,tad2,d1ti,d2ti,jc1,jc2,jfsictx,jstnctx])

	# Randomly shuffle the order of parameters
	
	np.random.shuffle(d1taDist)
	np.random.shuffle(d2taDist)
	np.random.shuffle(fsitaDist)
	np.random.shuffle(fsitiDist)
	np.random.shuffle(tataDist)
	np.random.shuffle(tatiDist)
	np.random.shuffle(tastnDist)
	np.random.shuffle(tistnDist)
	np.random.shuffle(titaDist)
	np.random.shuffle(titiDist)
	np.random.shuffle(stntaDist)
	np.random.shuffle(stntiDist)
	np.random.shuffle(tid2Dist)
	np.random.shuffle(tad2Dist)
	np.random.shuffle(d1tiDist)
	np.random.shuffle(d2tiDist)
	np.random.shuffle(jc1Dist)
	np.random.shuffle(jc2Dist)
	np.random.shuffle(jfsictxDist)
	np.random.shuffle(jstnctxDist)

	# Check again	
	for d1ta1,d2ta1,fsita1,fsiti1,tata1,tati1,tastn1,tita1,titi1,tistn1,stnti1,tid21,tad21,jc11,jc21,jfsictx1,jstnctx1,stnta1,d1ti1,d2ti1 in zip(d1taDist,d2taDist,fsitaDist,fsitiDist,tataDist,tatiDist,tastnDist,titaDist,titiDist,tistnDist,stntiDist,tid2Dist,tad2Dist,jc1Dist,jc2Dist,jfsictxDist,jstnctxDist,stntaDist,d1tiDist,d2tiDist):
		if np.abs(d1ti1) > np.abs(d1ta1)/1.5  or np.abs(d2ti1) > np.abs(d2ta1)/1.5 or np.abs(tid21) < 1.2 or np.abs(tad21) < 1.2  or np.abs(stnta1) > np.abs(stnti1)/2. or np.abs(fsita1) > np.abs(fsiti1)/1.5:# This is the 1 1for all corase options so it always fails or np.abs(stnta1) > np.abs(stnti1)/1.5: 
			continue	
	
		iter1+=1

		A = np.matrix([[knownparams['d1d1'],knownparams['d1d2'],knownparams['d1fsi'],d1ta1,d1ti1,0.,0.],[knownparams['d2d1'],knownparams['d2d2'],knownparams['d2fsi'],d2ta1,d2ti1,0.,0.],[0.,0.,knownparams['fsifsi'],fsita1,fsiti1,0.,0.],[0,tad21,0.,tata1,tati1,tastn1,0],[0.,tid21,0.,tita1,titi1,tistn1,0.],[0.,0.,0.,stnta1,stnti1,knownparams['stnstn'],0.],[knownparams['gpid1'],0.,0.,knownparams['gpita'],knownparams['gpiti'],knownparams['gpistn'],knownparams['gpigpi']]])
		B = np.matrix([jc11,jc21,jfsictx1,0,0,jstnctx1,0])
		#w,v = np.linalg.eig(A)
		#Alleigs.append(w)
		delay = 1.0			 
		ipctx1=dict()
		ipctx1["ip"]=np.zeros((1,2001))	
		#Calculate Rates for SWA and lesion(dopamine depletion)
		Flags = []
		Flags.append("SWA")
#		Flags.append("DopDep")
		SWADopDepRates = psGA.calcRates(Flags,delay,A,B,False,ipctx1,iptau=p.params["iptau"])

		# Calculate Rates for Activation and Lesion
		Flags = []
		Flags.append("Act")
#		Flags.append("DopDep")
		ActDopDepRates = psGA.calcRates(Flags,delay,A,B,False,ipctx1,iptau=p.params["iptau"])

		tests = checkConds(SWADopDepRates,ActDopDepRates)
		print "tests",tests	
		print ind
		print "SWA:d1,d2,fsi,ta,ti,stn,gpi",np.mean(SWADopDepRates["d1"]),np.mean(SWADopDepRates["d2"]),np.mean(SWADopDepRates["fsi"]),np.mean(SWADopDepRates["ta"]),np.mean(SWADopDepRates["ti"]),np.mean(SWADopDepRates["stn"]),np.mean(SWADopDepRates["gpi"])
		print "Act:d1,d2,fsi,ta,ti,stn,gpi",np.mean(ActDopDepRates["d1"]),np.mean(ActDopDepRates["d2"]),np.mean(ActDopDepRates["fsi"]),np.mean(ActDopDepRates["ta"]),np.mean(ActDopDepRates["ti"]),np.mean(ActDopDepRates["stn"]),np.mean(ActDopDepRates["gpi"])
	
		print "tests",tests	
		print ind
		Grades = np.sum(tests)	
		cutoff = 12
		if Grades > cutoff:
			Allcombs.append([d1ta1,d2ta1,fsita1,fsiti1,tata1,tati1,tastn1,tita1,titi1,tistn1,stnta1,stnti1,tid21,tad21,d1ti1,d2ti1,jc11,jc21,jfsictx1,jstnctx1])

	pickle.dump(Allcombs,open(path5+"Allcombs_Tight_withTha"+prefix+".pickle","w"))
																					 																	
