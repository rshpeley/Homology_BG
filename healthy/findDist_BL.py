import pickle
import numpy as np
#import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.cm as cm
from pylab import *
import paramsearchGA_DopDep_nonlinear_BL as psGA
# Since its a sparse matrix
import scipy.sparse.linalg as sp
import itertools
import scipy.cluster.hierarchy as sph
import os
import knownUnknownParams as p
knownparams = p.params["known"]
unknownparams = p.params["unknown"] 
leak = -0.05
Allcombs =[]
Allstability = []
AllAs = []
storage_home = "/home/j.bahuguna/homology/Mallet/wojc1jc2/slowNormParams/newParams/strictHealthy"

path5 = storage_home+'/output/' # Path for output 
def checkConds(SWARates,ActRates):
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
	if np.round(np.mean(SWARates['ta'])) >= 0 and np.round(np.mean(SWARates['ta'])) <= 5 :# and SWARates['taFF'] < 1.5:
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
	if SWARates['taFF'] < 0.9 and SWARates['tiFF'] < 0.9:
		tests[4] = 1 # Refer to paramSearch_nonlinear.py in ../ directory
	
	# Check if D1 >= D2, since healthy state
	if np.round(np.mean(ActRates['ti'])) >= 12 and np.round(np.mean(ActRates['ti'])) <=50:		# Mean around 14
	#if np.mean(ActRates['ti']) > np.mean(ActRates['ta']):		# Mean around 14
	#if np.mean(ActRates['ti']) - np.mean(ActRates['ta']) >=1./1.:		# Mean around 14
		tests[5] = 1
	# Sanity test , all rates are fairly above zero
	if np.mean(SWARates['d1'][1000:]) > 0.1 and np.mean(SWARates['d2'][100./p.params["dt"]:]) > 0.1 and np.mean(SWARates['fsi'][1000:]) > 0.1 and 			np.mean(SWARates['ta'][1000:]) > 0.5 and np.mean(SWARates['ti'][1000:]) > 0.5 and np.mean(SWARates['stn'][1000:]) > 0.1 and np.mean(SWARates['gpi'][100./p.params["dt"]:]) > 0.1:
		tests[6] = 1
		#	print "d1",np.mean(SWARates['d1'][1000:])
		#	print "d2",np.mean(SWARates['d2'][1000:])
		#	print "fsi",np.mean(SWARates['fsi'][1000:])
		#	print "ta",np.mean(SWARates['ta'][1000:])
		#	print "ti",np.mean(SWARates['ti'][1000:])
		#	print "stn",np.mean(SWARates['stn'][1000:])
		#	print "gpi",np.mean(SWARates['gpi'][1000:])

	#if GpeSWA > np.mean(SWARates['d1']) and GpeSWA > np.mean(SWARates['d2']):
	if np.mean(ActRates["ti"]) >np.mean(SWARates["ti"]) and np.mean(ActRates["ta"]) > np.mean(SWARates["ta"]):
		tests[7] = 1

	#if abs(np.mean(SWARates['stn']) - np.mean(ActRates['stn']))<=10 and np.mean(SWARates['stn'])<=40: # Additonal constraint since STN firing rates were going very high. Also consult supplementary figure in Mallet 2008
        if np.mean(SWARates['stn'])<=50 and  np.mean(SWARates['d1']) /np.mean(SWARates['d2']) > 1.5 and np.mean(ActRates['d1']) /np.mean(ActRates['d2']) >1.5  : # Additonal constraint since STN firing rates were going very high. Also consult supplementary figure in Mallet 2008 (np.mean(ActRates['d1']) /np.mean(ActRates['d2'])) >  (np.mean(SWARates['d1']) /np.mean(SWARates['d2'])) and
      		tests[8] = 1	
	if np.round(np.mean(ActRates['ta'])) >= 5 and np.round(np.mean(ActRates['ta'])) <=25: #and np.mean(SWARates['fsi'])*1.5 < np.mean(SWARates['d1']): # Mean around 19
	#if np.mean(SWARates['ti']) > np.mean(SWARates['ta']): # Mean around 19
	#if np.mean(SWARates['ti']) - np.mean(SWARates['ta']) >=1./1.: # Mean around 19
		tests[9] = 1

	return tests




def run(prefix,pars):
	iter1 = 0
	AllcombswoConds =[]
	AllcombswithConds =[]

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
	#jc2 = unknownparams['jc1'][ind[16]] # jc2 = jc1 !! The disbalance is not due to striatum 
	jfsictx = unknownparams['jfsictx'][ind[18]]
	jstnctx = unknownparams['jstnctx'][ind[19]]

	num = 50
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

	'''
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
	'''

	for d1ta1,d2ta1,fsita1,fsiti1,tata1,tati1,tastn1,tita1,titi1,tistn1,stnti1,tid21,tad21,jc11,jc21,jfsictx1,jstnctx1,stnta1,d1ti1,d2ti1 in zip(d1taDist,d2taDist,fsitaDist,fsitiDist,tataDist,tatiDist,tastnDist,titaDist,titiDist,tistnDist,stntiDist,tid2Dist,tad2Dist,jc1Dist,jc2Dist,jfsictxDist,jstnctxDist,stntaDist,d1tiDist,d2tiDist):
		#if np.abs(d1ti1) > np.abs(d1ta1)/np.random.uniform(1.2,2.5,1)  or np.abs(d2ti1) > np.abs(d2ta1)/np.random.uniform(1.2,2.5,1) or np.abs(tid21) > 1 or np.abs(tad21) > 1 or np.abs(stnta1) > np.abs(stnti1)/np.random.uniform(2.,2.5,1) or np.abs(fsita1) > np.abs(fsiti1)/np.random.uniform(1.2,2.5,1) : #or np.abs(jc21) > np.abs(jc11) : # or  # Try to restrict tid2 /tad2 to low values
		#	continue	
		#if np.abs(jc11) > np.abs(jc21)*1.2:
		print "jc11,jc21",jc11,jc21
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
		#Flags.append("DopDep")
		SWARates = psGA.calcRates(Flags,delay,A,B,True,ipctx1,iptau=p.params["iptau"])

		# Calculate Rates for Activation and Lesion
		Flags = []
		Flags.append("Act")
		#	Flags.append("DopDep")
		ActRates = psGA.calcRates(Flags,delay,A,B,True,ipctx1,iptau=p.params["iptau"])
		tests = checkConds(SWARates,ActRates)
		print "tests",tests	
		print ind
		print "SWA:d1,d2,fsi,ta,ti,stn,gpi",np.mean(SWARates["d1"]),np.mean(SWARates["d2"]),np.mean(SWARates["fsi"]),np.mean(SWARates["ta"]),np.mean(SWARates["ti"]),np.mean(SWARates["stn"]),np.mean(SWARates["gpi"])
		print "Act:d1,d2,fsi,ta,ti,stn,gpi",np.mean(ActRates["d1"]),np.mean(ActRates["d2"]),np.mean(ActRates["fsi"]),np.mean(ActRates["ta"]),np.mean(ActRates["ti"]),np.mean(ActRates["stn"]),np.mean(ActRates["gpi"])
		Grades = np.sum(tests)	
		cutoff = 9
		if Grades > cutoff:
			if np.abs(d1ti1)*np.random.uniform(1.2,2.5,1) < np.abs(d1ta1)  or np.abs(d2ti1)* np.random.uniform(1.2,2.5,1)> np.abs(d2ta1)  or np.abs(tid21) < 1.0 or np.abs(tad21) < 1.0  or np.abs(stnta1)* np.random.uniform(2.,2.5,1) < np.abs(stnti1) or np.abs(fsita1)*np.random.uniform(1.2,2.5,1) < np.abs(fsiti1):
				AllcombswithConds.append([d1ta1,d2ta1,fsita1,fsiti1,tata1,tati1,tastn1,tita1,titi1,tistn1,stnta1,stnti1,tid21,tad21,d1ti1,d2ti1,jc11,jc21,jfsictx1,jstnctx1])
			else:
				AllcombswoConds.append([d1ta1,d2ta1,fsita1,fsiti1,tata1,tati1,tastn1,tita1,titi1,tistn1,stnta1,stnti1,tid21,tad21,d1ti1,d2ti1,jc11,jc21,jfsictx1,jstnctx1])


	# Have at leats one more combination, the original one
	#if np.abs(d1ti)*np.random.uniform(1.2,2.5,1) <  np.abs(d1ta)  and np.abs(d2ti)*np.random.uniform(1.2,2.5,1) < np.abs(d2ta)  and np.abs(tid2) < 1 and np.abs(tad2) < 1. and np.abs(stnta1)*np.random.uniform(2.,2.5,1) < np.abs(stnti1) and np.abs(fsita)*np.random.uniform(1.2,2.5,1) < np.abs(fsiti) and np.abs(d2ta) > np.abs(d1ta):# and np.abs(jc1) > np.abs(jc2)*1.2 :    # and np.abs(jc1) >1.0*np.abs(jc2) :# This is the 1 1for all corase options so it always fails or np.abs(stnta1) > np.abs(stnti1)/1.5: 
		# The original combination
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
	#Flags.append("DopDep")
	SWARates = psGA.calcRates(Flags,delay,A,B,True,ipctx1,iptau=p.params["iptau"])

	# Calculate Rates for Activation and Lesion
	Flags = []
	Flags.append("Act")
	ActRates = psGA.calcRates(Flags,delay,A,B,True,ipctx1,iptau=p.params["iptau"])
	tests = checkConds(SWARates,ActRates)
	print "tests",tests	
	print ind
	print "SWA:d1,d2,fsi,ta,ti,stn,gpi",np.mean(SWARates["d1"]),np.mean(SWARates["d2"]),np.mean(SWARates["fsi"]),np.mean(SWARates["ta"]),np.mean(SWARates["ti"]),np.mean(SWARates["stn"]),np.mean(SWARates["gpi"])
	print "Act:d1,d2,fsi,ta,ti,stn,gpi",np.mean(ActRates["d1"]),np.mean(ActRates["d2"]),np.mean(ActRates["fsi"]),np.mean(ActRates["ta"]),np.mean(ActRates["ti"]),np.mean(ActRates["stn"]),np.mean(ActRates["gpi"])

	Grades = np.sum(tests)	
	cutoff = 9
	if Grades > cutoff:
		if np.abs(d1ti)*np.random.uniform(1.2,2.5,1) <  np.abs(d1ta)  and np.abs(d2ti)*np.random.uniform(1.2,2.5,1) < np.abs(d2ta)  and np.abs(tid2) < 1.0 and np.abs(tad2) < 1.0 and np.abs(stnti) > np.abs(stnta)*np.random.uniform(2.,2.5,1) and np.abs(fsita)*np.random.uniform(1.2,2.5,1) < np.abs(fsiti):
			AllcombswithConds.append([d1ta,d2ta,fsita,fsiti,tata,tati,tastn,tita,titi,tistn,stnta,stnti,tid2,tad2,d1ti,d2ti,jc1,jc2,jfsictx,jstnctx])
		else:
			AllcombswoConds.append([d1ta,d2ta,fsita,fsiti,tata,tati,tastn,tita,titi,tistn,stnta,stnti,tid2,tad2,d1ti,d2ti,jc1,jc2,jfsictx,jstnctx])


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

	for d1ta1,d2ta1,fsita1,fsiti1,tata1,tati1,tastn1,tita1,titi1,tistn1,stnti1,tid21,tad21,jc11,jc21,jfsictx1,jstnctx1,stnta1,d1ti1,d2ti1 in zip(d1taDist,d2taDist,fsitaDist,fsitiDist,tataDist,tatiDist,tastnDist,titaDist,titiDist,tistnDist,stntiDist,tid2Dist,tad2Dist,jc1Dist,jc2Dist,jfsictxDist,jstnctxDist,stntaDist,d1tiDist,d2tiDist):
		#if np.abs(d1ti1) > np.abs(d1ta1)/np.random.uniform(1.2,2.5,1)  or np.abs(d2ti1) > np.abs(d2ta1)/np.random.uniform(1.2,2.5,1) or np.abs(tid21) > 1 or np.abs(tad21) > 1. or np.abs(stnta1) > np.abs(stnti1)/np.random.uniform(1.2,2.5,1) or np.abs(fsita1) > np.abs(fsiti1)/np.random.uniform(1.2,2.5,1) or np.abs(d1ta1) > np.abs(d2ta1):# or np.abs(jc21) > np.abs(jc11) : #np.abs(stnta1) > np.abs(stnti1)/1.2 or  # Try to restrict tid2 /tad2 to low values
		#	continue	
		#if np.abs(jc11) > np.abs(jc21)*1.2:
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
		#	Flags.append("DopDep")
		SWARates = psGA.calcRates(Flags,delay,A,B,True,ipctx1,iptau=p.params["iptau"])

		# Calculate Rates for Activation and Lesion
		Flags = []
		Flags.append("Act")
		#Flags.append("DopDep")
		ActRates = psGA.calcRates(Flags,delay,A,B,True,ipctx1,iptau=p.params["iptau"])
		tests = checkConds(SWARates,ActRates)
		print "tests",tests	
		print ind
		print "SWA:d1,d2,fsi,ta,ti,stn,gpi",np.mean(SWARates["d1"]),np.mean(SWARates["d2"]),np.mean(SWARates["fsi"]),np.mean(SWARates["ta"]),np.mean(SWARates["ti"]),np.mean(SWARates["stn"]),np.mean(SWARates["gpi"])
		print "Act:d1,d2,fsi,ta,ti,stn,gpi",np.mean(ActRates["d1"]),np.mean(ActRates["d2"]),np.mean(ActRates["fsi"]),np.mean(ActRates["ta"]),np.mean(ActRates["ti"]),np.mean(ActRates["stn"]),np.mean(ActRates["gpi"])

		Grades = np.sum(tests)	
		cutoff = 9
		if Grades > cutoff:
			if np.abs(d1ti1)*np.random.uniform(1.2,2.5,1) < np.abs(d1ta1)  or np.abs(d2ti1)* np.random.uniform(1.2,2.5,1)> np.abs(d2ta1)  or np.abs(tid21) < 1.0 or np.abs(tad21) < 1.0  or np.abs(stnta1)* np.random.uniform(2.,2.5,1) < np.abs(stnti1) or np.abs(fsita1)*np.random.uniform(1.2,2.5,1) < np.abs(fsiti1):
				AllcombswithConds.append([d1ta1,d2ta1,fsita1,fsiti1,tata1,tati1,tastn1,tita1,titi1,tistn1,stnta1,stnti1,tid21,tad21,d1ti1,d2ti1,jc11,jc21,jfsictx1,jstnctx1])
			else:
				AllcombswoConds.append([d1ta1,d2ta1,fsita1,fsiti1,tata1,tati1,tastn1,tita1,titi1,tistn1,stnta1,stnti1,tid21,tad21,d1ti1,d2ti1,jc11,jc21,jfsictx1,jstnctx1])



	pickle.dump(AllcombswithConds,open(path5+"Allcombs_withConds_BL"+prefix+".pickle","w"))
	pickle.dump(AllcombswoConds,open(path5+"Allcombs_woConds_BL"+prefix+".pickle","w"))

																					 																	
